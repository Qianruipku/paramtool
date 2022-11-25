#include "function.h"
#include <iostream>
#include <cmath>
#include "const.h"
#include "tool.h"

double fd_integral(const double x, const double j);

void FUNC::conductivity()
{
    vector<double> mlist, zlist, nlist, denlist_i;
	read_elemets(mlist, zlist, nlist);
    double rho_i = read_density(); //g/cm3
	double T_eV = read_temperature();
    //--------------------------------------------------------
    double z_per_mol(0), mass_per_mol(0), n_per_mol(0);
	for(int i = 0 ; i < nlist.size(); ++i)
	{
        mass_per_mol += nlist[i] * mlist[i];
        z_per_mol += nlist[i] * zlist[i];
		n_per_mol += nlist[i];
	}
    double z_avg = z_per_mol / n_per_mol;
    double molden = rho_i / mass_per_mol * P_NA * n_per_mol; //unit: cm^-3
    for(int i = 0 ; i < nlist.size(); ++i)
	{
        denlist_i.push_back(molden * nlist[i] / n_per_mol); //unit: cm^-3
    }
    double density_e = molden * z_avg; //unit cm^-3
    double mu_eV = FEG_mu(density_e, T_eV);
    cout<<"density: "<<density_e<<yellow(" g/cm^3")<<" ; temperature: "<<T_eV<<yellow(" eV")<<endl;
    cout<<"Fermi energy: "<<mu_eV<<yellow(" eV")<<" Tf/T = "<<mu_eV/T_eV<<endl;

    //--------------------------------------------------------
    lee_more(T_eV, mu_eV, density_e, denlist_i, zlist);
    
}

//density_e: cm^-3; denlist_i: cm^-3
void FUNC:: lee_more(const double T_eV, const double mu_eV, const double density_e, const vector<double>& denlist_i, const vector<double>& zlist)
{
    // Hartree atomic unit
    double kT = T_eV / Ha2eV;
    double kTf = mu_eV / Ha2eV;
    double mu_kT = mu_eV / T_eV;
    double density_e_au = density_e*pow(P_bohr*1e-8, 3); 
    //----------------------------------
    double alpha, beta;
    double F1_2;
    double expmultiF12 = -1; //(1+exp(-mu_kT))*F1_2
    if(mu_kT < -20)
    {
        alpha = 32.0 /3.0 / M_PI;
        beta = 128.0 / 3.0 / M_PI;
        F1_2 = sqrt(M_PI) / 2.0 * exp(mu_kT);
        expmultiF12 = 1.0;
    }
    else if (mu_kT > 1e4)
    {
        alpha = 1.0;
        beta = M_PI * M_PI / 3.0;
        F1_2 = pow(mu_kT, 1.5)/1.5;
        expmultiF12 = F1_2;
    }
    else
    {
        F1_2 = fd_integral(mu_kT, 1.0/2.0);
        double F2 = fd_integral(mu_kT, 2.0);
        double F3 = fd_integral(mu_kT, 3.0);
        double F4 = fd_integral(mu_kT, 4.0);
        alpha = 4.0/3.0 * F2/((1 + exp(-mu_kT)) * pow(F1_2,2));
        beta = 20.0/9.0 * F4 * (1 - 16.0*pow(F3, 2)/(15.0*F4*F2))/((1+exp(-mu_kT))*pow(F1_2,2));
        expmultiF12 = (1+exp(-mu_kT))*F1_2;
        // cout<<F1_2<<" "<<F2<<" "<<F3<<" "<<F4<<endl;
        // cout<<alpha<<" "<<beta<<endl;
    }
    //----------------------------------
    double bmax,bmin;
    double Debye = 4*M_PI* density_e_au/ sqrt(pow(kT,2) + pow(kTf,2));
    bmin = M_PI/sqrt(3*kT);
    double tau_frac = 0;
    double density_i_tot_au = 0;
    for(int i = 0 ; i < denlist_i.size(); ++i)
    {
        double density_i_au = denlist_i[i] * pow(P_bohr*1e-8, 3);
        Debye += 4*M_PI* density_i_au *pow(zlist[i] ,2) / kT;
        bmin = max(zlist[i]/(3*kT), bmin);
        tau_frac += pow(zlist[i],2) * density_i_au;
        density_i_tot_au += density_i_au;
    }
    bmax = max(sqrt(1.0/Debye), pow(3.0/(4.0*M_PI*density_i_tot_au), 1.0/3.0));
    cout<<bmax<<" "<<bmin<<" "<<bmax/bmin<<endl;
    double cou_log = 1.0/2.0 * log(1.0 + pow(bmax/bmin, 2));
    cou_log = max(cou_log,2.0);
    double tau = 1.0/tau_frac * 3.0 * pow(kT,3.0/2.0)/(2.0*sqrt(2.0)*M_PI*cou_log)*expmultiF12;
    //----------------------------------
    double sigma_au = density_e_au * tau * alpha;
    double kappa_au = density_e_au * kT * tau * beta;
    //----------------------------------
    const double au2si_sigma = hau2A * hau2A * hau2s / hau2J / hau2m;
    const double au2si_kappa = hau2J /(hau2s * hau2m * hau2K);
    double sigma = sigma_au * au2si_sigma;
    double kappa = kappa_au * au2si_kappa;
    cout<<"Lee-More:"<<endl;
    cout<<"electrical conductivity: "<<sigma<<yellow(" S m^-1")<<endl;
    cout<<"thermal conductivity: "<<kappa<<yellow(" W (mK)^-1")<<endl;
}

double fd_integral(const double x, const double j)
{
    double de, thr;
    if(x <= 30)
    {
        de = 1e-4;
        thr = 1e-18;
    }
    else
    {
        de = 1e-2;
        thr = 1e-8;
    }
    double last_term = 1 + thr;
    auto func = [j,x](const double t) -> double{
        return pow(t,j) / (1+exp(t-x));
    };
    double e = 0;
	double sum = func(e);
    //when t > 2*j && t > x
    //dfunc/dt < 0
	while(last_term > thr || e <= 2*j || e <= x)
	{
		e += de;
		sum += 4 * func(e);
		e += de;
		last_term = func(e);
		sum += 2 * last_term;
	}
	sum += func(e+de);
	sum /= 3;
	return sum * de;
}
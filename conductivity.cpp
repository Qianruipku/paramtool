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
    double F1_2 = fd_integral(mu_kT, 1.0/2.0);
    double F2 = fd_integral(mu_kT, 2.0);
    double F3 = fd_integral(mu_kT, 3.0);
    double F4 = fd_integral(mu_kT, 4.0);
    double alpha = 4.0/3.0 * F3/((1 + exp(-mu_kT)) * pow(F1_2,2));
    double beta = 20.0/9.0 * F4 * (1 - 16.0*pow(F3, 2)/(15.0*F4*F2))/((1+exp(-mu_kT))*pow(F1_2,2));
    //----------------------------------
    double bmax,bmin;
    double Debye = 4*M_PI* density_e_au/ sqrt(pow(kT,2) + pow(kTf,2));
    bmin = M_PI/sqrt(3*kT);
    double tau_frac = 0;
    for(int i = 0 ; i < denlist_i.size(); ++i)
    {
        double density_i_au = denlist_i[i] * pow(P_bohr*1e-8, 3);
        Debye += 4*M_PI* density_i_au *pow(zlist[i] ,2) / kT;
        bmin = max(zlist[i]/(3*kT), bmin);
        tau_frac += 1.0/(pow(zlist[i],2) * density_i_au);
    }
    bmax = sqrt(1.0/Debye);
    double cou_log = 1.0/2.0 * log(1.0 + pow(bmax/bmin, 2));
    double tau = tau_frac * 3.0 * pow(kT,3.0/2.0)/(2.0*sqrt(2.0)*M_PI*cou_log)*(1+exp(-mu_kT))*F1_2;
    //----------------------------------
    double sigma_au = density_e * tau * alpha;
    double kappa_au = density_e * kT * tau * beta;
    //----------------------------------
    const double au2si_sigma = hau2A * hau2A * hau2s / hau2J / hau2m;
    const double au2si_kappa = hau2J /(hau2s * hau2m * hau2K);
    double sigma = sigma * au2si_sigma;
    double kappa = kappa_au * au2si_kappa;
    cout<<"Lee-More:"<<endl;
    cout<<"electrical conductivity: "<<sigma<<yellow(" S m^-1")<<endl;
    cout<<"thermal conductivity: "<<kappa<<yellow(" W (mK)^-1")<<endl;
}

double fd_integral(const double x, const double j)
{
    double de = 1e-2;
    double last_term = 1;
    double thr = 1e-4;
    auto func = [j,x](const double t) -> double{
        return pow(t,j) / (1+exp(t-x));
    };
	double sum = func(0);
	double e = 0;
	int i = 0;
	while(last_term > thr)
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
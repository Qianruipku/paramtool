#include "function.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "const.h"
#include "tool.h"

double fd_integral(const double x, const double j);
double getA_alpha(const double mu_kT);
double getA_beta(const double mu_kT);

void FUNC::conductivity()
{
    vector<double> mlist, zlist, nlist, denlist_i, zionlist;
	read_elemets(mlist, zlist, nlist, zionlist);
    double rho_i = read_density(); //g/cm3
	double T_eV = read_temperature();

    thomas_fermi_ionization(rho_i, T_eV, mlist, zlist, nlist, zionlist);
    //--------------------------------------------------------
    molecule mol0(mlist, zlist, nlist);
    molecule mol(mlist, zionlist, nlist);
    double ionization = mol.tot_z/mol0.tot_z;
    double den_mole = rho_i / (mol.avg_m/P_NA); //unit: cm^-3
    for(int i = 0 ; i < nlist.size(); ++i)
	{
        denlist_i.push_back(den_mole * nlist[i] / mol.tot_n); //unit: cm^-3
    }
    double density_e = den_mole * mol.avg_z; //unit cm^-3
    double mu_eV = FEG_mu(density_e, T_eV);
    cout<<"density: "<<density_e<<" "<<yellow("cm^-3")<<" ; temperature: "<<T_eV<<" "<<yellow("eV")<<endl;
    cout<<"Chemical potential: "<<mu_eV<<" "<<yellow("eV")<<" ; mu/T = "<<mu_eV/T_eV<<endl;
    cout<<"Ionization: "<<ionization*100<<"%"<<endl;
    // double  t(10), g(0.05);
    // cout<<pow(2/(3*M_PI*g*t),3)*4*(mol.avg_m/P_NA)/pow(P_bohr*1e-8, 3)<<" dd "<<2.0/(g*g*t*pow(9.0*M_PI/4, 2.0/3.0))*Ha2eV<<endl;
    cout<<"Coulomb-coupling parameter: "<<coupling_parameter(mol, T_eV, density_e)<<" ; Fermi-degeneracy parameter: "<<degeneracy_parameter(T_eV, density_e)<<endl;
    //--------------------------------------------------------
    cout<<std::left<<setw(20)<<"Model"<<setw(29)<<"Sigma"+yellow("(Sm^-1)")<<setw(29)<<"Kappa"+yellow("(W(mK)^-1)")<<setw(20)<<"Lorentz number"<<endl;
    spitzer (T_eV, mu_eV, density_e, denlist_i, zionlist);
    lee_more(T_eV, mu_eV, density_e, denlist_i, zionlist);
    // Ichimaru(T_eV, mu_eV, density_e, denlist_i, zlist);
}

void FUNC:: spitzer(const double T_eV, const double mu_eV, const double density_e, 
                        const vector<double>& denlist_i, const vector<double>& zionlist)
{
    // Hartree atomic unit
    double kT_au = T_eV / Ha2eV;
    double kTf_au = fermi_energy(density_e) / Ha2eV;
    double mu_kT = mu_eV / T_eV;
    double density_e_au = density_e*pow(P_bohr*1e-8, 3); 
    //---------mean ionic charge---------
    double Z_avg = 0;
    for(int i = 0 ; i < denlist_i.size(); ++i)
    {
        Z_avg += denlist_i[i] * zionlist[i] * zionlist[i];
    }
    Z_avg /= density_e;
    //----------Coulomb logarithm--------
    double Lambda = 3.0*pow(kT_au, 1.5) / (2*sqrt(M_PI*density_e_au)*Z_avg);
    double T_K = T_eV * eV2K;
    if(T_K > 4.2e5)
    {
        Lambda *= sqrt(4.2e5/T_K);
    }
    if(Lambda < 1 || density_e > 5e5*pow(T_K,3))
    {
        cout<<std::left<<setw(20)<<"Spitzer"<<setw(20)<<"---"<<setw(20)<<"---"<<setw(20)<<"---"<<endl;
        return;
    }
    double cou_log = log(Lambda);
    //-------------conductivity----------
    double sigma_au = 4.0*sqrt(2*M_PI)*pow(kT_au, 1.5) / (Z_avg*pow(M_PI, 2)*cou_log);
    double kappa_au = 40.0*sqrt(2*M_PI)*pow(kT_au, 2.5) / (Z_avg*pow(M_PI, 2)*cou_log);
    std::vector<double> ref_invZ{ 0.0, 1.0/16.0, 1.0/4.0, 1.0/2.0, 1.0 };
    std::vector<double> ref_gamma{1.0, 0.923, 0.785, 0.683, 0.582};
    std::vector<double> ref_delta{1.0, 0.791, 0.513, 0.356, 0.225};
    std::vector<double> ref_epsilon{0.4, 0.396, 0.401, 0.410, 0.419};
    std::vector<double> invZ(1), gamma_E(1), delta_T(1), epsilon(1);
    invZ[0] = 1.0 / Z_avg;
    NaturalSplineInterpolation(invZ, gamma_E, ref_invZ, ref_gamma);
    NaturalSplineInterpolation(invZ, delta_T, ref_invZ, ref_delta);
    NaturalSplineInterpolation(invZ, epsilon, ref_invZ, ref_epsilon);
    sigma_au *= gamma_E[0];
    kappa_au *= delta_T[0] * epsilon[0];
    //----------------------------------
    const double au2si_sigma = hau2A * hau2A * hau2s / hau2J / hau2m;
    const double au2si_kappa = hau2J /(hau2s * hau2m * hau2K);
    double sigma = sigma_au * au2si_sigma;
    double kappa = kappa_au * au2si_kappa;
    cout<<std::left<<setw(20)<<"Spitzer"<<setw(20)<<sigma<<setw(20)<<kappa<<setw(20)<<kappa_au/sigma_au/kT_au<<endl;
}

void FUNC:: lee_more(const double T_eV, const double mu_eV, const double density_e, const vector<double>& denlist_i, const vector<double>& zionlist)
{
    // Hartree atomic unit
    double kT_au = T_eV / Ha2eV;
    double kTf_au = fermi_energy(density_e) / Ha2eV;
    double mu_kT = mu_eV / T_eV;
    double density_e_au = density_e*pow(P_bohr*1e-8, 3); 
    //---------mean ionic charge---------
    double Z_avg = 0;
    for(int i = 0 ; i < denlist_i.size(); ++i)
    {
        Z_avg += denlist_i[i] * zionlist[i] * zionlist[i];
    }
    Z_avg /= density_e;
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
    //----Lee-more approximation--------
    // alpha = getA_alpha(mu_kT);
    // beta = getA_beta(mu_kT);
    //----------------------------------
    double bmax,bmin;
    double Debye = 4*M_PI* density_e_au/ sqrt(pow(kT_au,2) + pow(kTf_au,2));
    bmin = M_PI/sqrt(3*kT_au);
    double tau_frac = 0;
    double density_i_tot_au = 0;
    for(int i = 0 ; i < denlist_i.size(); ++i)
    {
        double density_i_au = denlist_i[i] * pow(P_bohr*1e-8, 3);
        Debye += 4*M_PI* density_i_au *pow(zionlist[i] ,2) / kT_au;
        tau_frac += pow(zionlist[i],2) * density_i_au;
        density_i_tot_au += density_i_au;
    }
    if(kT_au < kTf_au) //degenrate limit, bim = Ze^2/2E_F, E_F=3/2 kT_F
    { 
        // correction according to the description in the Lee-More paper.
        bmin = Z_avg / sqrt(pow(3*kT_au*exp(10*(kT_au/kTf_au-1)), 2) + pow(3*kTf_au, 2));
        // bmin = Z_avg / 3*kTf_au;
    }
    else
    {
        bmin = max(Z_avg/(3*kT_au), bmin);
    }
    bmax = max(sqrt(1.0/Debye), 2*pow(3.0/(4.0*M_PI*density_i_tot_au), 1.0/3.0));
    double cou_log = 1.0/2.0 * log(1.0 + pow(bmax/bmin, 2));
    cou_log = max(cou_log,2.0);
    double tau = 1.0/tau_frac * 3.0 * pow(kT_au,3.0/2.0)/(2.0*sqrt(2.0)*M_PI*cou_log)*expmultiF12;
    //----------------------------------
    double sigma_au = density_e_au * tau * alpha;
    double kappa_au = density_e_au * kT_au * tau * beta;
    //----------------------------------
    const double au2si_sigma = hau2A * hau2A * hau2s / hau2J / hau2m;
    const double au2si_kappa = hau2J /(hau2s * hau2m * hau2K);
    double sigma = sigma_au * au2si_sigma;
    double kappa = kappa_au * au2si_kappa;
    cout<<std::left<<setw(20)<<"Lee-More"<<setw(20)<<sigma<<setw(20)<<kappa<<setw(20)<<kappa_au/sigma_au/kT_au<<endl;
}

void FUNC::Ichimaru(const double T_eV, const double mu_eV, const double density_e, 
                        const vector<double>& denlist_i, const vector<double>& zionlist)
{
    // Hartree atomic unit
    double T = T_eV / Ha2eV;
    double Ef = fermi_energy(density_e) / Ha2eV;
    double mu_kT = mu_eV / T_eV;
    double density_e_au = density_e*pow(P_bohr*1e-8, 3); 
    //----------------------------------
    const double gamma = 0.5772156649;
    const double expgamma = exp(gamma);
    double rs = pow(3.0/4/M_PI/density_e_au, 1.0/3.0);
    double xb = sqrt(rs * tanh(sqrt(2*M_PI/T) * pow(density_e_au, 1.0/3.0)));
    double Gamma_e = 1.0/rs/T;
    double theta = T / Ef;
    double fz_E(0), fz_T(0);
    double density_i_tot_au = 0;
    for(int i = 0 ; i < zionlist.size() ; ++i)
    {
        double Z = zionlist[i];
        if(Z > 26) 
        {
            cout<<"Ichimaru model do not support Z > 26."<<endl;
            return;
        }
        double zeta_DH = pow(Z+1, 1+1.0/Z) * expgamma * Gamma_e / pow(12*M_PI*M_PI, 1.0/3.0) / theta;
        double zeta_Born = 1.0/(2.5 * pow(theta, 1.5) * pow(Z,4.0/3.0)) * exp(-1.47*pow(Z,1.0/3.0));
        double LE = 0.5 * log(1 + 1/zeta_DH + tanh(1/zeta_Born))
                * (1 + 0.42*xb*xb*exp(-6e-4*rs*rs) + 0.063*pow(xb*xb*exp(-6e-4*rs*rs),5));
        double LT = 0.5 * log(1 + 75/(13*M_PI*M_PI)*(1/zeta_DH + tanh(1/zeta_Born)))
                * (1 + 0.38*xb*xb*exp(-6e-4*rs*rs) + 0.049*pow(xb*xb*exp(-6e-4*rs*rs),5));
        fz_E += Z*Z*LE;
        fz_T += Z*Z*LT;
        density_i_tot_au += denlist_i[i] * pow(P_bohr*1e-8, 3);

    }
    double rho_E = 8.0/3.0*sqrt(M_PI/2) * fz_E * density_i_tot_au / density_e_au / pow(T, 1.5);
    double rho_T = 52.0*sqrt(2*M_PI)/75 * fz_T * density_i_tot_au / density_e_au / pow(T, 2.5);
    //----------------------------------
    double sigma_au = 1 / rho_E;
    double kappa_au = 1 / rho_T;
    //----------------------------------
    const double au2si_sigma = hau2A * hau2A * hau2s / hau2J / hau2m;
    const double au2si_kappa = hau2J /(hau2s * hau2m * hau2K);
    double sigma = sigma_au * au2si_sigma;
    double kappa = kappa_au * au2si_kappa;
    cout<<"Ichimaru:"<<endl;
    cout<<"electrical conductivity: "<<sigma<<" "<<yellow("Sm^-1")<<endl;
    cout<<"thermal conductivity: "<<kappa<<" "<<yellow("W(mK)^-1")<<endl;
    cout<<"Lorenz number: "<<kappa_au/sigma_au/T<<endl;

}

double FUNC:: degeneracy_parameter(const double T_eV, const double density_e)
{
    return T_eV/fermi_energy(density_e);
}

double FUNC:: coupling_parameter(molecule &mol, const double T_eV, const double density_e)
{
    double a = pow(mol.avg_z/density_e*3.0/4/M_PI,1.0/3)*1e8/P_bohr; //a.u.
    return mol.avg_z * mol.avg_z / (a * T_eV / Ha2eV);
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

double getA_alpha(const double mu_kT)
{
    double a1 = 3.39;
    double a2 = 0.347;
    double a3 = 0.129;
    double b2 = 0.511;
    double b3 = 0.124;
    double y = log(1+exp(mu_kT));
    
    return (a1 + a2 * y + a3 * y * y)/(1 + b2*y + b3 * y * y);
}

double getA_beta(const double mu_kT)
{
    double a1 = 13.5;
    double a2 = 0.976;
    double a3 = 0.437;
    double b2 = 0.510;
    double b3 = 0.126;
    double y = log(1+exp(mu_kT));
    
    return (a1 + a2 * y + a3 * y * y)/(1 + b2*y + b3 * y * y);
}
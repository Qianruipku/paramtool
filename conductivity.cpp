#include "function.h"
#include <iostream>
#include <cmath>
#include "const.h"
#include "tool.h"

double fd_integral(const double x, const double j);

void FUNC::conductivity()
{
    vector<double> mlist, qlist, nlist;
	read_elemets(mlist, qlist, nlist);
    double density = read_density();
	double temperature = read_temperature();
    
}

void lee_more(const double T_eV, const double mu_eV, const double density_e, const double density_i, const double z_val)
{
    double kT = T_eV;
    double kTf = mu_eV;
    double mu_kT = mu_eV;
    // double density_e , density_i;
    // double z_val;
    //----------------------------------
    double F1_2 = fd_integral(mu_kT, 1.0/2.0);
    double F2 = fd_integral(mu_kT, 2.0);
    double F3 = fd_integral(mu_kT, 3.0);
    double F4 = fd_integral(mu_kT, 4.0);
    double alpha = 4.0/3.0 * F3/((1 + exp(-mu_kT)) * pow(F1_2,2));
    double beta = 20.0/9.0 * F4 * (1 - 16.0*pow(F3, 2)/(15.0*F4*F2))/((1+exp(-mu_kT))*pow(F1_2,2));
    //----------------------------------
    double bmax,bmin;
    double Debye = 4*M_PI*density_e * P_qe * P_qe / sqrt(pow(kT,2) + pow(kTf,2)) + 4*M_PI*density_i*pow(P_qe*z_val ,2) / kT;
    bmax = sqrt(1.0/Debye);
    bmin = max(z_val * P_qe*P_qe/(3*kT), P_hbar*2*M_PI/2/sqrt(3*P_Me*kT));
    double cou_log = 1.0/2.0 * log(1.0 + pow(bmax/bmin, 2));
    double tau = 3.0*sqrt(P_Me)*pow(kT,3.0/2.0)/(2.0*sqrt(2.0)*M_PI*pow(z_val,2)*density_i*pow(P_qe,4)*cou_log)*(1+exp(-mu_kT))*F1_2;
    //----------------------------------
    double sigma = density_e*pow(P_qe,2)*tau/P_Me * alpha * mu_kT;
    double kappa = density_e*P_kB*kT*tau/P_Me * beta * mu_kT;
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
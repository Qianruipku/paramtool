#include "function.h"
#include "const.h"
#include "elements.h"
#include "tool.h"
#include "iomanip"
#include <iostream>
using namespace std;
#define NT 60 
double getmu(double mu_0, double T);
double calint(double fun(double e, double mu, double T), double mu, double T,double thr);
double funmu(double e, double mu, double T);
double getEcut(double mu, double T);
void FUNC::mdparameter()
{
	double na, nmol;
	double density;
	double z_per_mol;
	string name;
	double mass_per_mol, n_per_mol, min_mass;
	double temp;
	cout<<green("Input the name of molecule (e.g. Mg-O, C-O-2, H-2-O)")<<endl;
	cin>>name;
	if(!readstring(name, mass_per_mol, z_per_mol, n_per_mol, min_mass))	exit(0);
	// cout<<name<<" "<<mass_per_mol<<" "<<z_per_mol<<" "<<n_per_mol<<" "<<min_mass<<endl;
	cout<<green("Input the number of molecules:")<<endl;
	cin>>nmol;
	na = n_per_mol * nmol;
	cout<<green("Input the density (g/cm^3):")<<endl;
    cin>>density;
	cout<<green("Input the temperature (eV):")<<endl;
    cin>>temp;

	double mass_atom_g = mass_per_mol/NA/n_per_mol;
	double l_ang=pow(na*mass_atom_g/density,1.00/3)*1e8;
	double l_bohr=l_ang/bohr;
	double WSr_ang=pow(mass_atom_g/density*3/4/M_PI,1.00/3)*1e8;
	double WSr_bohr=WSr_ang/bohr;
	double dt_fs=WSr_ang*1e-10/NT/sqrt(temp*qe/(0.001*min_mass))*1e15;
	double dt_au=dt_fs/(au2s*1e15);
	//mu of free electrons, assume all electrons are ionized.
	//mu0=hbar^2/2m*(3pi^2N/V)^(2/3)
	//mu = mu0(1 - pi^2/12*(kT/mu0)^2) (mu/T >> 1)
	double mu0_Ry = pow(3.0 * pow(M_PI,2) * nmol*z_per_mol/pow(l_bohr,3), 2.0/3.0); //unit in Ry
	double mu_eV,mu_Ry;
	if(temp < 1)
	{
		mu_Ry = mu0_Ry * (1 - pow(M_PI,2)/12.0*pow(temp/Ry2eV/mu0_Ry,2)-pow(M_PI,4)*7.0/960*pow(temp/Ry2eV/mu0_Ry,4));
		mu_eV = mu_Ry * Ry2eV;
	}
	else if(temp > 1e4)
	{
		mu_eV = temp*(-log(6*pow(M_PI,2)) + 1.5*log(4*M_PI*mu0_Ry*Ry2eV/temp));
		mu_Ry = mu_eV / Ry2eV;
	}
	else
	{
		mu_eV = getmu(mu0_Ry * Ry2eV, temp);
		mu_Ry = mu_eV / Ry2eV;
	}
	double Ecut1_eV = getEcut( mu_eV, temp);
	double Ecut1_Ry = Ecut1_eV / Ry2eV;
	double Ecut2_eV = mu_eV + 5*log(10.0)*temp;
	double Ecut2_Ry = Ecut2_eV / Ry2eV;
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Lattice constant: "<<fixed<<setprecision(8)<<l_bohr<<yellow(" Bohr; ")<<l_ang<<yellow(" Angstrom")<<endl;
	cout<<"Wigner-Seitz radius: "<<WSr_bohr<<yellow(" Bohr;")<<"  (0.7*WS = "<<0.7*WSr_bohr<<")"<<endl;
	cout<<"Temperature: "<<temp/Ry2eV<<yellow(" Ry; ")<<temp*eV2K<<yellow(" K")<<endl;
	cout<<"dt: "<<dt_fs<<yellow(" fs; ")<<dt_au<<yellow(" a.u.")<<endl;
	cout<<"1/(40dt): "<<double(1.0)/40/dt_fs<<yellow(" fs^-1")<<endl;
	cout<<"FEG mu: "<<setprecision(3)<<mu_eV<<yellow(" eV; ")<<mu_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut1 (interror < 1e-3): "<<Ecut1_eV<<yellow(" eV; ")<<Ecut1_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut2 (f<1e-5):          "<<Ecut2_eV<<yellow(" eV; ")<<Ecut2_Ry<<yellow(" Ry")<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	return;
}

//4/3 * mu_0^(3/2) = \int_0^\infty \frac{e^(1/2)}{exp((e-mu)/kT)+1} de
double getmu(double mu0, double T)
{
	double ref = 2.0/3.0 * pow(mu0, 3.0/2.0);
	double Deltamu = 5 * T;
	double mu1 = -Deltamu;
	double mu2 = Deltamu;
	double com1 = calint(funmu, mu1, T, 1e-8);
	double com2 = calint(funmu, mu2, T, 1e-8);
	while (com1 > ref)
    {
        mu2 = mu1;
        mu1 -= Deltamu;
        com1 = calint(funmu, mu1, T, 1e-8);
        Deltamu *= 10;
    }
    while (com2 < ref)
    {
        mu1 = mu2;
        mu2 += Deltamu;
        com2 = calint(funmu, mu2, T, 1e-8);
        Deltamu *= 10;
    }
	double diff = 1000;
	double mu3, com3;
	double thr = 1e-8 * ref;
	while (diff > thr)
    {
        mu3 = (mu2 + mu1) / 2;
        com3 = calint(funmu, mu3, T, 1e-8);
        if (com3 < ref)
        {
            mu1 = mu3;
        }
        else if (com3 > ref)
        {
            mu2 = mu3;
        }
        diff = abs(ref - com3);
    }
	return (mu2 + mu1) / 2;
}

double calint(double fun(double e, double mu, double T), double mu, double T, double thr)
{
	double de = T * 1e-4;
	double result = 1;
	double sum = fun(0, mu, T);
	double e = 0;
	int i = 0;
	while(result > thr)
	{
		e += de;
		sum += 4 * fun(e, mu, T);
		e += de;
		result = fun(e, mu, T);
		sum += 2 * result;
	}
	sum += fun(e+de, mu, T);
	sum /= 3;
	return sum * de;
}

double funmu(double e, double mu, double T)
{
	return sqrt(e)/(exp((e-mu)/T)+1);
}

double getEcut(double mu, double T)
{
	double ref = calint(funmu, mu, T, 1e-20);
	double de = T * 1e-4;
	double result = 1;
	double sum = funmu(0, mu, T);
	double e = 0;
	double diff = 1;
	while(diff > 1e-3)
	{
		e += de;
		sum += 4 * funmu(e, mu, T);
		e += de;
		sum += 2 * funmu(e,mu,T);
		diff = (ref - sum/3*de) / ref;
	}
	return e;
}
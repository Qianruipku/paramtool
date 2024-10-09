#include "function.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "const.h"
#include "tool.h"
double fd_integral(const double x, const double j);

void FUNC::Cv()
{
    vector<double> mlist, zlist, nlist, denlist_i, zionlist;
	read_elemets(mlist, zlist, nlist, zionlist);
    double rho_i = read_density(); //g/cm3
	double T_eV = read_temperature();

    thomas_fermi_ionization(rho_i, T_eV, mlist, zlist, nlist, zionlist);
    molecule mol(mlist, zionlist, nlist);
    double den_mole = rho_i / (mol.avg_m/P_NA); //unit: cm^-3
    double density_e = den_mole * mol.avg_z; //unit cm^-3
    //--------------------------------------------------------
    double T1_eV = T_eV * 0.999;
    double T2_eV = T_eV * 1.001;
    double mu1_eV = FEG_mu(density_e, T1_eV);
    double mu2_eV = FEG_mu(density_e, T2_eV);
    double E1_Ha = FEG_E(density_e, den_mole, T1_eV, mu1_eV);
    double E2_Ha = FEG_E(density_e, den_mole, T2_eV, mu2_eV);
    double Cv = (E2_Ha - E1_Ha) / (T2_eV - T1_eV) * Ha2eV;
    cout.precision(10);
    cout<<"T1: "<<T1_eV<<" T2: "<<T2_eV<<endl;
    cout<<"E1: "<<E1_Ha<<" E2: "<<E2_Ha<<endl;
    cout<<"mu1: "<<mu1_eV<<" mu2: "<<mu2_eV<<endl;
    cout<<"Cv/N_i: "<<Cv<<" "<<yellow("k_b/atom")<<endl;
}

double FUNC::FEG_E(const double density_e, const double density_i, const double T_eV, const double mu_eV)
{
    const double n_i_au = density_i*pow(P_bohr*1e-8, 3);
    const double T_au = T_eV / Ha2eV;
    std::cout<<mu_eV/T_eV<<std::endl;
    return sqrt(2)/pow(M_PI,2) / n_i_au * pow(T_au, 2.5) * fd_integral(mu_eV/T_eV, 3.0/2.0);
}
#include "function.h"
#include "const.h"
#include "elements.h"
#include "tool.h"
#include "iomanip"
#include <iostream>

#define NT 60 

void FUNC::mdparameter()
{
	vector<double> mlist, zlist, nlist, zionlist;
	read_elemets(mlist, zlist, nlist, zionlist);
	double nmol = read_number_of_molecules();
	double density = read_density();
	double temp = read_temperature();
	//----------------------------------------------------------
	thomas_fermi_ionization(density, temp, mlist, zlist, nlist, zionlist);
	molecule mol(mlist, zionlist, nlist);
	molecule mol0(mlist, zlist, nlist);

	double z_per_mol = mol.tot_z;
	double n_per_mol = mol.tot_n;
	double min_mass = mol.min_m;
	double na = n_per_mol * nmol;
	double ne = z_per_mol * nmol;

	double ionization = mol.tot_z/mol0.tot_z;

	double mass_atom_g = mol.avg_m/P_NA; //averge mass of each atom: g
	double l_cm = pow(na*mass_atom_g/density,1.00/3);
	double l_ang= l_cm * 1e8;
	double l_bohr=l_ang/P_bohr;
	double WSr_ang=pow(mass_atom_g/density*3/4/M_PI,1.00/3)*1e8;
	double WSr_bohr=WSr_ang/P_bohr;
	double dt_fs=WSr_ang*1e-10/NT/sqrt(temp*P_qe/(0.001*min_mass/P_NA))*1e15;
	double dt_au=dt_fs/(rau2s*1e15);
	double mu_eV = FEG_mu(ne / pow(l_cm,3), temp);
	double mu_Ry = mu_eV / Ry2eV;
	double Ecut1_eV = FEG_ECUT1( mu_eV, temp, 1e-2);
	double Ecut1_Ry = Ecut1_eV / Ry2eV;
	double Ecut2_eV = FEG_ECUT2( mu_eV, temp, 1e-5);;
	double Ecut2_Ry = Ecut2_eV / Ry2eV;

	double mu0_eV = fermi_energy(ne / pow(l_cm,3));
	double degeneracy = temp / mu0_eV;
	double coupling = mol.avg_z / (2*WSr_bohr) / temp * Ha2eV;


	int nbands1 = this->calbands(Ecut1_eV, mu0_eV, mol0.tot_z*nmol, ionization);
	int nbands2 = this->calbands(Ecut2_eV, mu0_eV, mol0.tot_z*nmol, ionization);
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Lattice constant: "<<fixed<<setprecision(8)<<l_bohr<<yellow(" Bohr; ")<<l_ang<<yellow(" Angstrom")<<endl;
	cout<<"Wigner-Seitz radius: "<<WSr_bohr<<yellow(" Bohr;")<<"  (0.7*WS = "<<0.7*WSr_bohr<<")"<<endl;
	cout<<"Temperature: "<<temp/Ry2eV<<yellow(" Ry; ")<<temp*eV2K<<yellow(" K")<<endl;
	cout<<"dt: "<<dt_fs<<yellow(" fs; ")<<dt_au<<yellow(" Ry^-1")<<endl;
	cout<<"1/(40dt): "<<double(1.0)/40/dt_fs<<yellow(" fs^-1")<<endl;
	cout<<"FEG mu: "<<setprecision(3)<<mu_eV<<yellow(" eV; ")<<mu_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut1 (Delta N < 1e-3):   "<<Ecut1_eV<<yellow(" eV; ")<<Ecut1_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut2 (f or occ. < 1e-5): "<<Ecut2_eV<<yellow(" eV; ")<<Ecut2_Ry<<yellow(" Ry")<<endl;
	cout<<"Ionization: "<<ionization*100<<"%"<<endl;
	cout<<"Nbands1: "<<nbands1<<" ; Nbands2: "<<nbands2<<endl;
	cout<<"Coupling constant: "<<coupling<<" ; Degeneracy parameter: "<<degeneracy<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	return;
}
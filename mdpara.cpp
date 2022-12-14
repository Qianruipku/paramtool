#include "function.h"
#include "const.h"
#include "elements.h"
#include "tool.h"
#include "iomanip"
#include <iostream>

#define NT 60 

void FUNC::mdparameter()
{
	vector<double> mlist, zlist, nlist;
	read_elemets(mlist, zlist, nlist);
	double nmol = read_number_of_molecules();
	double density = read_density();
	double temp = read_temperature();
	//----------------------------------------------------------
	vector<double> zionlist = thomas_fermi_ionization(density, temp, mlist, zlist, nlist);
	molecule mol(mlist, zionlist, nlist);

	double z_per_mol = mol.tot_z;
	double n_per_mol = mol.tot_n;
	double min_mass = mol.min_m;
	double na = n_per_mol * nmol;
	double ne = z_per_mol * nmol;

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
	double Ecut1_eV = FEG_ECUT1( mu_eV, temp);
	double Ecut1_Ry = Ecut1_eV / Ry2eV;
	double Ecut2_eV = FEG_ECUT2( mu_eV, temp);;
	double Ecut2_Ry = Ecut2_eV / Ry2eV;
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Lattice constant: "<<fixed<<setprecision(8)<<l_bohr<<yellow(" P_bohr; ")<<l_ang<<yellow(" Angstrom")<<endl;
	cout<<"Wigner-Seitz radius: "<<WSr_bohr<<yellow(" P_bohr;")<<"  (0.7*WS = "<<0.7*WSr_bohr<<")"<<endl;
	cout<<"Temperature: "<<temp/Ry2eV<<yellow(" Ry; ")<<temp*eV2K<<yellow(" K")<<endl;
	cout<<"dt: "<<dt_fs<<yellow(" fs; ")<<dt_au<<yellow(" Ry.-a.u.")<<endl;
	cout<<"1/(40dt): "<<double(1.0)/40/dt_fs<<yellow(" fs^-1")<<endl;
	cout<<"FEG mu: "<<setprecision(3)<<mu_eV<<yellow(" eV; ")<<mu_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut1 (interror < 1e-3): "<<Ecut1_eV<<yellow(" eV; ")<<Ecut1_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut2 (f<1e-5):          "<<Ecut2_eV<<yellow(" eV; ")<<Ecut2_Ry<<yellow(" Ry")<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	return;
}
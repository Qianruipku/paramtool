#include "function.h"
#include "const.h"
#include "elements.h"
#include "tool.h"
#include "iomanip"
#include <iostream>

#define NT 60 

void FUNC::mdparameter()
{
	vector<double> mlist, qlist, nlist;
	read_elemets(mlist, qlist, nlist);
	double z_per_mol(0), mass_per_mol(0), n_per_mol(0), min_mass(1e5);
	for(int i = 0 ; i < nlist.size(); ++i)
	{
		min_mass = min_mass > mlist[i] ? mlist[i] : min_mass;
        mass_per_mol += nlist[i] * mlist[i];
        z_per_mol += nlist[i] * qlist[i];
		n_per_mol += nlist[i];
	}
	// cout<<" "<<mass_per_mol<<" "<<z_per_mol<<" "<<n_per_mol<<" "<<min_mass<<endl;
	
	double nmol = read_number_of_molecules();
	double na = n_per_mol * nmol;
	double ne = z_per_mol * nmol;

	double density = read_density();
	double temp = read_temperature();

	double mass_atom_g = mass_per_mol/P_NA/n_per_mol;
	double l_ang=pow(na*mass_atom_g/density,1.00/3)*1e8;
	double l_bohr=l_ang/P_bohr;
	double WSr_ang=pow(mass_atom_g/density*3/4/M_PI,1.00/3)*1e8;
	double WSr_bohr=WSr_ang/P_bohr;
	double dt_fs=WSr_ang*1e-10/NT/sqrt(temp*P_qe/(0.001*min_mass))*1e15;
	double dt_au=dt_fs/(au2s*1e15);

	double mu_eV = FEG_mu(pow(l_bohr,3), ne, temp);
	double mu_Ry = mu_eV / Ry2eV;
	double Ecut1_eV = FEG_ECUT1( mu_eV, temp);
	double Ecut1_Ry = Ecut1_eV / Ry2eV;
	double Ecut2_eV = FEG_ECUT2( mu_eV, temp);;
	double Ecut2_Ry = Ecut2_eV / Ry2eV;
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Lattice constant: "<<fixed<<setprecision(8)<<l_bohr<<yellow(" P_bohr; ")<<l_ang<<yellow(" Angstrom")<<endl;
	cout<<"Wigner-Seitz radius: "<<WSr_bohr<<yellow(" P_bohr;")<<"  (0.7*WS = "<<0.7*WSr_bohr<<")"<<endl;
	cout<<"Temperature: "<<temp/Ry2eV<<yellow(" Ry; ")<<temp*eV2K<<yellow(" K")<<endl;
	cout<<"dt: "<<dt_fs<<yellow(" fs; ")<<dt_au<<yellow(" a.u.")<<endl;
	cout<<"1/(40dt): "<<double(1.0)/40/dt_fs<<yellow(" fs^-1")<<endl;
	cout<<"FEG mu: "<<setprecision(3)<<mu_eV<<yellow(" eV; ")<<mu_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut1 (interror < 1e-3): "<<Ecut1_eV<<yellow(" eV; ")<<Ecut1_Ry<<yellow(" Ry")<<endl;
	cout<<"Guess Ecut2 (f<1e-5):          "<<Ecut2_eV<<yellow(" eV; ")<<Ecut2_Ry<<yellow(" Ry")<<endl;
	cout<<"---------------------------------------------------------"<<endl;
	return;
}
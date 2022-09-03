#include "function.h"
#include "const.h"
#include "tool.h"
#include "iomanip"
#include <iostream>
using namespace std;
#define NT 60 

void FUNC::mdparameter()
{
	int na;
	double density;
	int rank;
	double mass;
	double temp;
	cout<<green("Input the rank of element")<<endl;
	cin>>rank;
	switch(rank)
	{
		case 1:
			mass = ele_H;
			break;
		case 6:
			mass=ele_C;
			break;
		case 13:
			mass = ele_Al;
			break;
		case 14:
			mass=ele_Si;
			break;
		default:
			cout<<red("No such RANK yet")<<endl;
			exit(1);
	}
	cout<<green("Input the number of atoms:")<<endl;
	cin>>na;
	cout<<green("Input the density (g/cm^3):")<<endl;
    cin>>density;
	cout<<green("Input the temperature (eV):")<<endl;
    cin>>temp;

	double mass_atom_g = mass/NA;
	double l_ang=pow(na*mass_atom_g/density,1.00/3)*1e8;
	double l_bohr=l_ang/bohr;
	double WSr_ang=pow(mass_atom_g/density*3/4/M_PI,1.00/3)*1e8;
	double WSr_bohr=WSr_ang/bohr;
	double dt_fs=WSr_ang*1e-10/NT/sqrt(temp*qe/(0.001*mass_atom_g))*1e15;
	double dt_au=dt_fs/(au2s*1e15);
	cout<<"---------------------------------------------------------"<<endl;
	cout<<"Lattice constant: "<<fixed<<setprecision(6)<<l_bohr<<yellow(" Bohr; ")<<l_ang<<yellow(" Angstrom")<<endl;
	cout<<"Wigner-Seitz radius: "<<WSr_bohr<<yellow(" Bohr;")<<"  (0.7*WS = "<<0.7*WSr_bohr<<")"<<endl;
	cout<<"Temperature: "<<temp/Ry2eV<<yellow(" Ry; ")<<temp*eV2K<<yellow(" K")<<endl;
	cout<<"dt: "<<dt_fs<<yellow(" fs; ")<<dt_au<<yellow(" a.u.")<<endl;
	cout<<"1/(40dt): "<<double(1.0)/40/dt_fs<<yellow(" fs^-1")<<endl;
	cout<<"---------------------------------------------------------"<<endl;
		
	return;
}


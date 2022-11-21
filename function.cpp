#include "function.h"
#include "tool.h"
#include <iostream>
#include "elements.h"

double FUNC:: read_density()
{
    double density_gcm;
    cout<<green("Input the density (g/cm^3):")<<endl;
    cin>>density_gcm;
    return density_gcm;
}
double FUNC:: read_temperature()
{
    double temperature_eV;
    cout<<green("Input the temperature (eV):")<<endl;
    cin>>temperature_eV;
    return temperature_eV;
}

double FUNC:: read_number_of_molecules()
{
    double nmol;
    cout<<green("Input the number of molecules:")<<endl;
	cin>>nmol;
    return nmol;
}

double FUNC::read_elemets(vector<double> &mlist, vector<double> &zlist, vector<double> &nlist)
{
    string name;
    cout<<green("Input the name of molecule (e.g. Mg-O, C-O-2, H-2-O)")<<endl;
	cin>>name;
    readstring(name, mlist, zlist, nlist);
}
#ifndef MOLECULE_H
#define MOLECULE_H
#include <vector>
using namespace std;
class molecule
{
public:
    molecule(){};
    molecule(const vector<double> mlist_in, const vector<double> zlist_in, const vector<double> nlist_in);
    molecule(molecule& mol);
    ~molecule(){};
    //setup this class
    void setup();
    //relative atomic mass of each atom
    vector<double> mlist;
    //number of ionized electrons of each atom
    vector<double> zlist;
    //number of atoms
    vector<double> nlist;
    //totol mass, z, n
    double tot_m, tot_z, tot_n;
    //averge mass and z
    double avg_m, avg_z;
    //number of elements
    double nele;

    double min_m, max_m;
    molecule operator=(const molecule& mol);

};
#endif
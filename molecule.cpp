#include "molecule.h"
molecule::molecule(const vector<double> mlist_in, const vector<double> zlist_in, const vector<double> nlist_in)
{
    this->mlist = mlist_in;
    this->zlist = zlist_in;
    this->nlist = nlist_in;
    this->setup();
}

void molecule:: setup()
{
    this->nele = min(min(nlist.size(), zlist.size()), mlist.size());
    this->avg_m = this->avg_z = this->tot_m = this->tot_n = this->tot_z = 0;
    this->min_m = 1e5;
    this->max_m = 0;
    for(int i = 0 ; i < this->nele; ++i)
	{
        this->tot_m += nlist[i] * mlist[i];
        this->tot_z += nlist[i] * zlist[i];
        this->tot_n += nlist[i];
        this->min_m = min(min_m, mlist[i]);
        this->max_m = max(max_m, mlist[i]);
    }
    this->avg_m = this->tot_m / this->tot_n;
    this->avg_z = this->tot_z / this->tot_n;
}

molecule:: molecule(molecule& mol)
{
    this->operator=(mol);
}

molecule molecule:: operator=(const molecule mol)
{
    this->mlist = mol.mlist;
    this->zlist = mol.zlist;
    this->nlist = mol.nlist;
    this->nele = mol.nele;
    this->avg_m = mol.avg_m;
    this->avg_z = mol.avg_z;
    this->tot_m = mol.tot_m;
    this->tot_z = mol.tot_z;
    this->tot_n = mol.tot_n;
    this->min_m = mol.min_m;
    this->max_m = mol.max_m;
    return *this;
}

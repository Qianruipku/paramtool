#include "elements.h"
#include <vector>
#include "tool.h"
#include <iostream>
#include <algorithm>

void Stringsplit(const string& str, const char split, vector<string>& res);
void get_masscharge(const string in, double &mass, double &charge);


void readstring(const string name, vector<double> &mlist, vector<double> &zlist, vector<double> &nlist, vector<double> &zionlist)
{
	mlist.clear();
	zlist.clear();
	nlist.clear();
	zionlist.clear();
    vector<string> list;
    Stringsplit(name, '-', list);
    if(list.size() == 0)
    {
        cout<<red("No information")<<endl;
        exit(0);
    }
    double tmpmass, tmpchg, tmpzion;
    for(vector<string>::iterator it = list.begin(); it != list.end() ;++it)
    {
        double tmpit = str2dou(*it);
        if(tmpit > 1e-6)
        {
            cout<<red("Wrong format!")<<endl;
            exit(0);
        }
		
		string namepluszion = *it;
		size_t pos = namepluszion.find('_');
		string name;
		if( pos != namepluszion.npos)
		{
			name = namepluszion.substr(0, pos);
			tmpzion = str2dou(namepluszion.substr(pos+1, namepluszion.size()));
		}
		else
		{
			name = namepluszion;
			tmpzion = 0.0;
		}

        get_masscharge(name, tmpmass, tmpchg);
        double count = 1.0;
        ++it;
        if(it != list.end())
        {
            double numit = str2dou(*it);
            if(numit > 1e-6)
                count = numit;
            else
                --it;
        }
        else
            --it;
		mlist.push_back(tmpmass);
		zlist.push_back(tmpchg);
		nlist.push_back(count);
		zionlist.push_back(tmpzion);
		// cout<<tmpmass<<" "<<tmpchg<<" "<<count<<endl;
    }
}

void Stringsplit(const string& str, const char split, vector<string>& res)
{
	if (str == "")		return;
	string strs = str + split;
	size_t pos = strs.find(split);

	while (pos != strs.npos)
	{
		string temp = strs.substr(0, pos);
		res.push_back(temp);
		strs = strs.substr(pos + 1, strs.size());
		pos = strs.find(split);
	}
    return;
}

void get_masscharge(const string in, double &mass, double &charge)
{
    string small = in;
    transform(small.begin(),small.end(),small.begin(),::tolower);
    if(small == "h")
	{
		mass = MASS_H;
        charge = 1.0;
	}
	else if(small == "d")
	{
		mass = MASS_D;
        charge = 1.0;
	}
	else if(small == "t")
	{
		mass = MASS_T;
        charge = 1.0;
	}
	else if(small == "he")
	{
		mass = MASS_He;
        charge = 2.0;
	}
	else if(small == "li")
	{
		mass = MASS_Li;
        charge = 3.0;
	}
	else if(small == "be")
	{
		mass = MASS_Be;
        charge = 4.0;
	}
	else if(small == "b")
	{
		mass = MASS_B;
        charge = 5.0;
	}
	else if(small == "c")
	{
		mass = MASS_C;
        charge = 6.0;
	}
	else if(small == "n")
	{
		mass = MASS_N;
        charge = 7.0;
	}
	else if(small == "o")
	{
		mass = MASS_O;
        charge = 8.0;
	}
	else if(small == "f")
	{
		mass = MASS_F;
        charge = 9.0;
	}
	else if(small == "na")
	{
		mass = MASS_Na;
        charge = 11.0;
	}
	else if(small == "mg")
	{
		mass = MASS_Mg;
        charge = 12.0;
	}
	else if(small == "al")
	{
		mass = MASS_Al;
        charge = 13.0;
	}
	else if(small == "si")
	{
		mass = MASS_Si;
        charge = 14.0;
	}
	else if(small == "p")
	{
		mass = MASS_P;
        charge = 15.0;
	}
	else if(small == "s")
	{
		mass = MASS_S;
        charge = 16.0;
	}
	else if(small == "cl")
	{
		mass = MASS_Cl;
        charge = 17.0;
	}
	else if(small == "ar")
	{
		mass = MASS_Ar;
        charge = 18.0;
	}
	else if(small == "k")
	{
		mass = MASS_K;
        charge = 19.0;
	}
	else if(small == "ca")
	{
		mass = MASS_Ca;
        charge = 20.0;
	}
	else if(small == "sc")
	{
		mass = MASS_Sc;
        charge = 21.0;
	}
	else if(small == "ti")
	{
		mass = MASS_Ti;
        charge = 22.0;
	}
	else if(small == "v")
	{
		mass = MASS_V;
        charge = 23.0;
	}
	else if(small == "cr")
	{
		mass = MASS_Cr;
        charge = 24.0;
	}
	else if(small == "mn")
	{
		mass = MASS_Mn;
        charge = 25.0;
	}
	else if(small == "fe")
	{
		mass = MASS_Fe;
        charge = 26.0;
	}
	else if(small == "co")
	{
		mass = MASS_Co;
        charge = 27.0;
	}
	else if(small == "ni")
	{
		mass = MASS_Ni;
        charge = 28.0;
	}
	else if(small == "cu")
	{
		mass = MASS_Cu;
        charge = 29.0;
	}
	else if(small == "zn")
	{
		mass = MASS_Zn;
        charge = 30.0;
	}
	else if(small == "ga")
	{
		mass = MASS_Ga;
        charge = 31.0;
	}
	else if(small == "ge")
	{
		mass = MASS_Ge;
        charge = 32.0;
	}
	else if(small == "as")
	{
		mass = MASS_As;
        charge = 33.0;
	}
	else if(small == "se")
	{
		mass = MASS_Se;
        charge = 34.0;
	}
	else if(small == "br")
	{
		mass = MASS_Br;
        charge = 35.0;
	}
	else if(small == "kr")
	{
		mass = MASS_Kr;
        charge = 36.0;
	}
	else if(small == "rb")
	{
		mass = MASS_Rb;
        charge = 37.0;
	}
	else if(small == "sr")
	{
		mass = MASS_Sr;
        charge = 38.0;
	}
	else if(small == "y")
	{
		mass = MASS_Y;
        charge = 39.0;
	}
	else if(small == "zr")
	{
		mass = MASS_Zr;
        charge = 40.0;
	}
	else if(small == "nb")
	{
		mass = MASS_Nb;
        charge = 41.0;
	}
	else if(small == "mo")
	{
		mass = MASS_Mo;
        charge = 42.0;
	}
	else if(small == "tc")
	{
		mass = MASS_Tc;
        charge = 43.0;
	}
	else if(small == "ru")
	{
		mass = MASS_Ru;
        charge = 44.0;
	}
	else if(small == "rh")
	{
		mass = MASS_Rh;
        charge = 45.0;
	}
	else if(small == "pd")
	{
		mass = MASS_Pd;
        charge = 46.0;
	}
	else if(small == "ag")
	{
		mass = MASS_Ag;
        charge = 47.0;
	}
	else if(small == "cd")
	{
		mass = MASS_Cd;
        charge = 48.0;
	}
	else if(small == "in")
	{
		mass = MASS_In;
        charge = 49.0;
	}
	else if(small == "sn")
	{
		mass = MASS_Sn;
        charge = 50.0;
	}
	else if(small == "sb")
	{
		mass = MASS_Sb;
        charge = 51.0;
	}
	else if(small == "te")
	{
		mass = MASS_Te;
        charge = 52.0;
	}
	else if(small == "i")
	{
		mass = MASS_I;
        charge = 53.0;
	}
	else if(small == "xe")
	{
		mass = MASS_Xe;
        charge = 54.0;
	}
	else if(small == "cs")
	{
		mass = MASS_Cs;
        charge = 55.0;
	}
	else if(small == "ba")
	{
		mass = MASS_Ba;
        charge = 56.0;
	}
	else if(small == "la")
	{
		mass = MASS_La;
        charge = 57.0;
	}
	else if(small == "ce")
	{
		mass = MASS_Ce;
        charge = 58.0;
	}
	else if(small == "pr")
	{
		mass = MASS_Pr;
        charge = 59.0;
	}
	else if(small == "nd")
	{
		mass = MASS_Nd;
        charge = 60.0;
	}
	else if(small == "pm")
	{
		mass = MASS_Pm;
        charge = 61.0;
	}
	else if(small == "sm")
	{
		mass = MASS_Sm;
        charge = 62.0;
	}
	else if(small == "eu")
	{
		mass = MASS_Eu;
        charge = 63.0;
	}
	else if(small == "gd")
	{
		mass = MASS_Gd;
        charge = 64.0;
	}
	else if(small == "tb")
	{
		mass = MASS_Tb;
        charge = 65.0;
	}
	else if(small == "dy")
	{
		mass = MASS_Dy;
        charge = 66.0;
	}
	else if(small == "ho")
	{
		mass = MASS_Ho;
        charge = 67.0;
	}
	else if(small == "er")
	{
		mass = MASS_Er;
        charge = 68.0;
	}
	else if(small == "tm")
	{
		mass = MASS_Tm;
        charge = 69.0;
	}
	else if(small == "yb")
	{
		mass = MASS_Yb;
        charge = 70.0;
	}
	else if(small == "lu")
	{
		mass = MASS_Lu;
        charge = 71.0;
	}
	else if(small == "hf")
	{
		mass = MASS_Hf;
        charge = 72.0;
	}
	else if(small == "ta")
	{
		mass = MASS_Ta;
        charge = 73.0;
	}
	else if(small == "w")
	{
		mass = MASS_W;
        charge = 74.0;
	}
	else if(small == "re")
	{
		mass = MASS_Re;
        charge = 75.0;
	}
	else if(small == "os")
	{
		mass = MASS_Os;
        charge = 76.0;
	}
	else if(small == "ir")
	{
		mass = MASS_Ir;
        charge = 77.0;
	}
	else if(small == "pt")
	{
		mass = MASS_Pt;
        charge = 78.0;
	}
	else if(small == "au")
	{
		mass = MASS_Au;
        charge = 79.0;
	}
	else if(small == "hg")
	{
		mass = MASS_Hg;
        charge = 80.0;
	}
	else if(small == "tl")
	{
		mass = MASS_Tl;
        charge = 81.0;
	}
	else if(small == "pb")
	{
		mass = MASS_Pb;
        charge = 82.0;
	}
	else if(small == "bi")
	{
		mass = MASS_Bi;
        charge = 83.0;
	}
	else if(small == "po")
	{
		mass = MASS_Po;
        charge = 84.0;
	}
	else if(small == "at")
	{
		mass = MASS_At;
        charge = 85.0;
	}
	else if(small == "rn")
	{
		mass = MASS_Rn;
        charge = 86.0;
	}
	else if(small == "fr")
	{
		mass = MASS_Fr;
        charge = 87.0;
	}
	else if(small == "ra")
	{
		mass = MASS_Ra;
        charge = 88.0;
	}
	else
	{
		cout<<green("We do not store the mass of ")<<in<<green(" yet!")<<endl;
        cout<<green("Please give the relative molecular mass of ")<<in<<green(":")<<endl;
        cin>>mass;
        cout<<green("Please give the number of electrons in ")<<in<<green(":")<<endl;
        cin>>charge;
	}
}
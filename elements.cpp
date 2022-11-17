#include "elements.h"
#include <vector>
#include "tool.h"
#include <iostream>
#include <algorithm>

using namespace std;
void Stringsplit(const string& str, const char split, vector<string>& res);
void get_masscharge(const string in, double &mass, double &charge);


bool readstring(const string txt, double & mass, double & q_tot)
{
    vector<string> list;
    Stringsplit(txt, '-', list);
    if(list.size() == 0)
    {
        cout<<red("No information")<<endl;
        return false;
    }
    double tmpmass, tmpchg;
    mass = q_tot = 0.0;
    for(vector<string>::iterator it = list.begin(); it != list.end() ;++it)
    {
        double tmpit = str2dou(*it);
        if(tmpit > 1e-6)
        {
            cout<<red("Wrong format!")<<endl;
            return false;
        }
        get_masscharge(*it, tmpmass, tmpchg);
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
        mass += count * tmpmass;
        q_tot += count * tmpchg;
    }
    return true;
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
	else
	{
		cout<<green("We do not store the mass of ")<<in<<green(" yet!")<<endl;
        cout<<green("Please give the relative molecular mass of ")<<in<<green(":")<<endl;
        cin>>mass;
        cout<<green("Please give the number of electrons in ")<<in<<green(":")<<endl;
        cin>>charge;
	}
}
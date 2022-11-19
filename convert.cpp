#include "function.h"
#include "tool.h"
#include "const.h"
#include <iostream>
#include "iomanip"

void FUNC::unitcovert()
{
    int unittype = 0;
    cout<<green("Please choose: ")<<endl;
    cout<<green("1: temperature")<<endl;
    cin>>unittype;
    if(unittype==1)
		convertT();
	else
	{
		cout<<red("Wrong Input!")<<endl;
		exit(1);
	}

}

void FUNC::convertT()
{
    int Ttype = 0;
    double rawdata = 0; 
    cout<<green("Please choose the unit: ")<<endl;
    cout<<green("1: eV")<<endl;
    cout<<green("2: K")<<endl;
    cout<<green("3: Ry")<<endl;
    cin>>Ttype;
    cout<<green("Please input the number:  ");
    cin>>rawdata;
    
    double UeV, UK, URy;
    if(Ttype == 1)
    {
        UeV = rawdata;
        UK = UeV * eV2K;
        URy = UeV / Ry2eV;
    }
    else if(Ttype == 2)
    {
        UK = rawdata;
        UeV = UK / eV2K;
        URy = UeV / Ry2eV;
    }
    else if(Ttype == 3)
    {
        URy = rawdata;
        UeV = URy * Ry2eV;
        UK = UeV * eV2K;
    }
    else
    {
        cout<<red("Wrong Input!")<<endl;
		exit(1);
    }
    cout<<"The temperature is:"<<endl;
    cout<<setprecision(15)<<UK<<yellow(" K; ")<<UeV<<yellow(" eV; ")<<URy<<yellow(" Ry;")<<endl;

}


#include "function.h"
#include "tool.h"
#include <iostream>
using namespace std;

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
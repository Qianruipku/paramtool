#include <iostream>
#include "tool.h"
#include "function.h"

int main()
{
	int func=0;
	cout<<green("Please choose the tool ")<<endl;
	cout<<green("1: Unit convert")<<endl;
	cout<<green("2: calculate MD parameter")<<endl;
	cout<<green("3: Unit convert")<<endl;
	FUNC test;
	cin>>func;
	switch(func)
	{
		case 1:
			test.unitcovert();
			break;
		case 2:
			test.mdparameter();
			break;
		case 3:
			test.conductivity();
			break;
		default:
			cout<<red("Wrong Input!")<<endl;
			exit(1);
	}
	
	return 0;
}

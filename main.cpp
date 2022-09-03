#include <iostream>
#include "tool.h"
#include "function.h"

using namespace std;

int main()
{
	int func=0;
	cout<<green("Please choose the tool ")<<endl;
	cout<<green("1: calculate MD parameter")<<endl;
	cout<<green("2: Unit convert")<<endl;
	FUNC test;
	cin>>func;
	switch(func)
	{
		case 1:
			test.mdparameter();
			break;
		case 2:
			test.unitcovert();
			break;
		default:
			cout<<red("Wrong Input!")<<endl;
			exit(1);
	}
	
	return 0;
}

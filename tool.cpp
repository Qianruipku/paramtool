#include "tool.h"
#include <iostream>

string red(string txt)
{
	string out = "\033[31m"+txt+"\033[0m";
	return out;
}

string green(string txt)
{
	string out = "\033[32m"+txt+"\033[0m";
	return out;
}

string yellow(string txt)
{
	string out = "\033[33m"+txt+"\033[0m";
	return out;
}
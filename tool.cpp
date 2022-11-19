#include "tool.h"
#include <iostream>
#include <sstream>

string red(string txt)
{
#ifdef __COLOR
	string out = "\033[31m"+txt+"\033[0m";
	return out;
#else  
    return txt;
#endif
}

string green(string txt)
{
#ifdef __COLOR
	string out = "\033[32m"+txt+"\033[0m";
	return out;
#else  
    return txt;
#endif
}

string yellow(string txt)
{
#ifdef __COLOR
	string out = "\033[33m"+txt+"\033[0m";
	return out;
#else  
    return txt;
#endif
}

int str2int(const string str)
{
    stringstream ss;
    int y;
    ss<<str;
    ss>>y;
    return y;
}
string int2str(const int y)
{
    stringstream ss;
    string str;
    ss<<y;
    ss>>str;
    return str;
}
string dou2str(const double y)
{
    stringstream ss;
    string str;
    ss<<y;
    ss>>str;
    return str;
}
double str2dou(const string str)
{
    stringstream ss;
    double y;
    ss<<str;
    ss>>y;
    return y;
}
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

void tridiagsolver(std::vector<double>& x, const std::vector<double>& a, const std::vector<double>& b, 
             const std::vector<double>& c, const std::vector<double>& d)
{
    int num = d.size();
    double *l = new double[num];
    double *y = new double[num];
    double *u = new double[num - 1];

    // Forward reduction
    l[0] = b[0];
    y[0] = d[0] / l[0];
    u[0] = c[0] / l[0];
    for (int i = 1; i < num; i++)
    {
        l[i] = b[i] - a[i-1] * u[i-1];
        y[i] = (d[i] - y[i-1] * a[i-1]) / l[i];
        if (i < num - 1)
            u[i] = c[i] / l[i];
    }

    // Forward reduction
    x[num-1] = y[num-1];
    for (int i = num-2; i >= 0; i--)
    {
        x[i] = y[i] - u[i] * x[i+1];
    }

    delete [] l;
    delete [] y;
    delete [] u;
}

void NaturalSplineInterpolation(const std::vector<double>& x,           std::vector<double>& y, 
                                const std::vector<double>& ref_x, const std::vector<double>& ref_y) 
{
    int ref_dim = ref_x.size();
    int data_dim = x.size();
    // Check if the input vectors have the same size
    if (ref_dim != ref_y.size() || data_dim != y.size()) {
        throw std::runtime_error("Input vector sizes are not matching");
    }
    //  Define vectors h and delta
    std::vector<double> h(ref_dim), delta(ref_dim);
    for (int i = 0; i < ref_dim - 1; ++i) {
        h[i] = ref_x[i + 1] - ref_x[i];
        delta[i] = (ref_y[i + 1] - ref_y[i]) / h[i];
    }

    // Build the tridiagonal matrix
    // Tridiag(a,b,c)*M=D
    std::vector<double> a(ref_dim), b(ref_dim), c(ref_dim), d(ref_dim);
    c[0] = 0.0;
    b[0] = 1.0;
    d[0] = 0.0;
    
    for (int i = 1; i < ref_dim - 1; ++i) {
        a[i-1] = h[i - 1];
        c[i] = h[i];
        b[i] = 2 * (h[i - 1] + h[i]);
        d[i] = 6 * (delta[i] - delta[i - 1]);
    }

    a[ref_dim-1] = 0;
    b[ref_dim-1] = 1;
    d[ref_dim-1] = 0;

    // Solve the tridiagonal matrix to get the second derivatives
    std::vector<double> m(ref_dim);
    tridiagsolver(m,a,b,c,d);

    int index = 0;
    for(int i = 0; i < data_dim; ++i)
    {
        double tmpx = x[i];
        // Find the interpolation interval
        while (index < ref_dim - 1 && tmpx > ref_x[index + 1]) {
            ++index;
        }

        // Compute the interpolation result
        double t = tmpx - ref_x[index];
        double a0 = ref_y[index];
        double a1 = delta[index] - h[index] * (2 * m[index] + m[index + 1]) / 6.0;
        double a2 = m[index] / 2.0;
        double a3 = (m[index + 1] - m[index]) / (6.0 * h[index]);
        y[i] = a0 + a1 * t + a2 * t * t + a3 * t * t * t;

    }
}
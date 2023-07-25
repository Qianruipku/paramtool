#ifndef TOOL_H
#define TOOL_H
#include <string>
#include <vector>
using namespace std;

string red(string txt);

string green(string txt);

string yellow(string txt);

int str2int(const string str);
string int2str(const int y);
double str2dou(const string str);
string dou2str(const double y);

/**
 * @brief a=(a0,a2,...,an-2), b=(b0,b1,...,bn-1), c=(c0,c1,...,cn-2), A=diag(a,b,c)
 *                            d=(d0,d1,...,dn-1)
 *        Solve Ax=d, where A is a tridiagonal matrix with diagonals a, b, and c
 * 
 * @param x Solution vector
 * @param a Sub-diagonal elements of A
 * @param b Diagonal elements of A
 * @param c Super-diagonal elements of A
 * @param d Right-hand side vector
 */
void tridiagsolver(std::vector<double>& x, const std::vector<double>& a, const std::vector<double>& b, 
             const std::vector<double>& c, const std::vector<double>& d);

/**
 * @brief Natural Spline Interpolation function
 *        Free boundary (Natural), with both end second derivatives equal to 0. S“(ref_x_0 or ref_x_n-1) = 0
 * @param x [in]  x points
 * @param y [out] y points to be interpolated
 * @param ref_x [in] ref x data
 * @param ref_y [in] ref y data
 */
void NaturalSplineInterpolation(const std::vector<double>& x,           std::vector<double>& y, 
                                const std::vector<double>& ref_x, const std::vector<double>& ref_y);

#endif
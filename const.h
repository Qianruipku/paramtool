#include <complex>
#include <cmath>
using namespace std;

#ifndef CONST_H
#define CONST_H

//basis
const double bohr = 0.52917721067;
const double NA = 6.02214076e23;
const double qe = 1.6021766208e-19;
const double Me = 9.1093897e-31;
const double kB = 1.3806505e-23;
constexpr double hbar = 6.62607015e-34/2/M_PI;
//atomic unit
constexpr double Ha2eV = hbar*hbar/(bohr*bohr*1e-20)/Me/qe; //27.21136857564
constexpr double Ry2eV = Ha2eV/2; //13.60568428782
constexpr double au2s = hbar/Ry2eV/qe; //4.837771834548454e-17

constexpr double eV2K = qe/kB;  //11604.518026

#endif
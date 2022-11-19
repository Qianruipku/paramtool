#include <complex>
#include <cmath>
using namespace std;

#ifndef CONST_H
#define CONST_H

//basis
constexpr double P_bohr = 0.52917721067;
constexpr double P_NA = 6.02214076e23;
constexpr double P_qe = 1.6021766208e-19;
constexpr double P_Me = 9.1093897e-31;
constexpr double P_kB = 1.3806505e-23;
constexpr double P_h = 6.62607015e-34;
constexpr double P_hbar = P_h/2/M_PI; //1.054571817646e-34
//atomic unit
constexpr double Ha2eV = P_hbar*P_hbar/(P_bohr*P_bohr*1e-20)/P_Me/P_qe; //27.21136857564
constexpr double Ry2eV = Ha2eV/2; //13.60568428782
constexpr double au2s = P_hbar/Ry2eV/P_qe; //4.837771834548454e-17

constexpr double eV2K = P_qe/P_kB;  //11604.518026

#endif
#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <string>
#include <vector>
using namespace std;

#define MASS_H   1.00794
#define MASS_D   2.01410
#define MASS_T   3.016
#define MASS_He  4.002602
#define MASS_Li  6.941
#define MASS_Be  9.012182
#define MASS_B   10.811
#define MASS_C   12.0107
#define MASS_O   15.9994
#define MASS_F   18.9984032
#define MASS_Ne  20.1797
#define MASS_Na  22.98976928
#define MASS_Mg  24.3050
#define MASS_Al  26.9815386
#define MASS_Si  28.0855
#define MASS_P   30.973762
#define MASS_S   32.065
#define MASS_Cl  35.453
#define MASS_Ar  39.948
#define MASS_K   39.0983
#define MASS_Ca  40.078
#define MASS_Sc  44.955912
#define MASS_Ti  47.867
#define MASS_V   50.9415
#define MASS_Cr  51.9961
#define MASS_Mn  54.938045
#define MASS_Fe  55.845
#define MASS_Co  58.933195
#define MASS_Ni  58.6934
#define MASS_Cu  63.546
#define MASS_Zn  65.38
#define MASS_Ga  69.723
#define MASS_Ge  72.64
#define MASS_As  74.92160
#define MASS_Se  78.96
#define MASS_Br  79.904
#define MASS_Kr  83.798
#define MASS_Rb  85.4678
#define MASS_Sr  87.62
#define MASS_Y   88.90585
#define MASS_Zr  91.224
#define MASS_Nb  92.90638
#define MASS_Mo  95.94
#define MASS_Tc  97.9072
#define MASS_Ru  101.07
#define MASS_Rh  102.90550
#define MASS_Pd  106.42
#define MASS_Ag  107.8682
#define MASS_Cd  112.411
#define MASS_In  114.818
#define MASS_Sn  118.710
#define MASS_Sb  121.760
#define MASS_Te  127.60
#define MASS_I   126.90447
#define MASS_Xe  131.293
#define MASS_Cs  132.9054519
#define MASS_Ba  137.327
#define MASS_La  138.90547
#define MASS_Ce  140.116
#define MASS_Pr  140.90765
#define MASS_Nd  144.242
#define MASS_Pm  145
#define MASS_Sm  150.36
#define MASS_Eu  151.964
#define MASS_Gd  157.25
#define MASS_Tb  158.92535
#define MASS_Dy  162.500
#define MASS_Ho  164.93032
#define MASS_Er  167.259
#define MASS_Tm  168.93421
#define MASS_Yb  173.04
#define MASS_Lu  174.967
#define MASS_Hf  178.49
#define MASS_Ta  180.94788
#define MASS_W   183.84
#define MASS_Re  186.207
#define MASS_Os  190.23
#define MASS_Ir  192.217
#define MASS_Pt  195.084
#define MASS_Au  196.966569
#define MASS_Hg  200.59
#define MASS_Tl  204.3833
#define MASS_Pb  207.2
#define MASS_Bi  208.98040
#define MASS_Po  208.9824
#define MASS_At  209.9871
#define MASS_Rn  222.0176
#define MASS_Fr  223
#define MASS_Ra  226


void readstring(const string name, vector<double> &mlist, vector<double> &zlist, vector<double> &nlist);
#endif
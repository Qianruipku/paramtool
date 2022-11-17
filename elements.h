#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <string>

using namespace std;

#define MASS_H  1.00794
#define MASS_D  2.01410
#define MASS_T  3.016
#define MASS_Li  6.941
#define MASS_Be  9.0121831
#define MASS_B  10.811
#define MASS_C  12.011
#define MASS_O  15.9994
#define MASS_F  18.9984
#define MASS_Na  22.98976928
#define MASS_Mg  24.3050
#define MASS_Al  26.981539
#define MASS_Si  28.0855

bool readstring(const string txt, double & mass, double & q_tot);

#endif
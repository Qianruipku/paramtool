#include "function.h"
#include "const.h"
double getmu(double mu_0, double T);
double calint(double fun(double e, double mu, double T), double mu, double T,double thr);
double funmu(double e, double mu, double T);

vector<double> FUNC::thomas_fermi_ionization(const double density_gm, const double T_eV, const vector<double> &mlist, const vector<double> &zlist, const vector<double> &nlist)
{
	double alpha = 14.3139;
    double beta = 0.6624;
    double a1 = 0.003323;
    double a2 = 0.9718;
    double a3 = 9.26148e-5;
    double a4 = 3.10165;
    double b0 = -1.7630;
    double b1 = 1.43175;
    double b2 = 0.31546;
    double c1 = -0.366667;
    double c2 = 0.983333;

	double m_per = 0, n_per_mol = 0;
	for(int i = 0 ; i < nlist.size(); ++i)
	{
        m_per += nlist[i] * mlist[i];
		n_per_mol += nlist[i];
	}
	m_per /= n_per_mol;

	vector<double> zionlist;
	for(int i = 0 ; i< nlist.size(); ++i)
	{
		double Z = zlist[i];
		double T0 = T_eV / pow(Z, 4.0/3.0);
    	double R = density_gm / (Z*m_per);
    	double TF = T0 / (1 + T0);
    	double A = a1*pow(T0, a2) + a3*pow(T0,a4);
    	double B = -exp(b0 + b1*TF + b2*pow(TF,7));
    	double C = c1*TF + c2;
    	double Q1 = A * pow(R, B);
    	double Q = pow((pow(R, C)+pow(Q1, C)), (1.0/C));
    	double x = alpha * pow(Q, beta);
		zionlist.push_back(Z * x / (1 + x + sqrt(1.0 + 2.0*x)));
	}

    return zionlist;
}

//mu of free electrons, assume all electrons are ionized.
//mu0=P_hbar^2/2m*(3pi^2N/V)^(2/3)
//mu = mu0(1 - pi^2/12*(kT/mu0)^2) (mu/T >> 1)
//density_e: cm^-3  ; T_eV: eV
double FUNC:: FEG_mu(const double density_e, const double T_eV)
{
	double density_e_bohr = density_e*pow(P_bohr*1e-8, 3);
    double mu0_eV = pow(3.0 * pow(M_PI,2) * density_e_bohr, 2.0/3.0) * Ry2eV; //unit in Ry
    double mu_eV;
    if(T_eV <= 1)
	{
		mu_eV = mu0_eV * (1 - pow(M_PI,2)/12.0*pow(T_eV/mu0_eV,2)-pow(M_PI,4)*7.0/960*pow(T_eV/mu0_eV,4));
	}
	else if(T_eV >= 1e4)
	{
		mu_eV = T_eV*(-log(6*pow(M_PI,2)) + 1.5*log(4*M_PI*mu0_eV/T_eV));
	}
	else
	{
		mu_eV = getmu(mu0_eV, T_eV);
	}
    return mu_eV;
}

double FUNC:: FEG_ECUT1(double mu_eV, double T_eV)
{
	double ref = calint(funmu, mu_eV, T_eV, 1e-20);
	double de = T_eV * 1e-4;
	double result = 1;
	double sum = funmu(0, mu_eV, T_eV);
	double e = 0;
	double diff = 1;
	while(diff > 1e-3)
	{
		e += de;
		sum += 4 * funmu(e, mu_eV, T_eV);
		e += de;
		sum += 2 * funmu(e,mu_eV,T_eV);
		diff = (ref - sum/3*de) / ref;
	}
	return e;
}
double FUNC:: FEG_ECUT2(double mu_eV, double T_eV)
{
    return mu_eV + 5*log(10.0)*T_eV;
}

//4/3 * mu_0^(3/2) = \int_0^\infty \frac{e^(1/2)}{exp((e-mu)/kT)+1} de
double getmu(double mu0, double T)
{
	double ref = 2.0/3.0 * pow(mu0, 3.0/2.0);
	double Deltamu = 5 * T;
	double mu1 = -Deltamu;
	double mu2 = Deltamu;
	double com1 = calint(funmu, mu1, T, 1e-8);
	double com2 = calint(funmu, mu2, T, 1e-8);
	while (com1 > ref)
    {
        mu2 = mu1;
        mu1 -= Deltamu;
        com1 = calint(funmu, mu1, T, 1e-8);
        Deltamu *= 10;
    }
    while (com2 < ref)
    {
        mu1 = mu2;
        mu2 += Deltamu;
        com2 = calint(funmu, mu2, T, 1e-8);
        Deltamu *= 10;
    }
	double diff = 1000;
	double mu3, com3;
	double thr = 1e-8 * ref;
	while (diff > thr)
    {
        mu3 = (mu2 + mu1) / 2;
        com3 = calint(funmu, mu3, T, 1e-8);
        if (com3 < ref)
        {
            mu1 = mu3;
        }
        else if (com3 > ref)
        {
            mu2 = mu3;
        }
        diff = abs(ref - com3);
    }
	return (mu2 + mu1) / 2;
}

double calint(double fun(double e, double mu, double T), double mu, double T, double thr)
{
	double de = T * 1e-4;
	double result = 1;
	double sum = fun(0, mu, T);
	double e = 0;
	int i = 0;
	while(result > thr || e <= T)
	{
		e += de;
		sum += 4 * fun(e, mu, T);
		e += de;
		result = fun(e, mu, T);
		sum += 2 * result;
	}
	sum += fun(e+de, mu, T);
	sum /= 3;
	return sum * de;
}

double funmu(double e, double mu, double T)
{
	return sqrt(e)/(exp((e-mu)/T)+1);
}


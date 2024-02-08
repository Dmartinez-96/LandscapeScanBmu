#ifndef MIXING_HPP
#define MIXING_HPP

#include <vector>
#include <complex>

using namespace std;

// Sfermion angles, evaluated at SUSY scale
std::vector<double> ferm_angle_calc(double& vHiggs, double& mu, double& yt, double& yc, double& yu, double& yb, double& ys, double& yd, double& ytau, double& ymu, double& ye,
									double& at, double& ac, double& au, double& ab, double& as, double& ad, double& atau, double& amu, double& ae,
									double& tanb, double& gpr, double& g2, double& mQ1_2, double& mQ2_2, double& mQ3_2,
									double& mL1_2, double& mL2_2, double& mL3_2, double& mU1_2, double& mU2_2, double& mU3_2, double& mD1_2, double& mD2_2, double& mD3_2,
									double& mE1_2, double& mE2_2, double& mE3_2);
// Higgs mixing angle, evaluated at SUSY scale
double alpha_angle_calc(const double& tanb, const double& mA0_2);
// Neutralino mixing matrix, evaluated at SUSY scale
std::vector<std::vector<complex<double>>> N_neutralino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M1, const double& M2, const double& g2, const double& gp);// ,
															//const double& mQ1sq, const double& mL1sq, const double& mAsq);
// Chargino mixing matrices, evaluated at SUSY scale
std::vector<std::vector<double>> U_chargino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M2, const double& g2);
std::vector<std::vector<double>> V_chargino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M2, const double& g2);

#endif
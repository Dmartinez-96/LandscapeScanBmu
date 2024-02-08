#include <vector>
#include <cmath>
#include <iostream>
#include "Yukawa_routine.hpp"
#include "Gauge_routine.hpp"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double B0_PV_0mom(double m1, double m2) {
	double B0_M = max(m1, m2);
	double B0_m = min(m1, m2);
	return (((-1.0) * log(pow((B0_M / 91.1876), 2.0))) + 1.0 + ((pow(B0_m, 2.0) / (pow(B0_m, 2.0) - pow(B0_M, 2.0))) * (log(pow((B0_M / B0_m), 2.0)))));
}

double B1_PV_0mom(double m1, double m2) {
	double B1_M = max(m1, m2);
	double B1_m = min(m1, m2);
	double B1_x = pow((m2 / m1), 2.0);

	double condterm;
	
	if (B1_x > 1.0) {
		condterm = log(B1_x);
	}
	else {
		condterm = 0.0;
	}
	return (0.5 * (((-1.0) * log(pow((B1_M / 91.1876), 2.0)))
				   + 0.5 + (1.0 / (1.0 - B1_x)) + (log(B1_x) / (pow((1.0 - B1_x), 2.0))) - condterm));
}

double topMass_Solver(const double& tan_beta, const double& m_gluino, const double& m_stop_1, const double& m_stop_2,
					  const double& g2_mZ, const double& g3_mZ, const double& theta_stop) {
	// XtTanb_SUSY = (at/yt)*tanb-mu might need to be checked
	// Needed parameters
	double mtpole = 173.2;
	double g2_2 = pow(g2_mZ, 2.0);
	double g3_2 = pow(g3_mZ, 2.0);
	double zeta3 = 1.2020569031595942;

	// Get mt corrections from SM and SUSY
	double Deltamt_QCD_1l = (g3_2 / (12.0 * pow(M_PI, 2.0))) * (5.0 + (3.0 * log(pow((91.1876 / 173.2), 2.0))));
	double Deltamt_QCD_2l = ((pow(g3_2, 2.0) / (4608.0 * pow(M_PI, 4.0)))
		* ((396.0 * pow(log(pow((173.2 / 91.1876), 2.0)), 2.0))
			- (1476.0 * log(pow((173.2 / 91.1876), 2.0))) - (48.0 * zeta3) + 2011.0 + ((16.0 * pow(M_PI, 2.0)) * (1.0 + log(4.0))))) - (pow(Deltamt_QCD_1l, 2.0));
	double Deltamt_stop_gl = (((-1.0) * g3_2 / (12.0 * pow(M_PI, 2.0)))) * (B1_PV_0mom(m_gluino, m_stop_1) + B1_PV_0mom(m_gluino, m_stop_2)
		- (sin(2.0 * theta_stop) * (m_gluino / 173.2) * (B0_PV_0mom(m_gluino, m_stop_1) - B0_PV_0mom(m_gluino, m_stop_2))));
	// Calculate running mt and return
	double mt_solution = mtpole / (1.0 + Deltamt_QCD_1l + Deltamt_QCD_2l + Deltamt_stop_gl);
	return mt_solution;

}

struct topYukSolParams {
	double tan_beta;
	double m_gluino;
	double m_stop_1;
	double m_stop_2;
	double g2_mZ;
	double g3_mZ;
	double theta_stop;
};

double num_top_Solver(const double& tan_beta, const double& m_gluino, const double& m_stop_1, const double& m_stop_2,
					  const double& g2_mZ, const double& g3_mZ, const double& theta_stop) {
	// XtTanb_SUSY = (at/yt)*tanb-mu might need to be checked
	// Needed parameters
	double mtpole = 173.2;
	double g2_2 = pow(g2_mZ, 2.0);
	double g3_2 = pow(g3_mZ, 2.0);
	double zeta3 = 1.2020569031595942;

	// Get mt corrections from SM and SUSY
	double Deltamt_QCD_1l = (g3_2 / (12.0 * pow(M_PI, 2.0))) * (5.0 + (3.0 * log(pow((91.1876 / 173.2), 2.0))));
	double Deltamt_QCD_2l = ((pow(g3_2, 2.0) / (4608.0 * pow(M_PI, 4.0)))
		* ((396.0 * pow(log(pow((173.2 / 91.1876), 2.0)), 2.0))
			- (1476.0 * log(pow((173.2 / 91.1876), 2.0))) - (48.0 * zeta3) + 2011.0 + ((16.0 * pow(M_PI, 2.0)) * (1.0 + log(4.0))))) - (pow(Deltamt_QCD_1l, 2.0));
	double Deltamt_stop_gl = (((-1.0) * g3_2 / (12.0 * pow(M_PI, 2.0)))) * (B1_PV_0mom(m_gluino, m_stop_1) + B1_PV_0mom(m_gluino, m_stop_2)
		- (sin(2.0 * theta_stop) * (m_gluino / 173.2) * (B0_PV_0mom(m_gluino, m_stop_1) - B0_PV_0mom(m_gluino, m_stop_2))));
	// Calculate running mt and return
	double mt_solution = mtpole / (1.0 + Deltamt_QCD_1l + Deltamt_QCD_2l + Deltamt_stop_gl);
	return mt_solution;

}

struct botYukSolParams {
	double tan_beta;
	double m_gluino;
	double m_stop_1;
	double m_stop_2;
	double m_sbot_1;
	double m_sbot_2;
	double m_charg_1;
	double m_charg_2;
	double g2_mZ;
	double g3_mZ;
	double theta_stop;
	double theta_sbot;
};

double botMass_Solver(const double& tan_beta, const double& m_gluino, const double& m_stop_1, const double& m_stop_2, const double& m_sbot_1, const double& m_sbot_2, const double& m_chargino_1, const double& m_chargino_2,
				 	  const double& g2_mZ, const double& g3_mZ, const double& theta_stop, const double& theta_sbot, const double& yt, const double& atmZ, const double& mu) {
	double g2_2 = pow(g2_mZ, 2.0);
	double g3_2 = pow(g3_mZ, 2.0);
	double mb_DRbar_SM_mZ = 2.83;
	
	// Now get mb corrections from SM and SUSY
	double Deltamb_sbot_gl = (((-1.0) * g3_2 / (12.0 * pow(M_PI, 2.0)))) * (B1_PV_0mom(m_gluino, m_sbot_1) + B1_PV_0mom(m_gluino, m_sbot_2)
																			 - (sin(2.0 * theta_sbot) * (m_gluino / 173.2) * (B0_PV_0mom(m_gluino, m_sbot_1) - B0_PV_0mom(m_gluino, m_sbot_2))));
	double Deltamb_sbot_chargino = (((g2_2 / (16.0 * pow(M_PI, 2.0)))
									 * ((m_chargino_1 * m_chargino_2 * tan_beta / (pow(m_chargino_2, 2.0) - pow(m_chargino_1, 2.0)))
										* ((pow(cos(theta_stop), 2.0) * (B0_PV_0mom(m_chargino_1, m_stop_1) - B0_PV_0mom(m_chargino_2, m_stop_1)))
										   + (pow(sin(theta_stop), 2.0) * (B0_PV_0mom(m_chargino_1, m_stop_2) - B0_PV_0mom(m_chargino_2, m_stop_2)))))));
	Deltamb_sbot_chargino += ((yt / (16.0 * pow(M_PI, 2.0))) * m_chargino_2 * (((atmZ * tan_beta) + (yt * mu)) / (pow(m_stop_1, 2.0) - pow(m_stop_2, 2.0)))
							  * (B0_PV_0mom(m_chargino_2, m_stop_1) - B0_PV_0mom(m_chargino_2, m_stop_2)));
	double mb_solution = mb_DRbar_SM_mZ * (1.0 - Deltamb_sbot_gl - Deltamb_sbot_chargino);
	return mb_solution;
}

struct tauYukSolParams {
	double tan_beta;
	double m_chargino_1;
	double m_chargino_2;
	double m_tau_sneutrino;
	double g2_mZ;
	double g3_mZ;
};

double tauMass_Solver(const double& tan_beta, const double& g1_mZ, const double& g2_mZ, const double& g3_mZ, const double& m_chargino_1, const double& m_chargino_2, const double& m_tau_sneutrino) {
	double g1_2 = pow(g1_mZ, 2.0);
	double g2_2 = pow(g2_mZ, 2.0);
	double g3_2 = pow(g3_mZ, 2.0);
	double mtau_MSbar_SM_mZ = 1.7463;
	double DRbar_conv = (3.0 / (128.0 * pow(M_PI, 2.0))) * (g1_2 - g2_2);
	double Deltamtau = (g2_2 / (16.0 * pow(M_PI, 2.0))) * (m_chargino_1 * m_chargino_2 * tan_beta / (pow(m_chargino_2, 2.0) - pow(m_chargino_1, 2.0)))\
		* (B0_PV_0mom(m_chargino_1, m_tau_sneutrino) - B0_PV_0mom(m_chargino_2, m_tau_sneutrino));
	double mtau_solution = mtau_MSbar_SM_mZ * (1.0 + Deltamtau);
	return mtau_solution;
}

std::vector<double> get_init_yukawas(const double& tan_beta, const double& g1_mZ, const double& g2_mZ, const double& g3_mZ) {
	const double mtpole = 173.2;
	double gpr_2 = 3.0 * pow(g1_mZ, 2.0) / 5.0;
	double g2_2 = pow(g2_mZ, 2.0);
	double g3_2 = pow(g3_mZ, 2.0);
	double vHiggs_mZ = sqrt(2.0) * 91.1876 / sqrt(gpr_2 + g2_2);
	double beta_mZ = atan(tan_beta);
	double vu_mZ = vHiggs_mZ * sin(beta_mZ);
	double vd_mZ = vHiggs_mZ * cos(beta_mZ);
	double zeta3 = 1.2020569031595942;

	// Get mt corrections from SM
	double Deltamt_QCD_1l = (g3_2 / (12.0 * pow(M_PI, 2.0))) * (5.0 + (3.0 * log(pow((91.1876 / 173.2), 2.0))));
	double Deltamt_QCD_2l = ((pow(g3_2, 2.0) / (4608.0 * pow(M_PI, 4.0)))
		* ((396.0 * pow(log(pow((173.2 / 91.1876), 2.0)), 2.0))
			- (1476.0 * log(pow((173.2 / 91.1876), 2.0))) - (48.0 * zeta3) + 2011.0 + ((16.0 * pow(M_PI, 2.0)) * (1.0 + log(4.0))))) - (pow(Deltamt_QCD_1l, 2.0));
	// Calculate running mt and return
	double mt_solution = mtpole / (1.0 + Deltamt_QCD_1l + Deltamt_QCD_2l);
	double mb_DRbar_SM_mZ = 2.83;
	double mtau_MSbar_SM_mZ = 1.7463;

	const double mymc = 0.619;
	const double mymu = 1.27E-3;
	const double myms = 0.055;
	const double mymd = 2.9e-3;
	const double mymmu = 0.1027181359;
	const double myme = 4.86570161e-4;


	double yt_eval = mt_solution / vu_mZ;
	double yc_eval = mymc / vu_mZ;
	double yu_eval = mymu / vu_mZ;
	double yb_eval = mb_DRbar_SM_mZ / vd_mZ;
	double ys_eval = myms / vd_mZ;
	double yd_eval = mymd / vd_mZ;
	double ytau_eval = mtau_MSbar_SM_mZ / vd_mZ;
	double ymu_eval = mymmu / vd_mZ;
	double ye_eval = myme / vd_mZ;
	vector<double> inityuks = { yt_eval, yc_eval, yu_eval, yb_eval, ys_eval, yd_eval, ytau_eval, ymu_eval, ye_eval };
	return inityuks;
}

std::vector<double> get_yukawa_couplings(const double& tan_beta, const double& m_gluino, const double& m_stop_1, const double& m_stop_2, const double& m_sbot_1, const double& m_sbot_2,
										 const double& m_tau_sneutrino, const double& m_chargino_1, const double& m_chargino_2, const double& g1_mZ, const double& g2_mZ,
										 const double& g3_mZ, const double& theta_stop, const double& theta_sbot, const double& yt, const double& atmZ, const double& mu) {
	double currmt = topMass_Solver(tan_beta, m_gluino, m_stop_1, m_stop_2, g2_mZ, g3_mZ, theta_stop);
	double currmb = botMass_Solver(tan_beta, m_gluino, m_stop_1, m_stop_2, m_sbot_1, m_sbot_2, m_chargino_1, m_chargino_2, g2_mZ, g3_mZ, theta_stop, theta_sbot, yt, atmZ, mu);
	double currmtau = tauMass_Solver(tan_beta, g1_mZ, g2_mZ, g3_mZ, m_chargino_1, m_chargino_2, m_tau_sneutrino);
	const double mymc = 0.619;
	const double mymu = 1.27E-3;
	const double myms = 0.055;
	const double mymd = 2.9e-3;
	const double mymmu = 0.1027181359;
	const double myme = 4.86570161e-4;
	double gpr_2 = 3.0 * pow(g1_mZ, 2.0) / 5.0;
	double g2_2 = pow(g2_mZ, 2.0);
	double g3_2 = pow(g3_mZ, 2.0);
	double vHiggs_mZ = sqrt(2.0) * 91.1876 / sqrt(gpr_2 + g2_2);
	double beta_mZ = atan(tan_beta);
	double vu_mZ = vHiggs_mZ * sin(beta_mZ);
	double vd_mZ = vHiggs_mZ * cos(beta_mZ);

	double yt_eval = currmt / vu_mZ;
	double yc_eval = mymc / vu_mZ;
	double yu_eval = mymu / vu_mZ;
	double yb_eval = currmb / vd_mZ;
	double ys_eval = myms / vd_mZ;
	double yd_eval = mymd / vd_mZ;
	double ytau_eval = currmtau / vd_mZ;
	double ymu_eval = mymmu / vd_mZ;
	double ye_eval = myme / vd_mZ;
	return {yt_eval, yc_eval, yu_eval, yb_eval, ys_eval, yd_eval, ytau_eval, ymu_eval, ye_eval};
}

//int main() {
//	/*std::vector<double> yukawa_solutions = get_yukawa_couplings(50.0, 2.56384017e+03, 1.79083697e+03, 2.09111236e+03, 2.01968721e+03, 2.05558868e+03, 7.99752212e+02,
//																1.32682332e+03, 2.17219938e+03, sqrt(5.0 / 3.0) * 0.353683, 0.63399, 1.11297, 1.52291, 0.0261675);
//	cout << "yt = " << yukawa_solutions[0] << endl << "yc = " << yukawa_solutions[1] << endl << "yu = " << yukawa_solutions[2] << endl << "yb = " << yukawa_solutions[3] << endl << "ys = " << yukawa_solutions[4] << endl << "yd = " << yukawa_solutions[5] << endl;
//	cout << "ytau = " << yukawa_solutions[6] << endl << "ymu = " << yukawa_solutions[7] << endl << "ye = " << yukawa_solutions[8] << endl << endl << endl;
//
//	double currmt = 173.2;
//	double currmb = 2.83;
//	double currmtau = 1.74624;
//	double gpr_2 = pow((0.353683), 2.0);
//	double g2_2 = pow(0.63399, 2.0);
//	double vHiggs_mZ = sqrt(2.0) * 91.1876 / sqrt(gpr_2 + g2_2);
//	double beta_mZ = atan(10.0);
//	double vu_mZ = vHiggs_mZ * sin(beta_mZ);
//	double vd_mZ = vHiggs_mZ * cos(beta_mZ);
//
//	double yt_eval = currmt / vu_mZ;
//	double yb_eval = currmb / vd_mZ;
//	double ytau_eval = currmtau / vd_mZ;
//	cout << "no thresh yt, yb, ytau: " << endl;
//	cout << "yt = " << yt_eval << endl << "yb = " << yb_eval << endl << "ytau = " << ytau_eval << endl;
//	return 0;*/
//	vector<double> testalphacouplings = get_gauge_couplings(2.56383979e+03, 8.38326330e+02, 1.80173518e+03, 2.04782871e+03, 2.04996947e+03, 2.15958524e+03, 2.04996996e+03, 2.15959357e+03,
//		1.89672637e+03, 2.02802923e+03, 2.03448557e+03, 2.16114659e+03, 2.03450366e+03, 2.16115491e+03, 4.71481397e+02, 8.80012830e+02,
//		4.77195390e+02, 8.05807809e+02, 4.75975476e+02, 8.05444680e+02, 1.33043414e+03, 2.14429404e+03);
//	/*vector<double> testalphacouplings = get_gauge_couplings(1.16158869e+03, 7.23808395e+02, 8.07416231e+02, 1.01619685e+03, 1.01777045e+03, 1.05672617e+03, 1.01777304e+03, 1.05672845e+03,
//															9.70036389e+02, 1.01108521e+03, 1.01418376e+03, 1.05984264e+03, 1.01418597e+03, 1.05984490e+03, 2.22503135e+02, 3.62572234e+02,
//															2.29779862e+02, 3.61388872e+02, 2.29790484e+02, 3.61392200e+02, 3.86464683e+02, 6.50867662e+02);*/
//
//															//cout << "Tree level:" << endl;
//	vector<double> treetestgaugecouplings = first_run_gauge_couplings();
//	//cout << 3.0 * pow(treetestgaugecouplings[0], 2.0) / (5.0 * 4.0 * M_PI * 80.404 / 91.1876) << endl;
//	//cout << "g'(mZ, no thresh) = " << sqrt(3.0 / 5.0) * treetestgaugecouplings[0] << endl;
//	//cout << "g2(mZ, no thresh) = " << treetestgaugecouplings[1] << endl;
//	//cout << "g3(mZ, no thresh) = " << treetestgaugecouplings[2] << endl;
//
//	const double epsabs = 1e-8;
//	const double epsrel = 1e-8;
//	const int max_iter = 100;
//	int status, iter = 0;
//	double r = 0.0;
//
//	std::vector<double> pole_higgs_masses = { 119.303387, 834.505765, 834.910073, 838.721155 };
//	std::vector<double> pole_squark_slep_masses = { 1.80173518e+03, 2.04782871e+03, 2.04996947e+03, 2.15958524e+03, 2.04996996e+03, 2.15959357e+03,
//													1.89672637e+03, 2.02802923e+03, 2.03448557e+03, 2.16114659e+03, 2.03450366e+03, 2.16115491e+03, 4.71481397e+02, 8.80012830e+02,
//													4.77195390e+02, 8.05807809e+02, 4.75975476e+02, 8.05444680e+02, 8.45072132e+02, 8.01345307e+02, 8.01272618e+02 };
//	std::vector<double> neutralino_masses = { 9.75773643e+02, 1.33027600e+03 + 800.0, -2.14038273e+03, 2.14404607e+03 };
//	std::vector<double> chargino_masses = { 1.33043414e+03, 2.14429404e+03 };
//	std::vector<double> SUSYscale_soft_trilins = { -1.95268377e+03, -9.78023426, -1.97485558e-02, -1.51587202e+03, -4.69315410e+01, -2.29779672, -7.40641273e+02, -4.83132559e+01, -2.33771633e-01 };
//	std::vector<double> yukawas = { 8.27486425e-01, 3.13097972e-03, 6.32202671e-06, 5.15044690e-01, 1.36762829e-02, 6.69561754e-04, 5.48238006e-01, 2.92755148e-02, 1.41580626e-04 };
//	std::vector<double> running_gauge_couplings = { sqrt(5.0 / 3.0) * 3.63839068e-01, 6.37022576e-01, 1.02490785e+00 };
//	std::vector<double> running_squark_slep_masses_sq = { 4.29759696e+06, 4.29742448e+06, 3.77064283e+06, 6.16299345e+05, 6.16417177e+05, 6.88472415e+05, 3.87746471e+06, 3.87744945e+06, 3.17568224e+06, 3.81567120e+06,
//														  3.81533444e+06, 3.43067759e+06, 2.12694905e+05, 2.12917324e+05, 3.51766873e+05 };
//	std::vector<double> running_gaugino_masses = { 9.80721418e+02, 1.31162910e+03, 2.50318413e+03 };
//	double tempalphaem = testalphacouplings[0];
//	std::cout << "alphaem = " << tempalphaem << endl;
//	double tempalphas = testalphacouplings[1];
//	std::cout << "alphas = " << tempalphas << endl;
//	double mu = 2.16440506e+03;
//	const double vHiggs = 174.1035847379182;
//	double tanbMZ = 50.0;
//	double tanbSUSY = 4.92411376e+01;
//
//	double a = 0.45;
//
//	double b = 0.46;
//	double intervalStart = findIntervalStart(a, b, tempalphaem, tempalphas, pole_squark_slep_masses,
//		neutralino_masses, pole_higgs_masses, chargino_masses, vHiggs, mu,
//		yukawas, SUSYscale_soft_trilins, running_squark_slep_masses_sq,
//		running_gaugino_masses, running_gauge_couplings, tanbMZ, tanbSUSY);
//
//	double x_lo = intervalStart, x_hi = intervalStart + 0.01;
//
//	GCSolParams params = { pole_higgs_masses, pole_squark_slep_masses, neutralino_masses, chargino_masses, SUSYscale_soft_trilins, yukawas, running_gauge_couplings, running_squark_slep_masses_sq,
//						   running_gaugino_masses, tempalphaem, tempalphas, vHiggs, mu, tanbMZ, tanbSUSY };
//
//	gsl_function F;
//	F.function = &get_sinthW;
//	F.params = &params;
//
//	gsl_root_fsolver* s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
//	gsl_root_fsolver_set(s, &F, x_lo, x_hi);
//
//	//std::cout << "Using GSL Brent's method to find root" << std::endl;
//	std::cout << "Iter\t[x_lo, x_hi]\t\tRoot" << std::endl;
//
//	do {
//		iter++;
//		status = gsl_root_fsolver_iterate(s);
//		r = gsl_root_fsolver_root(s);
//		x_lo = gsl_root_fsolver_x_lower(s);
//		x_hi = gsl_root_fsolver_x_upper(s);
//		status = gsl_root_test_interval(x_lo, x_hi, epsabs, epsrel);
//
//		std::cout << iter << "\t[" << x_lo << ", " << x_hi << "]\t" << r << std::endl;
//	} while (status == GSL_CONTINUE && iter < max_iter);
//
//	gsl_root_fsolver_free(s);
//
//	std::cout << "Root found at sin(thetaW) = " << r << std::endl;
//
//	std::cout << "g'(mZ) = " << sqrt(4.0 * M_PI * tempalphaem) / sqrt(1.0 - pow(r, 2.0)) << endl;
//	std::cout << "g2(mZ) = " << sqrt(4.0 * M_PI * tempalphaem) / r << endl;
//	std::cout << "g3(mZ) = " << sqrt(4.0 * M_PI * tempalphas) << endl;
//
//	double m_gluino = 2.56383979e+03;
//	double m_stop_1 = 1.80173518e+03;
//	double m_stop_2 = 2.04782871e+03;
//	double g1_mZ = sqrt(5.0 / 3.0) * sqrt(4.0 * M_PI * tempalphaem) / sqrt(1.0 - pow(r, 2.0));
//	double g2_mZ = sqrt(4.0 * M_PI * tempalphaem) / r;
//	double g3_mZ = sqrt(4.0 * M_PI * tempalphas);
//	double theta_stop = 1.06617;
//	double currmt = topMass_Solver(tanbMZ, m_gluino, m_stop_1, m_stop_2, g2_mZ, g3_mZ, theta_stop);
//	double gpr_2 = 3.0 * pow(g1_mZ, 2.0) / 5.0;
//	double g2_2 = pow(g2_mZ, 2.0);
//	double g3_2 = pow(g3_mZ, 2.0);
//	double vHiggs_mZ = sqrt(2.0) * 91.1876 / sqrt(gpr_2 + g2_2);
//	double beta_mZ = atan(tanbMZ);
//	double vu_mZ = vHiggs_mZ * sin(beta_mZ);
//	double vd_mZ = vHiggs_mZ * cos(beta_mZ);
//	double curr_yt = currmt / vu_mZ;
//	double curr_atmZ = (-2.58498115e+03) * curr_yt;
//	double curr_mumZ = 2080.2967;
//
//	std::vector<double> yukawa_solutions = get_yukawa_couplings(50.0, 2.56384017e+03, 1.79083697e+03, 2.09111236e+03, 2.01968721e+03, 2.05558868e+03, 7.99752212e+02,
//		1.32682332e+03, 2.17219938e+03, sqrt(5.0 / 3.0) * sqrt(4.0 * M_PI * tempalphaem) / sqrt(1.0 - pow(r, 2.0)), sqrt(4.0 * M_PI * tempalphaem) / r, sqrt(4.0 * M_PI * tempalphas), 1.06617, 0.499946,
//		curr_yt, curr_atmZ, curr_mumZ);
//	cout << "yt = " << yukawa_solutions[0] << endl << "yc = " << yukawa_solutions[1] << endl << "yu = " << yukawa_solutions[2] << endl << "yb = " << yukawa_solutions[3] << endl << "ys = " << yukawa_solutions[4] << endl << "yd = " << yukawa_solutions[5] << endl;
//	cout << "ytau = " << yukawa_solutions[6] << endl << "ymu = " << yukawa_solutions[7] << endl << "ye = " << yukawa_solutions[8] << endl << endl << endl;
//
//	/*double currmt = 173.2;
//	double currmb = 2.83;
//	double currmtau = 1.74624;
//	double gpr_2 = pow((0.353683), 2.0);
//	double g2_2 = pow(0.63399, 2.0);
//	double vHiggs_mZ = sqrt(2.0) * 91.1876 / sqrt(gpr_2 + g2_2);
//	double beta_mZ = atan(10.0);
//	double vu_mZ = vHiggs_mZ * sin(beta_mZ);
//	double vd_mZ = vHiggs_mZ * cos(beta_mZ);
//
//	double yt_eval = currmt / vu_mZ;
//	double yb_eval = currmb / vd_mZ;
//	double ytau_eval = currmtau / vd_mZ;
//	cout << "no thresh yt, yb, ytau: " << endl;
//	cout << "yt = " << yt_eval << endl << "yb = " << yb_eval << endl << "ytau = " << ytau_eval << endl;*/
//	return 0;
//}
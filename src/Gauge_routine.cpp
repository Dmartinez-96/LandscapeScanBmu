#include <vector>
#include <cmath>
#include <iostream>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "Gauge_routine.hpp"
//#include "pole_mass_calcs.hpp"
//#include "mixing.hpp"

using namespace std;

//#ifndef M_PIl
//#define M_PIl 3.14159265358979323846
//#endif

// double findIntervalStart(double a, double b, double& tempalphaem, double& tempalphas, std::vector<double>& pole_squark_slep_masses,
// 						 std::vector<double>& neutralino_masses, std::vector<double>& pole_higgs_masses, std::vector<double>& chargino_masses, double& vHiggs, double& mu,
// 						 std::vector<double>& yukawas, std::vector<double>& SUSYscale_soft_trilins, std::vector<double>& running_squark_slep_masses_sq,
// 						 std::vector<double>& running_gaugino_masses, std::vector<double>& running_gauge_couplings, double& tanbMZ, double& tanbSUSY) {
// 	cout << "Testing interval [" << a << ", " << b << "]:" << endl;
// 	double fa = get_thetaW_RHS(a, tempalphaem, tempalphas, pole_squark_slep_masses,
// 		neutralino_masses, pole_higgs_masses, chargino_masses, vHiggs, mu,
// 		yukawas, SUSYscale_soft_trilins, running_squark_slep_masses_sq,
// 		running_gaugino_masses, running_gauge_couplings, tanbMZ, tanbSUSY);// - (pow(a, 2.0) * (1.0 - pow(a, 2.0)));
// 	cout << "fa: " << fa << endl;
// 	double fb = get_thetaW_RHS(b, tempalphaem, tempalphas, pole_squark_slep_masses,
// 		neutralino_masses, pole_higgs_masses, chargino_masses, vHiggs, mu,
// 		yukawas, SUSYscale_soft_trilins, running_squark_slep_masses_sq,
// 		running_gaugino_masses, running_gauge_couplings, tanbMZ, tanbSUSY);// - (pow(b, 2.0) * (1.0 - pow(b, 2.0)));
// 	cout << "fb: " << fb << endl << endl;

// 	while (fa * fb > 0) { // Same sign, no root in [a, b]
// 		// Adjust 'a' and 'b' as needed, within [0, 1]
// 		// Example: move the interval to the right
// 		if (b == 0.6) {
// 			//cout << "Couldn't find it" << endl;
// 			return 0.485;
// 		}
// 		a = b;
// 		b = std::min(b + 0.01, 0.6); // Ensure it doesn't go beyond 0.53 or so
// 		cout << "Testing interval [" << a << ", " << b << "]:" << endl;
// 		fa = fb;
// 		cout << "fa: " << fa << endl;
// 		fb = get_thetaW_RHS(b, tempalphaem, tempalphas, pole_squark_slep_masses,
// 			neutralino_masses, pole_higgs_masses, chargino_masses, vHiggs, mu,
// 			yukawas, SUSYscale_soft_trilins, running_squark_slep_masses_sq,
// 			running_gaugino_masses, running_gauge_couplings, tanbMZ, tanbSUSY);// - (pow(b, 2.0) * (1.0 - pow(b, 2.0)));
// 		cout << "fb: " << fb << endl << endl;
// 	}
// 	return a; // 'a' is the start of the interval that straddles y=0
// }


// double get_sinthW(double x, void* params) {
// 	GCSolParams* p = (GCSolParams*)params;

// 	double sinthetaW_eval = get_thetaW_RHS(x, p->tempalphaem, p->tempalphas, p->pole_squark_slep_masses,
// 		p->neutralino_masses, p->pole_higgs_masses, p->chargino_masses, p->vHiggs, p->mu,
// 		p->yukawas, p->SUSYscale_soft_trilins, p->running_squark_slep_masses_sq,
// 		p->running_gaugino_masses, p->running_gauge_couplings, p->tanbMZ, p->tanbSUSY);// - (pow(x, 2.0) * (1.0 - pow(x, 2.0)));
// 	return sinthetaW_eval;
// }

std::vector<double> get_gauge_couplings(double& m_gluino, double& m_Hpm, double& m_stop_1, double& m_stop_2, double& m_scharm_1, double& m_scharm_2,
										double& m_sup_1, double& m_sup_2, double& m_sbottom_1, double& m_sbottom_2, double& m_sstrange_1, double& m_sstrange_2,
										double& m_sdown_1, double& m_sdown_2, double& m_stau_1, double& m_stau_2, double& m_smu_1, double& m_smu_2,
										double& m_selectron_1, double& m_selectron_2, double& m_chargino_1, double& m_chargino_2) {
	
	double alphaem_MSbar_MZ = 1.0 / 137.036;
	double Deltaalphaem = ((1.0 / 3.0) + (7.0 * log(80.404 / 91.1876)) - ((16.0 / 9.0) * log(173.2 / 91.1876))
						   - ((1.0 / 3.0) * log(m_Hpm / 91.1876)) - ((4.0 / 9.0) * (log(m_stop_1 / 91.1876) + log(m_stop_2 / 91.1876) + log(m_scharm_1 / 91.1876)
																					+ log(m_scharm_2 / 91.1876) + log(m_sup_1 / 91.1876) + log(m_sup_2 / 91.1876)))
						   - ((1.0 / 9.0) * (log(m_sbottom_1 / 91.1876) + log(m_sbottom_2 / 91.1876) + log(m_sstrange_1 / 91.1876) + log(m_sstrange_2 / 91.1876)
						 	 				 + log(m_sdown_1 / 91.1876) + log(m_sdown_2 / 91.1876)))
						   - ((1.0 / 3.0) * (log(m_stau_1 / 91.1876) + log(m_stau_2 / 91.1876) + log(m_smu_1 / 91.1876) + log(m_smu_2 / 91.1876) + log(m_selectron_1 / 91.1876)
											 + log(m_selectron_2 / 91.1876)))
						   - ((4.0 / 3.0) * (log(m_chargino_1 / 91.1876) + log(m_chargino_2 / 91.1876))));
	//cout << "2 alphaem_MSbar_MZ Deltaalphaem = " << 2.0 * alphaem_MSbar_MZ * Deltaalphaem << endl;
	double alphaem_DRbar_MZ;
	
	/*if (Deltaalphaem > 0) {
		alphaem_DRbar_MZ = (M_PIl + (sqrt(M_PIl) * sqrt(M_PIl - (2.0 * alphaem_MSbar_MZ * Deltaalphaem)))) / Deltaalphaem;
	}
	else {
		alphaem_DRbar_MZ = (M_PIl - (sqrt(M_PIl) * sqrt(M_PIl - (2.0 * alphaem_MSbar_MZ * Deltaalphaem)))) / Deltaalphaem;
	}*/
	alphaem_DRbar_MZ = alphaem_MSbar_MZ / (1.0 - 0.0682 - ((alphaem_MSbar_MZ / (2.0 * M_PIl)) * Deltaalphaem));
	
	
	//cout << "alphaem_DRbar_MZ = " << alphaem_DRbar_MZ << endl;
	//cout << "alphaem_DRbar_MZ(alternative) = " << alphaem_MSbar_MZ / (1.0 - ((alphaem_MSbar_MZ / (2.0 * M_PIl)) * Deltaalphaem)) << endl;

	// Strong couplings
	double alphas_MSbar_MZ = 0.1185;
	double Deltaalphas = (0.5 - (2.0 * log(173.2 / 91.1876) / 3.0) - (2.0 * log(m_gluino / 91.1876))
							- ((1.0 / 6.0) * (log(m_stop_1 / 91.1876) + log(m_stop_2 / 91.1876) + log(m_scharm_1 / 91.1876)
								+ log(m_scharm_2 / 91.1876) + log(m_sup_1 / 91.1876) + log(m_sup_2 / 91.1876)
								+ log(m_sbottom_1 / 91.1876) + log(m_sbottom_2 / 91.1876) + log(m_sstrange_1 / 91.1876) + log(m_sstrange_2 / 91.1876)
								+ log(m_sdown_1 / 91.1876) + log(m_sdown_2 / 91.1876))));
	//cout << "Deltaalphas = " << Deltaalphas << endl;
	double alphas_DRbar_MZ;
	if (Deltaalphas > 0) {
		alphas_DRbar_MZ = (M_PIl + (sqrt(M_PIl) * sqrt(M_PIl - (2.0 * alphas_MSbar_MZ * Deltaalphas)))) / Deltaalphas;
	}
	else {
		alphas_DRbar_MZ = (M_PIl - (sqrt(M_PIl) * sqrt(M_PIl - (2.0 * alphas_MSbar_MZ * Deltaalphas)))) / Deltaalphas;
	}
	//alphas_DRbar_MZ = alphas_MSbar_MZ / (1.0 - ((alphas_MSbar_MZ / (2.0 * M_PIl)) * Deltaalphas));
	std::vector<double> retalphs = { alphaem_DRbar_MZ, alphas_DRbar_MZ };
	return retalphs;
}

std::vector<double> first_run_gauge_couplings() {
	double alphaem_MSbar_MZ = (1.0 / 137.036);// / (1.0 - 0.0682 - ((1.0 / (137.036 * M_PIl)) * ((1.0 / 3.0) + (7.0 * log(80.404 / 91.1876)) - ((16.0 / 9.0) * log(173.2 / 91.1876)))));
	//cout << "alphaem_DRbar_MZ(no thresh) = " << alphaem_MSbar_MZ << endl;
	double alphas_MSbar_MZ = 0.1185;// / (1.0 - ((0.1185 / (2.0 * M_PIl)) * (0.5 - (2.0 * log(173.2 / 91.1876)))));
	double thetaW_OS = asin(0.486);//asin(sqrt(1.0 - pow((80.404 / 91.1876), 2.0)));
	double evaluated_g1_MZ = sqrt(20.0 * M_PIl * alphaem_MSbar_MZ / 3.0) / cos(thetaW_OS);
	double evaluated_g2_MZ = sqrt(4.0 * M_PIl * alphaem_MSbar_MZ) / sin(thetaW_OS);
	double evaluated_g3_MZ = sqrt(4.0 * M_PIl * alphas_MSbar_MZ);
	std::vector<double> returngauges = { evaluated_g1_MZ, evaluated_g2_MZ, evaluated_g3_MZ };
	return returngauges;
}

//int main() {
//	vector<double> testalphacouplings = get_gauge_couplings(2.56383979e+03, 8.38326330e+02, 1.80173518e+03, 2.04782871e+03, 2.04996947e+03, 2.15958524e+03, 2.04996996e+03, 2.15959357e+03,
//											   				1.89672637e+03, 2.02802923e+03, 2.03448557e+03, 2.16114659e+03, 2.03450366e+03, 2.16115491e+03, 4.71481397e+02, 8.80012830e+02,
//															4.77195390e+02, 8.05807809e+02, 4.75975476e+02, 8.05444680e+02, 1.33043414e+03, 2.14429404e+03);
//	/*vector<double> testalphacouplings = get_gauge_couplings(1.16158869e+03, 7.23808395e+02, 8.07416231e+02, 1.01619685e+03, 1.01777045e+03, 1.05672617e+03, 1.01777304e+03, 1.05672845e+03,
//															9.70036389e+02, 1.01108521e+03, 1.01418376e+03, 1.05984264e+03, 1.01418597e+03, 1.05984490e+03, 2.22503135e+02, 3.62572234e+02,
//															2.29779862e+02, 3.61388872e+02, 2.29790484e+02, 3.61392200e+02, 3.86464683e+02, 6.50867662e+02);*/
//	
//	//cout << "Tree level:" << endl;
//	vector<double> treetestgaugecouplings = first_run_gauge_couplings();
//	//cout << 3.0 * pow(treetestgaugecouplings[0], 2.0) / (5.0 * 4.0 * M_PI * 80.404 / 91.1876) << endl;
//	//cout << "g'(mZ, no thresh) = " << sqrt(3.0 / 5.0) * treetestgaugecouplings[0] << endl;
//	//cout << "g2(mZ, no thresh) = " << treetestgaugecouplings[1] << endl;
//	//cout << "g3(mZ, no thresh) = " << treetestgaugecouplings[2] << endl;
//
//	double epsabs = 1e-8;
//	double epsrel = 1e-8;
//	int max_iter = 100;
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
//	double vHiggs = 174.1035847379182;
//	double tanbMZ = 50.0;
//	double tanbSUSY = 4.92411376e+01;
//
//	double a = 0.45;
//
//	double b = 0.46;
//	double intervalStart = findIntervalStart(a, b, tempalphaem, tempalphas, pole_squark_slep_masses,
//										     neutralino_masses, pole_higgs_masses, chargino_masses, vHiggs, mu,
//										     yukawas, SUSYscale_soft_trilins, running_squark_slep_masses_sq,
//										     running_gaugino_masses, running_gauge_couplings, tanbMZ, tanbSUSY);
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
//	return 0;
//}
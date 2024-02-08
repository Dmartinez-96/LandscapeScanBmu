#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <iostream>
#include "mixing.hpp"
//#include "Pass_Velt.hpp"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Sfermion angles, evaluated at SUSY scale
std::vector<double> ferm_angle_calc(double& vHiggs, double& mu, double& yt, double& yc, double& yu, double& yb, double& ys, double& yd, double& ytau, double& ymu, double& ye,
									double& at, double& ac, double& au, double& ab, double& as, double& ad, double& atau, double& amu, double& ae,
									double& tanb, double& gpr, double& g2, double& mQ1_2, double& mQ2_2, double& mQ3_2,
									double& mL1_2, double& mL2_2, double& mL3_2, double& mU1_2, double& mU2_2, double& mU3_2, double& mD1_2, double& mD2_2, double& mD3_2,
									double& mE1_2, double& mE2_2, double& mE3_2) {
	double vu = vHiggs * sin(atan(tanb));
	double vd = vHiggs * cos(atan(tanb));
	double theta_weinb = atan(gpr / g2);
	std::vector<double> ret_angles;
	double tan2th_t = ((2.0 * vu * (at - (mu * yt / tanb))) / (mQ3_2 - mU3_2 + ((0.5 - ((4.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double tan2th_c = ((2.0 * vu * (ac - (mu * yc / tanb))) / (mQ2_2 - mU2_2 + ((0.5 - ((4.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double tan2th_u = ((2.0 * vu * (au - (mu * yu / tanb))) / (mQ1_2 - mU1_2 + ((0.5 - ((4.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double theta_2t = atan2((2.0 * vu * (at - (mu * yt / tanb))), (mQ3_2 - mU3_2 + ((0.5 - ((4.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2t = (theta_2t >= 0) ? theta_2t : theta_2t + (2.0 * M_PI);
	double theta_2c = atan2((2.0 * vu * (ac - (mu * yc / tanb))), (mQ2_2 - mU2_2 + ((0.5 - ((4.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2c = (theta_2c >= 0) ? theta_2c : theta_2c + (2.0 * M_PI);
	double theta_2u = atan2((2.0 * vu * (au - (mu * yu / tanb))), (mQ1_2 - mU1_2 + ((0.5 - ((4.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2u = (theta_2u >= 0) ? theta_2u : theta_2u + (2.0 * M_PI);
	double theta_t = theta_2t / 2.0;
	double theta_c = theta_2c / 2.0;
	double theta_u = theta_2u / 2.0;

	double tan2th_b = ((2.0 * vd * (ab - (mu * yb * tanb))) / (mQ3_2 - mD3_2 + (((-0.5) + ((2.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double tan2th_s = ((2.0 * vd * (as - (mu * ys * tanb))) / (mQ2_2 - mD2_2 + (((-0.5) + ((2.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double tan2th_d = ((2.0 * vd * (ad - (mu * yd * tanb))) / (mQ1_2 - mD1_2 + (((-0.5) + ((2.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double theta_2b = atan2((2.0 * vd * (ab - (mu * yb * tanb))), (mQ3_2 - mD3_2 + (((-0.5) + ((2.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2b = (theta_2b >= 0) ? theta_2b : theta_2b + (2.0 * M_PI);
	double theta_2s = atan2((2.0 * vd * (as - (mu * ys * tanb))), (mQ2_2 - mD2_2 + (((-0.5) + ((2.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2s = (theta_2s >= 0) ? theta_2s : theta_2s + (2.0 * M_PI);
	double theta_2d = atan2((2.0 * vd * (ad - (mu * yd * tanb))), (mQ1_2 - mD1_2 + (((-0.5) + ((2.0 / 3.0) * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2d = (theta_2d >= 0) ? theta_2d : theta_2d + (2.0 * M_PI);
	double theta_b = theta_2b / 2.0;
	double theta_s = theta_2s / 2.0;
	double theta_d = theta_2d / 2.0;

	double tan2th_tau = ((2.0 * vd * (atau - (mu * ytau * tanb))) / (mL3_2 - mE3_2 + (((-0.5) + (2.0 * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double tan2th_mu = ((2.0 * vd * (amu - (mu * ymu * tanb))) / (mL2_2 - mE2_2 + (((-0.5) + (2.0 * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double tan2th_e = ((2.0 * vd * (ae - (mu * ye * tanb))) / (mL1_2 - mE1_2 + (((-0.5) + (2.0 * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	double theta_2tau = atan2((2.0 * vd * (atau - (mu * ytau * tanb))), (mL3_2 - mE3_2 + (((-0.5) + (2.0 * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2tau = (theta_2tau >= 0) ? theta_2tau : theta_2tau + (2.0 * M_PI);
	double theta_2mu = atan2((2.0 * vd * (amu - (mu * ymu * tanb))), (mL2_2 - mE2_2 + (((-0.5) + (2.0 * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2mu = (theta_2mu >= 0) ? theta_2mu : theta_2mu + (2.0 * M_PI);
	double theta_2e = atan2((2.0 * vd * (ae - (mu * ye * tanb))), (mL1_2 - mE1_2 + (((-0.5) + (2.0 * pow(sin(theta_weinb), 2.0))) * (pow(91.1876, 2.0) * cos(2.0 * atan(tanb))))));
	theta_2e = (theta_2e >= 0) ? theta_2e : theta_2e + (2.0 * M_PI);
	double theta_tau = 0.5 * theta_2tau;
	double theta_mu = 0.5 * theta_2mu;
	double theta_e = 0.5 * theta_2e;

	/*cout << "theta_t = " << theta_t << endl << "theta_b = " << theta_b << endl;
	cout << "sin(theta_t) = " << sin(theta_t) << "\tcos(theta_t) = " << cos(theta_t) << endl;
	cout << "sin(theta_b) = " << sin(theta_b) << "\tcos(theta_b) = " << cos(theta_b) << endl;
	cout << "sin(theta_tau) = " << sin(theta_tau) << "\tcos(theta_tau) = " << cos(theta_tau) << endl;*/
	ret_angles = { theta_t, theta_c, theta_u, theta_b, theta_s, theta_d, theta_tau, theta_mu, theta_e };
	return ret_angles;
}

// Higgs mixing angle, evaluated at SUSY scale
double alpha_angle_calc(const double& tanb, const double& mA0_2) {
	/*double tan_2alpha = (((mA0_2 + pow(91.1876, 2.0)) / (mA0_2 - pow(91.1876, 2.0))) * tan(2.0 * atan(tanb)));
	double alpha_2 = ((tan_2alpha >= 0) ? atan(tan_2alpha) : )*/
	double alpha_eval = 0.5 * atan(((mA0_2 + pow(91.1876, 2.0)) / (mA0_2 - pow(91.1876, 2.0))) * tan(2.0 * atan(tanb)));
	return alpha_eval;
}

// Neutralino mixing matrix, evaluated at SUSY scale
std::vector<std::vector<complex<double>>> N_neutralino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M1, const double& M2, const double& g2, const double& gp) {//,
															//const double& mQ1sq, const double& mL1sq, const double& mAsq) {
	Eigen::MatrixXcd neutMassMat(4, 4);
	double betaval = atan(tanb);
	const double loc_thW = atan(gp / g2);
	double sW = sqrt(1.0 - pow((80.404 / 91.1876), 2.0));//sW = sin(loc_thW);// ;
	double cW = 80.404 / 91.1876; //cW = cos(loc_thW);// 
	double mZ = 91.1876;
	/*complex<double> dM1 = ((-1.0) / (16.0 * pow(M_PI, 2.0))) * (pow(gp, 2.0)) * (((11.0 * PV_B1_genp(pow(M1, 2.0), 0.0, mQ1sq)) + (9.0 * PV_B1_genp(pow(M1, 2.0), 0.0, mL1sq)))
																				 + PV_B1_genp(pow(M1, 2.0), pow(mu, 2.0), mAsq) + PV_B1_genp(pow(M1, 2.0), pow(mu, 2.0), pow(mZ, 2.0))
																				 + ((mu / M1) * sin(2.0 * betaval) * (PV_B0_genp(pow(M1, 2.0), pow(mu, 2.0), mAsq) - PV_B0_genp(pow(M1, 2.0), pow(mu, 2.0), pow(mZ, 2.0)))));
	complex<double> dM2 = ((-1.0) / (16.0 * pow(M_PI, 2.0))) * (pow(g2, 2.0)) * (((9.0 * PV_B1_genp(pow(M2, 2.0), 0.0, mQ1sq)) + (3.0 * PV_B1_genp(pow(M2, 2.0), 0.0, mL1sq)))
																				 + PV_B1_genp(pow(M2, 2.0), pow(mu, 2.0), mAsq) + PV_B1_genp(pow(M2, 2.0), pow(mu, 2.0), pow(mZ, 2.0))
																				 + ((mu / M2) * sin(2.0 * betaval) * (PV_B0_genp(pow(M2, 2.0), pow(mu, 2.0), mAsq) - PV_B0_genp(pow(M2, 2.0), pow(mu, 2.0), pow(mZ, 2.0)))));*/
	/*neutMassMat << M1 / (1.0 - dM1), 0.0, (-1.0)* cos(betaval)* sW* mZ, sin(betaval)* sW* mZ,
					0.0, M2 / (1.0 - dM2), cos(betaval)* cW* mZ, (-1.0)* sin(betaval)* cW* mZ,
					(-1.0) * cos(betaval) * sW * mZ, cos(betaval) * cW * mZ, 0.0, mu,
					sin(betaval) * sW * mZ, (-1.0) * sin(betaval) * cW * mZ, mu, 0.0;*/
	neutMassMat << M1, 0.0, (-1.0) * gp * cos(betaval) * vHiggs / sqrt(2.0), gp * sin(betaval) * vHiggs / sqrt(2.0),
					0.0, M2, g2 * cos(betaval) * vHiggs / sqrt(2.0), (-1.0) * g2 * sin(betaval) * vHiggs / sqrt(2.0),
					(-1.0) * gp * cos(betaval) * vHiggs / sqrt(2.0), g2 * vHiggs * cos(betaval) / sqrt(2.0), 0.0, mu,
					gp * sin(betaval) * vHiggs / sqrt(2.0), (-1.0) * g2 * vHiggs * sin(betaval) / sqrt(2.0), mu, 0.0;
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(neutMassMat);
	Eigen::VectorXcd eigenvalues = solver.eigenvalues();
	//cout << "mN =\n" << eigenvalues << endl;
	for (int i = 0; i < 4; ++i) {
		//cout << "mN" << i + 1 << " = " << eigenvalues[i] << endl;
	}
	Eigen::MatrixXcd eigenvectors = solver.eigenvectors();
	//cout << "Nmix =\n" << eigenvectors << endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			//cout << "mNvec[" << i + 1 << "][" << j+1 << "] = " << eigenvectors(i, j) << endl;
		}
	}
	std::vector<int> indices(eigenvalues.size());
	std::iota(indices.begin(), indices.end(), 0);

	std::sort(indices.begin(), indices.end(), [&](int i, int j) {
		return abs(eigenvalues(i)) < abs(eigenvalues(j));
		});

	Eigen::MatrixXcd sortedEigenvectors(eigenvectors.rows(), eigenvectors.cols());
	for (size_t i = 0; i < indices.size(); ++i) {
		sortedEigenvectors.col(i) = eigenvectors.col(indices[i]);
	}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			//cout << "mNvec[" << i + 1 << "][" << j + 1 << "] = " << eigenvectors(i, j) << endl;
		}
	}



	//Eigen::MatrixXcd neutMassMatSq(4, 4);
	//neutMassMatSq = neutMassMat * neutMassMat.transpose();
	//Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver2(neutMassMatSq);
	//
	//Eigen::VectorXcd sqeigenvalues = solver2.eigenvalues();
	//for (int i = 0; i < 4; ++i) {
	//	cout << "sqrt(mNsq" << i << ") = " << sqrt(sqeigenvalues[i]) << endl;
	//}
	//Eigen::MatrixXcd eigenvectors2 = solver2.eigenvectors();
	//std::vector<int> indices2(sqeigenvalues.size());
	//std::iota(indices2.begin(), indices2.end(), 0);

	//std::sort(indices2.begin(), indices2.end(), [&](int i, int j) {
	//	return abs(sqeigenvalues(i)) < abs(sqeigenvalues(j));
	//	});

	//Eigen::MatrixXcd sortedEigenvectors2(eigenvectors2.rows(), eigenvectors2.cols());
	//for (size_t i = 0; i < indices2.size(); ++i) {
	//	sortedEigenvectors2.col(i) = eigenvectors2.col(indices2[i]);
	//}
	//for (int i = 0; i < 4; ++i) {
	//	for (int j = 0; j < 4; ++j) {
	//		cout << "mNsqvec[" << i + 1 << "][" << j + 1 << "] = " << eigenvectors2(i, j) << endl;
	//	}
	//}
	///*Eigen::MatrixXcd testNmix(4, 4);
	//testNmix << 0.99572706, -0.00311724224, 0.0262588883, -0.0124562534,
	//			0.00504890103, 0.997737036, -0.0570087125, 0.0352891197,
	//			-0.0097069806, 0.0154217646, 0.70679796, 0.707180732,
	//			-0.0271051248, 0.0653701663, 0.704625519, -0.706041735;*/
	//Eigen::MatrixXcd diagMatrix = sortedEigenvectors2.transpose() * neutMassMatSq * sortedEigenvectors2;
	////Eigen::MatrixXcd diagMatrix = testNmix.transpose() * neutMassMat * testNmix;

	//std::cout << "Diagonalized neutralinosq mat:\n" << diagMatrix << std::endl;

	std::vector<std::vector<complex<double>>> ret_neut_mix_mat(sortedEigenvectors.rows());

	Eigen::MatrixXcd invSortedEVs = sortedEigenvectors.inverse();
	for (int i = 0; i < sortedEigenvectors.rows(); ++i) {
		ret_neut_mix_mat[i].resize(sortedEigenvectors.cols());
		for (int j = 0; j < sortedEigenvectors.cols(); ++j) {
			ret_neut_mix_mat[i][j] = invSortedEVs(i, j);
			//cout << "nmix[" << i << "][" << j << "] = " << ret_neut_mix_mat[i][j] << endl;
		}
	}
	return ret_neut_mix_mat;
}

// Signed neutralino masses
//std::vector<double> Neutralino_masses(const double& vHiggs, const double& mu, const double& tanb, const double& M1, const double& M2, const double& g2, const double& gp) {//,
//	//const double& mQ1sq, const double& mL1sq, const double& mAsq) {
//	Eigen::MatrixXcd neutMassMat(4, 4);
//	double betaval = atan(tanb);
//	const double loc_thW = atan(gp / g2);
//	double sW = sqrt(1.0 - pow((80.404 / 91.1876), 2.0));//sW = sin(loc_thW);// ;
//	double cW = 80.404 / 91.1876; //cW = cos(loc_thW);// 
//	double mZ = 91.1876;
//	/*complex<double> dM1 = ((-1.0) / (16.0 * pow(M_PI, 2.0))) * (pow(gp, 2.0)) * (((11.0 * PV_B1_genp(pow(M1, 2.0), 0.0, mQ1sq)) + (9.0 * PV_B1_genp(pow(M1, 2.0), 0.0, mL1sq)))
//																				 + PV_B1_genp(pow(M1, 2.0), pow(mu, 2.0), mAsq) + PV_B1_genp(pow(M1, 2.0), pow(mu, 2.0), pow(mZ, 2.0))
//																				 + ((mu / M1) * sin(2.0 * betaval) * (PV_B0_genp(pow(M1, 2.0), pow(mu, 2.0), mAsq) - PV_B0_genp(pow(M1, 2.0), pow(mu, 2.0), pow(mZ, 2.0)))));
//	complex<double> dM2 = ((-1.0) / (16.0 * pow(M_PI, 2.0))) * (pow(g2, 2.0)) * (((9.0 * PV_B1_genp(pow(M2, 2.0), 0.0, mQ1sq)) + (3.0 * PV_B1_genp(pow(M2, 2.0), 0.0, mL1sq)))
//																				 + PV_B1_genp(pow(M2, 2.0), pow(mu, 2.0), mAsq) + PV_B1_genp(pow(M2, 2.0), pow(mu, 2.0), pow(mZ, 2.0))
//																				 + ((mu / M2) * sin(2.0 * betaval) * (PV_B0_genp(pow(M2, 2.0), pow(mu, 2.0), mAsq) - PV_B0_genp(pow(M2, 2.0), pow(mu, 2.0), pow(mZ, 2.0)))));*/
//																				 /*neutMassMat << M1 / (1.0 - dM1), 0.0, (-1.0)* cos(betaval)* sW* mZ, sin(betaval)* sW* mZ,
//																								 0.0, M2 / (1.0 - dM2), cos(betaval)* cW* mZ, (-1.0)* sin(betaval)* cW* mZ,
//																								 (-1.0) * cos(betaval) * sW * mZ, cos(betaval) * cW * mZ, 0.0, mu,
//																								 sin(betaval) * sW * mZ, (-1.0) * sin(betaval) * cW * mZ, mu, 0.0;*/
//	neutMassMat << M1, 0.0, (-1.0)* gp* cos(betaval)* vHiggs / sqrt(2.0), gp* sin(betaval)* vHiggs / sqrt(2.0),
//		0.0, M2, g2* cos(betaval)* vHiggs / sqrt(2.0), (-1.0)* g2* sin(betaval)* vHiggs / sqrt(2.0),
//		(-1.0)* gp* cos(betaval)* vHiggs / sqrt(2.0), g2* vHiggs* cos(betaval) / sqrt(2.0), 0.0, mu,
//		gp* sin(betaval)* vHiggs / sqrt(2.0), (-1.0)* g2* vHiggs* sin(betaval) / sqrt(2.0), mu, 0.0;
//	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(neutMassMat);
//	Eigen::VectorXcd eigenvalues = solver.eigenvalues();
//	std::vector<double> 
//}

static std::vector<std::vector<std::vector<double>>> UV_chargino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M2, const double& g2) {
	Eigen::MatrixXd chargmat(2, 2);
	chargmat << M2, g2 * vHiggs * sin(atan(tanb)),
				g2 * vHiggs * cos(atan(tanb)), (-1.0) * mu;
	// Singular value decomposition
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(chargmat, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXd U = svd.matrixU().inverse().normalized();
	Eigen::MatrixXd V = svd.matrixV().transpose().normalized();
	//cout << "matU =\n" << U << endl << "matV =\n" << V << endl;
	//cout << "diagonal chargino mat =\n" << U * chargmat * V.transpose() << endl << endl;
	//cout << "diagonal chargino mat2 =\n" << V.transpose() * chargmat * U << endl << endl;
	//cout << "---------------------------------------------------------------------------------------" << endl;
	std::vector<std::vector<double>> matU = { { U(0, 0), U(0, 1) }, { U(1, 0), U(1, 1) } };
	std::vector<std::vector<double>> matV = { { V(0, 0), V(0, 1) }, { V(1, 0), V(1, 1) } };
	return { matU, matV };
}

// Chargino mixing matrices, evaluated at SUSY scale
std::vector<std::vector<double>> U_chargino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M2, const double& g2) {
	double betaval = atan(tanb);
	/*std::vector<std::vector<std::vector<double>>> UV_arr = UV_chargino_calc(vHiggs, mu, tanb, M2, g2);
	std::vector<std::vector<double>> ret_Umat = UV_arr[0];*/
	Eigen::MatrixXcd XXdag(2, 2);
	XXdag << pow(M2, 2.0) + (pow(g2 * vHiggs * sin(betaval), 2.0)), g2 * vHiggs * ((M2 * cos(betaval)) + (mu * sin(betaval))),
			 g2 * vHiggs * ((M2 * cos(betaval)) + (mu * sin(betaval))), pow(mu, 2.0) + (pow(g2 * vHiggs * cos(betaval), 2.0));
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(XXdag);
	Eigen::VectorXcd eigenvalues = solver.eigenvalues().transpose();
	for (int i = 0; i < 2; ++i) {
		//cout << "mcharg" << i + 1 << " = " << sqrt(eigenvalues[i]) << endl;
	}
	Eigen::MatrixXcd eigenvectors = solver.eigenvectors().normalized();
	//cout << "mC_vecs = \n" << eigenvectors << endl;
	std::vector<int> indices(eigenvalues.size());
	std::iota(indices.begin(), indices.end(), 0);

	std::sort(indices.begin(), indices.end(), [&](int i, int j) {
		return abs(eigenvalues(i)) < abs(eigenvalues(j));
		});

	Eigen::MatrixXcd sortedEigenvectors(eigenvectors.rows(), eigenvectors.cols());
	for (size_t i = 0; i < indices.size(); ++i) {
		sortedEigenvectors.col(i) = eigenvectors.col(indices[i]).normalized();
	}

	Eigen::MatrixXcd diagMatrix = sortedEigenvectors.inverse() * XXdag * sortedEigenvectors;

	//std::cout << "Diagonalized chargino with Umat:\n" << diagMatrix << std::endl;

	std::vector<std::vector<double>> ret_Umat(sortedEigenvectors.rows());

	Eigen::MatrixXcd invSortedEVs = sortedEigenvectors.inverse();
	for (int i = 0; i < sortedEigenvectors.rows(); ++i) {
		ret_Umat[i].resize(sortedEigenvectors.cols());
		for (int j = 0; j < sortedEigenvectors.cols(); ++j) {
			ret_Umat[i][j] = real(invSortedEVs(i, j));
			//cout << "Umix[" << i << "][" << j << "] = " << ret_Umat[i][j] << endl;
		}
	}

	/*double xi1 = sqrt(4.0 + pow((XXdag[1][1] - XXdag[0][0] + sqrt((4.0 * XXdag[0][1] * XXdag[1][0]) + pow((XXdag[0][0] - XXdag[1][1]), 2.0))) / XXdag[1][0], 2.0));
	double xi2 = sqrt(4.0 + pow((XXdag[0][0] - XXdag[1][1] + sqrt((4.0 * XXdag[0][1] * XXdag[1][0]) + pow((XXdag[0][0] - XXdag[1][1]), 2.0))) / XXdag[1][0], 2.0));
	double U11, U12, U21, U22;
	U11 = (1.0 / (XXdag[1][0] * xi1)) * (XXdag[0][0] - XXdag[1][1] - sqrt((4.0 * pow(XXdag[0][1], 2.0)) + pow((XXdag[0][0] - XXdag[1][1]), 2.0)));
	U12 = ((XXdag[0][0] - XXdag[1][1] + sqrt((4.0 * XXdag[0][1] * XXdag[1][0]) + pow(XXdag[0][0] - XXdag[1][1], 2.0))) / (XXdag[1][0] * xi2));
	U21 = (2.0 / xi1);
	U22 = (2.0 / xi2);*//*
	std::cout << "Umat:" << endl << U11 << "\t" << U12 << endl << U21 << "\t" << U22 << endl;
	std::vector<std::vector<double>> Umat = { {U11, U12}, {U21, U22} };*/
	return ret_Umat;
}

std::vector<std::vector<double>> V_chargino_calc(const double& vHiggs, const double& mu, const double& tanb, const double& M2, const double& g2) {
	double betaval = atan(tanb);
	/*std::vector<std::vector<std::vector<double>>> UV_arr = UV_chargino_calc(vHiggs, mu, tanb, M2, g2);
	std::vector<std::vector<double>> ret_Vmat = UV_arr[1];*/
	////double betaval = atan(tanb);
	Eigen::MatrixXcd XdagX(2, 2);
	XdagX << pow(M2, 2.0) + (pow(g2 * vHiggs * cos(betaval), 2.0)), g2 * vHiggs * ((M2 * sin(betaval)) + (mu * cos(betaval))),
			 g2 * vHiggs * ((M2 * sin(betaval)) + (mu * cos(betaval))), pow(mu, 2.0) + (pow(g2 * vHiggs * sin(betaval), 2.0));
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(XdagX);
	Eigen::VectorXcd eigenvalues = solver.eigenvalues().transpose();
	/*for (int i = 0; i < 2; ++i) {
		cout << "mcharg" << i + 1 << " = " << sqrt(eigenvalues[i]) << endl;
	}*/
	Eigen::MatrixXcd eigenvectors = solver.eigenvectors().normalized();
	//cout << "mC_vecs = \n" << eigenvectors << endl;
	std::vector<int> indices(eigenvalues.size());
	std::iota(indices.begin(), indices.end(), 0);

	std::sort(indices.begin(), indices.end(), [&](int i, int j) {
		return abs(eigenvalues(i)) < abs(eigenvalues(j));
		});

	Eigen::MatrixXcd sortedEigenvectors(eigenvectors.rows(), eigenvectors.cols());
	for (size_t i = 0; i < indices.size(); ++i) {
		sortedEigenvectors.col(i) = eigenvectors.col(indices[i]).normalized();
	}

	Eigen::MatrixXcd diagMatrix = sortedEigenvectors.inverse() * XdagX * sortedEigenvectors;

	//std::cout << "Diagonalized chargino with Vmat:\n" << diagMatrix << std::endl;

	std::vector<std::vector<double>> ret_Vmat(sortedEigenvectors.rows());

	Eigen::MatrixXcd invSortedEVs = sortedEigenvectors.inverse();
	for (int i = 0; i < sortedEigenvectors.rows(); ++i) {
		ret_Vmat[i].resize(sortedEigenvectors.cols());
		for (int j = 0; j < sortedEigenvectors.cols(); ++j) {
			ret_Vmat[i][j] = real(invSortedEVs(i, j));
			//cout << "Vmix[" << i << "][" << j << "] = " << ret_Vmat[i][j] << endl;
		}
	}
	/*double xip1 = sqrt(4.0 + pow((XdagX[1][1] - XdagX[0][0] + sqrt((4.0 * XdagX[0][1] * XdagX[1][0]) + pow((XdagX[0][0] - XdagX[1][1]), 2.0))) / XdagX[1][0], 2.0));
	double xip2 = sqrt(4.0 + pow((XdagX[0][0] - XdagX[1][1] + sqrt((4.0 * XdagX[0][1] * XdagX[1][0]) + pow((XdagX[0][0] - XdagX[1][1]), 2.0))) / XdagX[1][0], 2.0));
	double V11, V12, V21, V22;
	V11 = (1.0 / (XdagX[1][0] * xip1)) * (XdagX[0][0] - XdagX[1][1] - sqrt((4.0 * pow(XdagX[0][1], 2.0)) + pow((XdagX[0][0] - XdagX[1][1]), 2.0)));
	V12 = ((XdagX[0][0] - XdagX[1][1] + sqrt((4.0 * XdagX[0][1] * XdagX[1][0]) + pow(XdagX[0][0] - XdagX[1][1], 2.0))) / (XdagX[1][0] * xip2));
	V21 = (2.0 / xip1);
	V22 = (2.0 / xip2);
	std::cout << "Vmat:" << endl << V11 << "\t" << V12 << endl << V21 << "\t" << V22 << endl;
	std::vector<std::vector<double>> Vmat = { {V11, V12}, {V21, V22} };*/
	return ret_Vmat;
}

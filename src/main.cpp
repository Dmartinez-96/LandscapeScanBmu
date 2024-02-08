#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <random>
#include <algorithm>
#include <cstdlib>
//#include "main.hpp"
#include "constants.hpp"
#include "Gauge_routine.hpp"
#include "Yukawa_routine.hpp"
#include "MSSM_RGE_solver.hpp"
#include "MSSM_RGE_solver_with_U3Q3finder.hpp"
#include "MSSM_RGE_solver_with_stopfinder.hpp"
#include "tree_mass_calc.hpp"
#include "mixing.hpp"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::vector<double> beta_g1g2(const double& g1val, const double& g2val, const double& g3val,
                              const double& ytval, const double& ycval, const double& yuval,
                              const double& ybval, const double& ysval, const double& ydval,
                              const double& ytauval, const double& ymuval, const double& yeval) {
    const double loop_fac = 1.0 / (16.0 * pow(M_PIl, 2.0));
    const double loop_fac_sq = pow(loop_fac, 2.0);
    const std::vector<double> b_1 = { 33.0 / 5.0, 1.0, -3.0 };

    const std::vector<std::vector<double>> b_2 = {
        {199.0 / 25.0, 27.0 / 5.0, 88.0 / 5.0},
        {9.0 / 5.0, 25.0, 24.0},
        {11.0 / 5.0, 9.0, 14.0}
    };

    const std::vector<std::vector<double>> c_2 = {
        {26.0 / 5.0, 14.0 / 5.0, 18.0 / 5.0},
        {6.0, 6.0, 2.0},
        {4.0, 4.0, 0.0}
    };
    double dg1_dt_1 = b_1[0] * pow(g1val, 3.0);
    double dg2_dt_1 = b_1[1] * pow(g2val, 3.0);
    double dg1_dt_2 = (pow(g1val, 3.0)
        * ((b_2[0][0] * pow(g1val, 2.0))
            + (b_2[0][1] * pow(g2val, 2.0))
            + (b_2[0][2] * pow(g3val, 2.0))// Tr(Yu^2)
            - (c_2[0][0] * (pow(ytval, 2.0)
                + pow(ycval, 2.0)
                + pow(yuval, 2.0)))// end trace, begin Tr(Yd^2)
            - (c_2[0][1] * (pow(ybval, 2.0)
                + pow(ysval, 2.0)
                + pow(ydval, 2.0)))// end trace, begin Tr(Ye^2)
            - (c_2[0][2] * (pow(ytauval, 2.0)
                + pow(ymuval, 2.0)
                + pow(yeval, 2.0)))));// end trace

    double dg2_dt_2 = (pow(g2val, 3.0)
        * ((b_2[1][0] * pow(g1val, 2.0))
            + (b_2[1][1] * pow(g2val, 2.0))
            + (b_2[1][2] * pow(g3val, 2.0))// Tr(Yu^2)
            - (c_2[1][0] * (pow(ytval, 2.0)
                + pow(ycval, 2.0)
                + pow(yuval, 2.0)))// end trace, begin Tr(Yd^2)
            - (c_2[1][1] * (pow(ybval, 2.0)
                + pow(ysval, 2.0)
                + pow(ydval, 2.0)))// end trace, begin Tr(Ye^2)
            - (c_2[1][2] * (pow(ytauval, 2.0)
                + pow(ymuval, 2.0)
                + pow(yeval, 2.0)))));// end trace
    double dg1_dt = (1.0) * ((loop_fac * dg1_dt_1)
        + (loop_fac_sq * dg1_dt_2));

    double dg2_dt = (1.0) * ((loop_fac * dg2_dt_1)
        + (loop_fac_sq * dg2_dt_2));
    std::vector<double> g1g2_derivs = { dg1_dt, dg2_dt };
    return g1g2_derivs;
}

double PowerLawSample(const double& minX, const double& maxX, const double& unifVar, const int& powerDraw) {
    return pow(((unifVar * (pow(maxX, (powerDraw + 1.0)) - pow(minX, (powerDraw + 1.0))))) + pow(minX, (powerDraw + 1.0)), (1.0 / (powerDraw + 1.0)));
}

double tanb_tree_calc(const double& curr_mHusq, const double& curr_mHdsq, const double& curr_mu, const double& curr_b) {
    double treeTanbEval = ((1.0 / (2.0 * curr_b))
                           * (curr_mHusq + curr_mHdsq + (2.0 * pow(curr_mu, 2.0))
                           + sqrt(pow((curr_mHusq + curr_mHdsq + (2.0 * pow(curr_mu, 2.0))), 2.0)
                                  - (4.0 * pow(curr_b, 2.0)))));
    double treeTanbEv_2 = ((1.0 / (2.0 * curr_b))
                           * (curr_mHusq + curr_mHdsq + (2.0 * pow(curr_mu, 2.0))
                           - sqrt(pow((curr_mHusq + curr_mHdsq + (2.0 * pow(curr_mu, 2.0))), 2.0)
                                  - (4.0 * pow(curr_b, 2.0)))));
    //cout << "tanb(+) = " << treeTanbEval << endl;
    //cout << "tanb(-) = " << treeTanbEv_2 << endl;
    return treeTanbEval;
}

vector<double> getRandNUHM3_BCs(const int& nPower) {
    random_device rd;
    mt19937 gen(rd());
    bernoulli_distribution binD(0.5);
    uniform_real_distribution<> dis(0.0, 1.0);
    int signBinaryRand = binD(gen);
    int signA0;
    if (signBinaryRand == 0) {
        signA0 = -1;
    } else {
        signA0 = 1;
    }
    double randmHusq, randmHdsq, randm012, randm03, randmhf, randA0, randB;
    randmHusq = pow(PowerLawSample(0.1, 10.0, dis(gen), nPower) * 1000.0, 2.0);
    randmHdsq = pow(PowerLawSample(0.1, 10.0, dis(gen), nPower) * 1000.0, 2.0);
    randm012 = PowerLawSample(0.1, 10.0, dis(gen), nPower) * 1000.0;
    randm03 = PowerLawSample(0.1, 5.0, dis(gen), nPower) * 1000.0;
    randmhf = PowerLawSample(0.1, 5.0, dis(gen), nPower) * 1000.0;
    randA0 = signA0 * PowerLawSample(0.0, 15.0, dis(gen), nPower) * 1000.0;
    randB = PowerLawSample(0.0, 10.0, dis(gen), nPower) * 1000.0;
    //cout << "mHusqGUT, mHdsqGUT, m012, m03, mhf, A0, B\n(" << sqrt(randmHusq) / 1000.0 << ")^2, (" << sqrt(randmHdsq) / 1000.0 << ")^2, " << randm012 / 1000.0 << ", " << randm03 / 1000.0 << ", " << randmhf / 1000.0 << ", " << randA0 / 1000.0 << ", " << randB / 1000.0 << endl;
    vector<double> randBCs_ret = {randmHusq, randmHdsq, randm012, randm03, randmhf, randA0, randB};
    return randBCs_ret;
}

int main() {//main(int argc, char* argv[]) { // Code for point-wise tester of NUHM3 scanner

    //cout << setprecision(20);
        //uniform_int_distribution<> randfileno(1, 1000000000);
        //string filename = "Landscape_scan_" + to_string(randfileno(gen)) + ".csv";
        //ofstream outputFile(filename);
        //mutex file_mutex;

    // The order of randBCs is (in units of GeV):
    /*
        {0: mHu^2(GUT), 1: mHd^2(GUT), 2: m0(1,2), 3: m0(3), 4: mhf, 5: A0, 6: B}
    */
    double fixedmu = 200.0;
    vector<double> randBCs = getRandNUHM3_BCs(1);

        // // Higgsino mass mu is set to a common value of 200 GeV
        // vector<double> mZgauges = first_run_gauge_couplings();
        // double initGuessTanb = atof(argv[1]);
        // vector<double> mZyuks = get_init_yukawas(initGuessTanb, mZgauges[0], mZgauges[1], mZgauges[2]);

        // cout << "g1 = " << mZgauges[0] << endl << "g2 = " << mZgauges[1] << endl << "g3 = " << mZgauges[2] << endl;
        // cout << "yt = " << mZyuks[0] << endl << "yc = " << mZyuks[1] << endl << "yu = " << mZyuks[2] << endl << "yb = " << mZyuks[3] << endl << "ys = " << mZyuks[4] << endl << "yd = " << mZyuks[5] << endl;
        // cout << "ytau = " << mZyuks[6] << endl << "ymu = " << mZyuks[7] << endl << "ye = " << mZyuks[8] << endl;
        // // cout << "Input power law draw n: ";
        // // cin >> nPower;
        // //outputFile << "mHusqGUT, mHdsqGUT, m012, m03, mhf, A0, B\n";
        // //const size_t max_threads = 10;
        // //vector<future<void>> futures;
        // vector<double> dummymZBCs = {mZgauges[0], mZgauges[1], mZgauges[2], 1000.0, 1000.0, 1000.0, 200.0, mZyuks[0], mZyuks[1], mZyuks[2],
        //                              mZyuks[3], mZyuks[4], mZyuks[5], mZyuks[6], mZyuks[7], mZyuks[8], 1000.0, 100.0, 10.0, 1000.0, 100.0, 10.0,
        //                              1000.0, 100.0, 10.0, -1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 
        //                              1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 15.0};
        // vector<double> firstGUTBCs = solveODEs(dummymZBCs, log(91.1876), log(2.5e16), 0.001);
        // cout << "g1(GUT) = " << firstGUTBCs[0] << endl << "g2(GUT) = " << firstGUTBCs[1] << endl << "g3(GUT) = " << firstGUTBCs[2] << endl;
        // cout << "yt(GUT) = " << firstGUTBCs[7] << endl << "yc(GUT) = " << firstGUTBCs[8] << endl << "yu(GUT) = " << firstGUTBCs[9] << endl << "yb(GUT) = " << firstGUTBCs[10] << endl << "ys(GUT) = " << firstGUTBCs[11] << endl << "yd(GUT) = " << firstGUTBCs[12] << endl;
        // cout << "ytau(GUT) = " << firstGUTBCs[13] << endl << "ymu(GUT) = " << firstGUTBCs[14] << endl << "ye(GUT) = " << firstGUTBCs[15] << endl;
    vector<double> initGUTBCs = {firstGaugeGUTs[0], firstGaugeGUTs[1], firstGaugeGUTs[2], randBCs[4], randBCs[4], randBCs[4],
                                 fixedmu, firstYukGUTs[0], firstYukGUTs[1], firstYukGUTs[2], firstYukGUTs[3], firstYukGUTs[4], 
                                 firstYukGUTs[5], firstYukGUTs[6], firstYukGUTs[7], firstYukGUTs[8], randBCs[5] * firstYukGUTs[0],
                                 randBCs[5] * firstYukGUTs[1], randBCs[5] * firstYukGUTs[2], randBCs[5] * firstYukGUTs[3], randBCs[5] * firstYukGUTs[4], 
                                 randBCs[5] * firstYukGUTs[5], randBCs[5] * firstYukGUTs[6], randBCs[5] * firstYukGUTs[7], randBCs[5] * firstYukGUTs[8],
                                 randBCs[0], randBCs[1], pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0),
                                 pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), 
                                 pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), randBCs[6] * fixedmu, 15.0};
    // cout << "GUT values: " << endl;
    // for (double value : initGUTBCs) {
    //     cout << value << endl;
    // }
    // Minimal acceptable tanb: 2.75
    double minSafeTanb = 2.75;
    // Maximal acceptable tanb: 60.0
    double maxSafeTanb = 60.0;
    // Can update these when I get to the threshold corrections, perhaps.

    // Find approximate stop mass scale and RGE solution at that scale
    double t_target = log(500.0);
    std::vector<RGEStruct2> weaksol_init_struct;
    double QGUT = log(2.5e16);
    double MZ_OU = 91.1876;
    weaksol_init_struct = solveODEstoapproxMSUSY(initGUTBCs, QGUT, -1.0e-3, t_target);
    double QSUSY_init = weaksol_init_struct[0].SUSYscale_eval;
    std::vector<double> weaksol_init = weaksol_init_struct[0].RGEsolvec;
    // cout << "QSUSY values: " << endl;
    // for (double value : weaksol_init) {
    //     cout << value << endl;
    // }
    //cout << "QSUSY(approx) = " << exp(QSUSY_init) << endl;
    //cout << "SUSY scale RGE solutions: " << endl;
    
    // Check for CCB minima and terminate if present
    bool CCBminCheck = false;
    vector<int> CCBindexVals = {27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41};
    for (double idxvalue : CCBindexVals) {
        if (weaksol_init[idxvalue] < 0.0) {
            cout << "CCB minima on index " << idxvalue << endl;
            CCBminCheck = true;
            return 1;
        }
    }
    
    // Check at tree-level for EWSB UFB and terminate if so.
    bool firstPotlUFB_Check = (abs(2.0 * weaksol_init[42])
                               > abs(weaksol_init[25] + weaksol_init[26]
                                     + (2.0 * pow(weaksol_init[6], 2.0))));
    if (firstPotlUFB_Check) {
        cout << "Potential UFB" << endl;
        return 1;
    }
    // Check for if b and mHu^2 + mHd^2 + 2mu^2 have same sign. If not, will get negative tanb.
    bool sign_Bmu = (weaksol_init[42] > 0);
    bool sign_otherTerm = ((weaksol_init[25] + weaksol_init[26] + (2.0 * pow(weaksol_init[6], 2.0))) > 0);
    //cout << "sign(Bmu): " << sign_Bmu << endl;
    //cout << "sign_otherTerm: " << sign_otherTerm << endl;
    bool sign_bothPos = (sign_Bmu && sign_otherTerm);
    bool sign_bothNeg = (!(sign_Bmu) && !(sign_otherTerm));
    bool sign_difSigns = !(sign_bothPos || sign_bothNeg);
    if (sign_difSigns == 1) {
        cout << "tanb negative" << endl;
        return 1;
    }
    // cout << "Potential UFB? " << PotlUFB_Check << endl;
    // cout << "2b = " << (2.0 * weaksol_init[42]) << endl;
    // cout << "mHuhat^2 + mHdhat^2 = " << (weaksol_init[25] + weaksol_init[26]
    //                                      + (2.0 * pow(weaksol_init[6], 2.0))) << endl;
    double tree_tanb_init = tanb_tree_calc(weaksol_init[25], weaksol_init[26], weaksol_init[6], weaksol_init[42]); 
    if (tree_tanb_init < 0) {
        cout << "tanb negative" << endl;
        return 1;
    } else if (tree_tanb_init < minSafeTanb) {
        cout << "tanb(tree) too small, non-perturbative" << endl;
        return 1;
    } else if (tree_tanb_init > maxSafeTanb) {
        cout << "tanb(tree) too large, non-perturbative" << endl;
        return 1;
    } else {
        cout << "tanb(tree) = " << tree_tanb_init << endl;
    }

    // Now get mass spectrum with original tanb (not new one yet)
    double vHiggs_init = sqrt(2.0 * pow(91.1876, 2.0) / (pow(weaksol_init[1], 2.0) + (0.6 * pow(weaksol_init[0], 2.0))));
    //cout << "vHiggs: " << vHiggs_init << endl;
    double expQSUSY_init = exp(QSUSY_init);
    double beta_INIT = atan(15.0);
    vector<double> massesSq_init = TreeMassCalculator(expQSUSY_init, vHiggs_init, weaksol_init[6], beta_INIT, weaksol_init[7],
                                                      weaksol_init[8], weaksol_init[9], weaksol_init[10], weaksol_init[11], weaksol_init[12],
                                                      weaksol_init[13], weaksol_init[14], weaksol_init[15], weaksol_init[0], weaksol_init[1],
                                                      weaksol_init[2], weaksol_init[29], weaksol_init[28], weaksol_init[27], weaksol_init[32], 
                                                      weaksol_init[31], weaksol_init[30], weaksol_init[35], weaksol_init[34], weaksol_init[33],
                                                      weaksol_init[38], weaksol_init[37], weaksol_init[36], weaksol_init[41], weaksol_init[40],
                                                      weaksol_init[39], weaksol_init[3], weaksol_init[4], weaksol_init[5], weaksol_init[25], weaksol_init[26],
                                                      weaksol_init[16], weaksol_init[17], weaksol_init[18], weaksol_init[19], weaksol_init[20],
                                                      weaksol_init[21], weaksol_init[22], weaksol_init[23], weaksol_init[24], weaksol_init[42]);
    vector<double> masses_init = {sqrt(abs(massesSq_init[0])), sqrt(abs(massesSq_init[1])), sqrt(abs(massesSq_init[2])), sqrt(abs(massesSq_init[3])), sqrt(abs(massesSq_init[4])),
                                  sqrt(abs(massesSq_init[5])), sqrt(abs(massesSq_init[6])), sqrt(abs(massesSq_init[7])), sqrt(abs(massesSq_init[8])), sqrt(abs(massesSq_init[9])),
                                  sqrt(abs(massesSq_init[10])), sqrt(abs(massesSq_init[11])), sqrt(abs(massesSq_init[12])), sqrt(abs(massesSq_init[13])), sqrt(abs(massesSq_init[14])),
                                  sqrt(abs(massesSq_init[15])), sqrt(abs(massesSq_init[16])), sqrt(abs(massesSq_init[17])), sqrt(abs(massesSq_init[18])), sqrt(abs(massesSq_init[19])),
                                  sqrt(abs(massesSq_init[20])), sqrt(abs(massesSq_init[21])), sqrt(abs(massesSq_init[22])), sqrt(abs(massesSq_init[23])), massesSq_init[24],
                                  sqrt(abs(massesSq_init[25])), sqrt(abs(massesSq_init[26])), sqrt(abs(massesSq_init[27])), sqrt(abs(massesSq_init[28])), sqrt(abs(massesSq_init[29])),
                                  sqrt(abs(massesSq_init[30]))};
                                  /* Output order:
     {0: mst1^2, 1: mst2^2, 2: msc1^2, 3: msc2^2, 4: msu1^2, 5: msu2^2, 6: msb1^2, 7: msb2^2,
      8: mss1^2, 9: mss2^2, 10: msd1^2, 11: msd2^2, 12: mstau1^2, 13: mstau2^2, 14: msmu1^2, 15: msmu2^2,
      16: mse1^2, 17: mse2^2, 18: m_Chargino_1, 19: m_chargino_2, 20: m_neutralino_1, 21: m_neutralino_2,
      22: m_neutralino_3, 23: m_neutralino_4, 24: m_gluino, 25: mA0^2, 26: mH0^2, 27: mHpm^2,
      28: m_snu_tau^2, 29: m_snu_mu^2, 30: m_snu_e^2}
    */
    cout << "Tree masses:" << endl;
    int printIdx = 0;
    for (double value : masses_init) {
        cout << printIdx << ": " << value << endl;
        printIdx++;
    }

    // Evolve back to mZ^OU=91.1876 GeV
    weaksol_init[6] = fixedmu;
    double tanb_init = 15.0;
    vector<double> starting_mZBCs = solveODEs(weaksol_init, QSUSY_init, log(91.1876), -0.001);
    
    vector<double> starting_mZalphas = get_gauge_couplings(masses_init[24], masses_init[27], masses_init[0], masses_init[1], masses_init[2], masses_init[3],
                                                           masses_init[4], masses_init[5], masses_init[6], masses_init[7], masses_init[8], masses_init[9], masses_init[10],
                                                           masses_init[11], masses_init[12], masses_init[13], masses_init[14], masses_init[15], masses_init[16], masses_init[17],
                                                           masses_init[18], masses_init[19]);
    vector<double> starting_gauges = {sqrt(20.0 * M_PI * starting_mZalphas[0] / 3.0) / cos(asin(0.486)), sqrt(4 * M_PI * starting_mZalphas[0]) / 0.486, sqrt(4 * M_PI * starting_mZalphas[1])};
    //double starting_vHiggs = sqrt(2.0 * pow(91.1876, 2.0) / (pow(starting_gauges[1], 2.0) + (0.6 * pow(starting_gauges[0], 2.0))));
    double starting_gpr = sqrt(3.0 / 5.0) * starting_gauges[0];
    vector<double> sferm_mixing_angles_init = ferm_angle_calc(vHiggs_init, weaksol_init[6], weaksol_init[7], weaksol_init[8], weaksol_init[9], weaksol_init[10],
                                                              weaksol_init[11], weaksol_init[12], weaksol_init[13], weaksol_init[14], weaksol_init[15], weaksol_init[16],
                                                              weaksol_init[17], weaksol_init[18], weaksol_init[19], weaksol_init[20], weaksol_init[21], weaksol_init[22],
                                                              weaksol_init[23], weaksol_init[24], tanb_init, starting_gpr, weaksol_init[1],
                                                              weaksol_init[27], weaksol_init[28], weaksol_init[29], weaksol_init[30], weaksol_init[31], weaksol_init[32],
                                                              weaksol_init[33], weaksol_init[34], weaksol_init[35], weaksol_init[36], weaksol_init[37], weaksol_init[38],
                                                              weaksol_init[39], weaksol_init[40], weaksol_init[41]);
    vector<double> starting_yuks = get_yukawa_couplings(tree_tanb_init, masses_init[24], masses_init[0], masses_init[1], masses_init[6], masses_init[7],
                                                        masses_init[28], masses_init[18], masses_init[19], starting_gauges[0], starting_gauges[1], starting_gauges[2],
                                                        sferm_mixing_angles_init[0], sferm_mixing_angles_init[3], starting_mZBCs[7], starting_mZBCs[16], starting_mZBCs[6]);

    cout << "g1(MZ) = " << starting_gauges[0] << endl << "g2(MZ) = " << starting_gauges[1] << endl << "g3(MZ) = " << starting_gauges[2] << endl;
    cout << "yt(MZ) = " << starting_yuks[0] << endl << "yc(MZ) = " << starting_yuks[1] << endl << "yu(MZ) = " << starting_yuks[2] << endl << "yb(MZ) = " << starting_yuks[3] << endl << "ys(MZ) = " << starting_yuks[4] << endl << "yd(MZ) = " << starting_yuks[5] << endl;
    cout << "ytau(MZ) = " << starting_yuks[6] << endl << "ymu(MZ) = " << starting_yuks[7] << endl << "ye(MZ) = " << starting_yuks[8] << endl;
    
    // Set up first threshold corrected set of mZ BCs, evolve to GUT scale
    vector<double> starting_threshCorr_mZBCs = {starting_gauges[0], starting_gauges[1], starting_gauges[2], 1000.0, 1000.0, 1000.0, fixedmu, starting_yuks[0], starting_yuks[1], starting_yuks[2],
                                                starting_yuks[3], starting_yuks[4], starting_yuks[5], starting_yuks[6], starting_yuks[7], starting_yuks[8], 1000.0, 100.0, 10.0, 1000.0, 100.0, 10.0,
                                                1000.0, 100.0, 10.0, -1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 
                                                1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, tree_tanb_init};
    vector<double> starting_GUTBCs = solveODEs(starting_threshCorr_mZBCs, log(91.1876), QGUT, 0.0001);
    starting_GUTBCs = {starting_GUTBCs[0], starting_GUTBCs[1], starting_GUTBCs[2], randBCs[4], randBCs[4], randBCs[4],
                       starting_GUTBCs[6], starting_GUTBCs[7], starting_GUTBCs[8], starting_GUTBCs[9], starting_GUTBCs[10], starting_GUTBCs[11], 
                       starting_GUTBCs[12], starting_GUTBCs[13], starting_GUTBCs[14], starting_GUTBCs[15], randBCs[5] * starting_GUTBCs[7],
                       randBCs[5] * starting_GUTBCs[8], randBCs[5] * starting_GUTBCs[9], randBCs[5] * starting_GUTBCs[10], randBCs[5] * starting_GUTBCs[11], 
                       randBCs[5] * starting_GUTBCs[12], randBCs[5] * starting_GUTBCs[13], randBCs[5] * starting_GUTBCs[14], randBCs[5] * starting_GUTBCs[15],
                       randBCs[0], randBCs[1], pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0),
                       pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), 
                       pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), randBCs[6] * starting_GUTBCs[6], tree_tanb_init};
    cout << "GUT values after thresholds: " << endl;
    for (double value : starting_GUTBCs) {
        cout << value << endl;
    }
    vector<double> betag1g2GUT = beta_g1g2(starting_GUTBCs[0], starting_GUTBCs[1], starting_GUTBCs[2], starting_GUTBCs[7], starting_GUTBCs[8], starting_GUTBCs[9], starting_GUTBCs[10], starting_GUTBCs[11],
                                           starting_GUTBCs[12], starting_GUTBCs[13], starting_GUTBCs[14], starting_GUTBCs[15]);
    double newQGUT = log(exp(QGUT) * exp((starting_GUTBCs[1] - starting_GUTBCs[0]) / (betag1g2GUT[0] - betag1g2GUT[1])));
    double newgGUT = sqrt(starting_GUTBCs[0] * starting_GUTBCs[1]);
    starting_GUTBCs[0] = newgGUT;
    starting_GUTBCs[1] = newgGUT;
    cout << "New Q(GUT): " << exp(newQGUT) << endl;
    double QGUTcurr = newQGUT;
    // Now iterate until tree level tanb is converged...
    double curr_lsq = pow((15.0 - tree_tanb_init), 2.0);
    double goal_lsq = 1.0e-6;
    int currIter = 0;
    double currTreeTanb = tree_tanb_init;
    double newTreeTanb = 0.0;
    while ((curr_lsq > goal_lsq) && (currIter < 100)) {
        vector<RGEStruct> weaksol_curr_struct = solveODEstoMSUSY(starting_GUTBCs, QGUTcurr, -1.0e-3, t_target, MZ_OU);
        double QSUSY_curr = weaksol_curr_struct[0].SUSYscale_eval;
        vector<double> weaksol_curr = weaksol_curr_struct[0].RGEsolvec;
        for (double idxvalue : CCBindexVals) {
        if (weaksol_curr[idxvalue] < 0.0) {
            cout << "CCB minima in convergence loop, index " << idxvalue << endl;
            CCBminCheck = true;
            return 1;
            }
        }

        bool currPotlUFB_Check = (abs(2.0 * weaksol_curr[42])
                                  > abs(weaksol_curr[25] + weaksol_curr[26]
                                        + (2.0 * pow(weaksol_curr[6], 2.0))));
        if (currPotlUFB_Check) {
            cout << "Potential UFB in convergence loop" << endl;
            return 1;
        }
        bool sign_Bmu = (weaksol_curr[42] > 0);
        bool sign_otherTerm = ((weaksol_curr[25] + weaksol_curr[26] + (2.0 * pow(weaksol_curr[6], 2.0))) > 0);
        //cout << "sign(Bmu): " << sign_Bmu << endl;
        //cout << "sign_otherTerm: " << sign_otherTerm << endl;
        bool sign_bothPos = (sign_Bmu && sign_otherTerm);
        bool sign_bothNeg = (!(sign_Bmu) && !(sign_otherTerm));
        bool sign_difSigns = !(sign_bothPos || sign_bothNeg);
        if (sign_difSigns == 1) {
            cout << "tanb negative in convergence loop" << endl;
            return 1;
        }
        ////////////////////////////////////////////////////////////////////////////////
        newTreeTanb = tanb_tree_calc(weaksol_curr[25], weaksol_curr[26], weaksol_curr[6], weaksol_curr[42]); 
        if (newTreeTanb < 0) {
            cout << "tanb negative" << endl;
            return 1;
        } else if (newTreeTanb < minSafeTanb) {
            cout << "tanb(tree) too small, non-perturbative" << endl;
            return 1;
        } else if (newTreeTanb > maxSafeTanb) {
            cout << "tanb(tree) too large, non-perturbative" << endl;
            return 1;
        } else {
            cout << "tanb(tree) = " << newTreeTanb << endl;
        }
        // Now get mass spectrum with original tanb (not new one yet)
        double vHiggs_curr = sqrt(2.0 * pow(91.1876, 2.0) / (pow(weaksol_curr[1], 2.0) + (0.6 * pow(weaksol_curr[0], 2.0))));
        //cout << "vHiggs: " << vHiggs_curr << endl;
        double expQSUSY_curr = exp(QSUSY_curr);
        double beta_TREE = atan(currTreeTanb);
        vector<double> massesSq_curr = TreeMassCalculator(expQSUSY_curr, vHiggs_curr, weaksol_curr[6], beta_TREE, weaksol_curr[7],
                                                        weaksol_curr[8], weaksol_curr[9], weaksol_curr[10], weaksol_curr[11], weaksol_curr[12],
                                                        weaksol_curr[13], weaksol_curr[14], weaksol_curr[15], weaksol_curr[0], weaksol_curr[1],
                                                        weaksol_curr[2], weaksol_curr[29], weaksol_curr[28], weaksol_curr[27], weaksol_curr[32], 
                                                        weaksol_curr[31], weaksol_curr[30], weaksol_curr[35], weaksol_curr[34], weaksol_curr[33],
                                                        weaksol_curr[38], weaksol_curr[37], weaksol_curr[36], weaksol_curr[41], weaksol_curr[40],
                                                        weaksol_curr[39], weaksol_curr[3], weaksol_curr[4], weaksol_curr[5], weaksol_curr[25], weaksol_curr[26],
                                                        weaksol_curr[16], weaksol_curr[17], weaksol_curr[18], weaksol_curr[19], weaksol_curr[20],
                                                        weaksol_curr[21], weaksol_curr[22], weaksol_curr[23], weaksol_curr[24], weaksol_curr[42]);
        vector<double> masses_curr = {sqrt(abs(massesSq_curr[0])), sqrt(abs(massesSq_curr[1])), sqrt(abs(massesSq_curr[2])), sqrt(abs(massesSq_curr[3])), sqrt(abs(massesSq_curr[4])),
                                      sqrt(abs(massesSq_curr[5])), sqrt(abs(massesSq_curr[6])), sqrt(abs(massesSq_curr[7])), sqrt(abs(massesSq_curr[8])), sqrt(abs(massesSq_curr[9])),
                                      sqrt(abs(massesSq_curr[10])), sqrt(abs(massesSq_curr[11])), sqrt(abs(massesSq_curr[12])), sqrt(abs(massesSq_curr[13])), sqrt(abs(massesSq_curr[14])),
                                      sqrt(abs(massesSq_curr[15])), sqrt(abs(massesSq_curr[16])), sqrt(abs(massesSq_curr[17])), sqrt(abs(massesSq_curr[18])), sqrt(abs(massesSq_curr[19])),
                                      sqrt(abs(massesSq_curr[20])), sqrt(abs(massesSq_curr[21])), sqrt(abs(massesSq_curr[22])), sqrt(abs(massesSq_curr[23])), massesSq_curr[24],
                                      sqrt(abs(massesSq_curr[25])), sqrt(abs(massesSq_curr[26])), sqrt(abs(massesSq_curr[27])), sqrt(abs(massesSq_curr[28])), sqrt(abs(massesSq_curr[29])),
                                      sqrt(abs(massesSq_curr[30]))};
                                    /* Output order:
        {0: mst1^2, 1: mst2^2, 2: msc1^2, 3: msc2^2, 4: msu1^2, 5: msu2^2, 6: msb1^2, 7: msb2^2,
        8: mss1^2, 9: mss2^2, 10: msd1^2, 11: msd2^2, 12: mstau1^2, 13: mstau2^2, 14: msmu1^2, 15: msmu2^2,
        16: mse1^2, 17: mse2^2, 18: m_Chargino_1, 19: m_chargino_2, 20: m_neutralino_1, 21: m_neutralino_2,
        22: m_neutralino_3, 23: m_neutralino_4, 24: m_gluino, 25: mA0^2, 26: mH0^2, 27: mHpm^2,
        28: m_snu_tau^2, 29: m_snu_mu^2, 30: m_snu_e^2}
        */
        // cout << "Tree masses:" << endl;
        // int printIdx = 0;
        // for (double value : masses_curr) {
        //     cout << printIdx << ": " << value << endl;
        //     printIdx++;
        // }

        // Evolve back to mZ^OU=91.1876 GeV
        weaksol_curr[6] = fixedmu;
        double tanb_curr = 15.0;
        vector<double> starting_mZBCs = solveODEs(weaksol_curr, QSUSY_curr, log(91.1876), -0.001);
        
        vector<double> starting_mZalphas = get_gauge_couplings(masses_curr[24], masses_curr[27], masses_curr[0], masses_curr[1], masses_curr[2], masses_curr[3],
                                                            masses_curr[4], masses_curr[5], masses_curr[6], masses_curr[7], masses_curr[8], masses_curr[9], masses_curr[10],
                                                            masses_curr[11], masses_curr[12], masses_curr[13], masses_curr[14], masses_curr[15], masses_curr[16], masses_curr[17],
                                                            masses_curr[18], masses_curr[19]);
        vector<double> starting_gauges = {sqrt(20.0 * M_PI * starting_mZalphas[0] / 3.0) / cos(asin(0.486)), sqrt(4 * M_PI * starting_mZalphas[0]) / 0.486, sqrt(4 * M_PI * starting_mZalphas[1])};
        //double starting_vHiggs = sqrt(2.0 * pow(91.1876, 2.0) / (pow(starting_gauges[1], 2.0) + (0.6 * pow(starting_gauges[0], 2.0))));
        double starting_gpr = sqrt(3.0 / 5.0) * starting_gauges[0];
        vector<double> sferm_mixing_angles_curr = ferm_angle_calc(vHiggs_curr, weaksol_curr[6], weaksol_curr[7], weaksol_curr[8], weaksol_curr[9], weaksol_curr[10],
                                                                weaksol_curr[11], weaksol_curr[12], weaksol_curr[13], weaksol_curr[14], weaksol_curr[15], weaksol_curr[16],
                                                                weaksol_curr[17], weaksol_curr[18], weaksol_curr[19], weaksol_curr[20], weaksol_curr[21], weaksol_curr[22],
                                                                weaksol_curr[23], weaksol_curr[24], tanb_curr, starting_gpr, weaksol_curr[1],
                                                                weaksol_curr[27], weaksol_curr[28], weaksol_curr[29], weaksol_curr[30], weaksol_curr[31], weaksol_curr[32],
                                                                weaksol_curr[33], weaksol_curr[34], weaksol_curr[35], weaksol_curr[36], weaksol_curr[37], weaksol_curr[38],
                                                                weaksol_curr[39], weaksol_curr[40], weaksol_curr[41]);
        vector<double> starting_yuks = get_yukawa_couplings(currTreeTanb, masses_curr[24], masses_curr[0], masses_curr[1], masses_curr[6], masses_curr[7],
                                                            masses_curr[28], masses_curr[18], masses_curr[19], starting_gauges[0], starting_gauges[1], starting_gauges[2],
                                                            sferm_mixing_angles_curr[0], sferm_mixing_angles_curr[3], starting_mZBCs[7], starting_mZBCs[16], starting_mZBCs[6]);

        cout << "g1(MZ) = " << starting_gauges[0] << endl << "g2(MZ) = " << starting_gauges[1] << endl << "g3(MZ) = " << starting_gauges[2] << endl;
        cout << "yt(MZ) = " << starting_yuks[0] << endl << "yc(MZ) = " << starting_yuks[1] << endl << "yu(MZ) = " << starting_yuks[2] << endl << "yb(MZ) = " << starting_yuks[3] << endl << "ys(MZ) = " << starting_yuks[4] << endl << "yd(MZ) = " << starting_yuks[5] << endl;
        cout << "ytau(MZ) = " << starting_yuks[6] << endl << "ymu(MZ) = " << starting_yuks[7] << endl << "ye(MZ) = " << starting_yuks[8] << endl;
        
        // Set up first threshold corrected set of mZ BCs, evolve to GUT scale
        vector<double> starting_threshCorr_mZBCs = {starting_gauges[0], starting_gauges[1], starting_gauges[2], 1000.0, 1000.0, 1000.0, fixedmu, starting_yuks[0], starting_yuks[1], starting_yuks[2],
                                                    starting_yuks[3], starting_yuks[4], starting_yuks[5], starting_yuks[6], starting_yuks[7], starting_yuks[8], 1000.0, 100.0, 10.0, 1000.0, 100.0, 10.0,
                                                    1000.0, 100.0, 10.0, -1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 
                                                    1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, currTreeTanb};
        vector<double> starting_GUTBCs = solveODEs(starting_threshCorr_mZBCs, log(91.1876), QGUTcurr, 0.0001);
        starting_GUTBCs = {starting_GUTBCs[0], starting_GUTBCs[1], starting_GUTBCs[2], randBCs[4], randBCs[4], randBCs[4],
                        starting_GUTBCs[6], starting_GUTBCs[7], starting_GUTBCs[8], starting_GUTBCs[9], starting_GUTBCs[10], starting_GUTBCs[11], 
                        starting_GUTBCs[12], starting_GUTBCs[13], starting_GUTBCs[14], starting_GUTBCs[15], randBCs[5] * starting_GUTBCs[7],
                        randBCs[5] * starting_GUTBCs[8], randBCs[5] * starting_GUTBCs[9], randBCs[5] * starting_GUTBCs[10], randBCs[5] * starting_GUTBCs[11], 
                        randBCs[5] * starting_GUTBCs[12], randBCs[5] * starting_GUTBCs[13], randBCs[5] * starting_GUTBCs[14], randBCs[5] * starting_GUTBCs[15],
                        randBCs[0], randBCs[1], pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0),
                        pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), 
                        pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), pow(randBCs[2], 2.0), pow(randBCs[2], 2.0), pow(randBCs[3], 2.0), randBCs[6] * starting_GUTBCs[6], currTreeTanb};
        // cout << "GUT values after thresholds: " << endl;
        // for (double value : starting_GUTBCs) {
        //     cout << value << endl;
        // }
        vector<double> betag1g2GUT = beta_g1g2(starting_GUTBCs[0], starting_GUTBCs[1], starting_GUTBCs[2], starting_GUTBCs[7], starting_GUTBCs[8], starting_GUTBCs[9], starting_GUTBCs[10], starting_GUTBCs[11],
                                            starting_GUTBCs[12], starting_GUTBCs[13], starting_GUTBCs[14], starting_GUTBCs[15]);
        newQGUT = log(exp(QGUTcurr) * exp((starting_GUTBCs[1] - starting_GUTBCs[0]) / (betag1g2GUT[0] - betag1g2GUT[1])));
        newgGUT = sqrt(starting_GUTBCs[0] * starting_GUTBCs[1]);
        starting_GUTBCs[0] = newgGUT;
        starting_GUTBCs[1] = newgGUT;
        cout << "New Q(GUT): " << exp(newQGUT) << endl;
        /////////////////////////////////////////////////////////////////////
        curr_lsq = pow((newTreeTanb - (currTreeTanb)), 2.0) + pow((newQGUT - QGUTcurr), 2.0); // could add in convergence on QSUSY
        cout << "Current least squares: " << curr_lsq << endl;
        QGUTcurr = newQGUT;
        currTreeTanb = newTreeTanb;
        currIter++;
        cout << "------------------------------------------------------------------" << endl;
    }
    if (curr_lsq < goal_lsq) {
        cout << "Successful convergence!" << endl;
    } else {
        cout << "Convergence failed." << endl;
        cout << "Reason: ran out of iteration attempts. Final lsq: " << curr_lsq << endl;
    }
    return 0;
}
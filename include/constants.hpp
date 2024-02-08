// CONSTANTS_HPP

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <vector>
#include <cmath>

// Define important constants
const double loop_fac = 1.0 / (16.0 * std::pow(M_PIl, 2.0));
const double loop_fac_sq = std::pow(loop_fac, 2.0);


// Define constant arrays
const std::vector<double> b_1l = { 33.0 / 5.0, 1.0, -3.0 };
const std::vector<double> firstYukGUTs = {0.48969085021845482064,
                                          0.0013118849863529838048,
                                          2.6915776127101603488e-6,
                                          0.076635069382742310173,
                                          0.0013303490408246554898,
                                          7.0145307635106774195e-5,
                                          0.10178079098510167677,
                                          0.0059282213267479759097,
                                          2.8080709208362459667e-05};
const std::vector<double> firstGaugeGUTs = {0.67832304347034433345,
                                            0.69366133058521628474,
                                            0.72618655795720066237};

const std::vector<std::vector<double>> b_2l = {
    {199.0 / 25.0, 27.0 / 5.0, 88.0 / 5.0},
    {9.0 / 5.0, 25.0, 24.0},
    {11.0 / 5.0, 9.0, 14.0}
};

const std::vector<std::vector<double>> c_2l = {
    {26.0 / 5.0, 14.0 / 5.0, 18.0 / 5.0},
    {6.0, 6.0, 2.0},
    {4.0, 4.0, 0.0}
};

#endif

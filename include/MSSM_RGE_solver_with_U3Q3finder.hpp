// SC_MSSM_RGE_solver_with_U3Q3finder.hpp

#ifndef SC_MSSM_U3Q3_RGE_SOLVER_HPP
#define SC_MSSM_U3Q3_RGE_SOLVER_HPP

#include <vector>
#include <cmath>

struct RGEStruct2 {
    std::vector<double> RGEsolvec;
    double SUSYscale_eval;

    RGEStruct2(const std::vector<double>& RGEsolvec, double SUSYscale_eval) : RGEsolvec(RGEsolvec), SUSYscale_eval(SUSYscale_eval) {}
};
void SCMSSM_approx_RGESolver(const std::vector<double>& x, std::vector<double>& dxdt, const double t);
//struct MyObserver2;

//class RGEapprox {};
    

std::vector<RGEStruct2> solveODEstoapproxMSUSY(std::vector<double> initialConditions, double startTime, double timeStep, double& t_target);

#endif // SC_MSSM_RGE_SOLVER_HPP
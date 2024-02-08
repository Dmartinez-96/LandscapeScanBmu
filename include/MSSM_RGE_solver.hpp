// SC_MSSM_RGE_solver.hpp

#ifndef SC_MSSM_RGE_SOLVER_HPP
#define SC_MSSM_RGE_SOLVER_HPP

#include <vector>

// struct SCMSSMRGESolver {
//     void operator()(const std::vector<double> &x, std::vector<double> &dxdt, const double t) const;
// };

void SCMSSMRGESolver(const std::vector<double>& x, std::vector<double>& dxdt, const double t);

std::vector<double> solveODEs(std::vector<double> initialConditions, double startTime, double endTime, double timeStep);

#endif // SC_MSSM_RGE_SOLVER_HPP
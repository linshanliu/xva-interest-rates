#ifndef EXPOSURE_HPP
#define EXPOSURE_HPP

#include <cstddef>
#include <vector>

#include "MtMCube.hpp"
#include "NettingSet.hpp"
#include "ScenarioSet.hpp"

// Exposure profile of one netting set on the simulation grid.
//
// With V(p,k) = sum_{i in NS} q_i V_i(p,k) the netted MtM on path p at t_k and
// D(p,k) = D(0,t_k) the stochastic discount factor:
//   EE[k]   = mean_p  D(p,k) * max(V, 0)      discounted expected (positive) exposure -> CVA
//   ENE[k]  = mean_p  D(p,k) * min(V, 0)      discounted expected negative exposure (<= 0) -> DVA
//   EEu[k]  = mean_p  max(V, 0)               undiscounted EE (for reporting / plots)
//   PFE[k]  = q-quantile_p of max(V, 0)       undiscounted potential future exposure
// Identity: EE[k] + ENE[k] = mean_p D(p,k) V(p,k).
struct ExposureProfile {
    std::vector<double> times;
    std::vector<double> EE;
    std::vector<double> ENE;
    std::vector<double> EEu;
    std::vector<double> PFE;
    double pfeQuantile = 0.95;
};

// Netted MtM of a netting set, laid out [path][k]: values[p * nSteps + k].
// Kept public: exposure allocation and incremental CVA work from these path values.
std::vector<double> nettedValues(const MtMCube& cube, const NettingSet& ns);

// Exposure statistics of a netting set. Throws std::invalid_argument on inconsistent inputs.
ExposureProfile computeExposure(const MtMCube& cube,
    const ScenarioSet& scen,
    const NettingSet& ns,
    double pfeQuantile = 0.95);

#endif // !EXPOSURE_HPP
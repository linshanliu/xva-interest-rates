#ifndef XVA_HPP
#define XVA_HPP

#include <vector>

#include "CreditCurve.hpp"
#include "Exposure.hpp"

// Unilateral CVA / DVA from a netting-set exposure profile.
//
// Assumptions:
//   - default times independent of rates (no wrong-way risk)
//   - DVA ignores first-to-default (each side's default treated separately)
//   - EE / ENE in the profile are already discounted: EE(t) = E[D(0,t) V(t)^+]
//
// Discretisation on the profile grid 0 = t_0 < t_1 < ... < t_n:
//   CVA = LGD_C * sum_k  EE*_k  * [ S_C(t_{k-1}) - S_C(t_k) ]
//   DVA = LGD_B * sum_k (-ENE*_k) * [ S_B(t_{k-1}) - S_B(t_k) ]
// where X*_k is the exposure attributed to bucket (t_{k-1}, t_k] (see BucketRule).
//
// Sign convention: CVA >= 0 is a cost, DVA >= 0 is a benefit,
// bilateral adjustment BCVA = CVA - DVA (amount to subtract from the risk-free value).

enum class BucketRule {
    Left,      // X*_k = X(t_{k-1})
    Right,     // X*_k = X(t_k)
    Average    // X*_k = (X(t_{k-1}) + X(t_k)) / 2   (trapezoid; default)
};

struct XVAResult {
    double cva = 0.0;
    double dva = 0.0;
    double bcva = 0.0;                 // cva - dva
    std::vector<double> cvaByBucket;   // contribution of (t_{k-1}, t_k], index k = 1..n (entry 0 is 0)
    std::vector<double> dvaByBucket;
};

double computeCVA(const ExposureProfile& e, const CreditCurve& counterparty, double lgd,
    BucketRule rule = BucketRule::Average);

double computeDVA(const ExposureProfile& e, const CreditCurve& own, double lgd,
    BucketRule rule = BucketRule::Average);

XVAResult computeXVA(const ExposureProfile& e,
    const CreditCurve& counterparty, double lgdCounterparty,
    const CreditCurve& own, double lgdOwn,
    BucketRule rule = BucketRule::Average);

#endif // !XVA_HPP
#include "XVA.hpp"

#include <stdexcept>

namespace {

    void validate(const ExposureProfile& e, double lgd) {
        if (!(lgd >= 0.0 && lgd <= 1.0))
            throw std::invalid_argument("XVA: LGD must be in [0, 1]");
        const std::size_t n = e.times.size();
        if (n < 2) throw std::invalid_argument("XVA: profile needs at least two times");
        if (e.EE.size() != n || e.ENE.size() != n)
            throw std::invalid_argument("XVA: profile vectors have inconsistent sizes");
        if (e.times.front() != 0.0) throw std::invalid_argument("XVA: profile must start at t = 0");
        for (std::size_t k = 1; k < n; ++k)
            if (e.times[k] <= e.times[k - 1])
                throw std::invalid_argument("XVA: profile times must be strictly increasing");
    }

    double bucketValue(const std::vector<double>& x, std::size_t k, BucketRule rule) {
        switch (rule) {
        case BucketRule::Left:  return x[k - 1];
        case BucketRule::Right: return x[k];
        default:                return 0.5 * (x[k - 1] + x[k]);
        }
    }

    // lgd * sum_k sign * X*_k * [S(t_{k-1}) - S(t_k)], with per-bucket contributions
    double integrate(const std::vector<double>& times, const std::vector<double>& x, double sign,
        const CreditCurve& curve, double lgd, BucketRule rule,
        std::vector<double>* byBucket)
    {
        const std::size_t n = times.size();
        if (byBucket) byBucket->assign(n, 0.0);

        double total = 0.0;
        double sPrev = curve.getSurvivalProbability(times[0]);
        for (std::size_t k = 1; k < n; ++k) {
            const double s = curve.getSurvivalProbability(times[k]);
            const double contrib = lgd * sign * bucketValue(x, k, rule) * (sPrev - s);
            if (byBucket) (*byBucket)[k] = contrib;
            total += contrib;
            sPrev = s;
        }
        return total;
    }

} // namespace

double computeCVA(const ExposureProfile& e, const CreditCurve& counterparty, double lgd, BucketRule rule) {
    validate(e, lgd);
    return integrate(e.times, e.EE, +1.0, counterparty, lgd, rule, nullptr);
}

double computeDVA(const ExposureProfile& e, const CreditCurve& own, double lgd, BucketRule rule) {
    validate(e, lgd);
    return integrate(e.times, e.ENE, -1.0, own, lgd, rule, nullptr);   // ENE <= 0, so -ENE >= 0
}

XVAResult computeXVA(const ExposureProfile& e,
    const CreditCurve& counterparty, double lgdCounterparty,
    const CreditCurve& own, double lgdOwn,
    BucketRule rule)
{
    validate(e, lgdCounterparty);
    validate(e, lgdOwn);
    XVAResult r;
    r.cva = integrate(e.times, e.EE, +1.0, counterparty, lgdCounterparty, rule, &r.cvaByBucket);
    r.dva = integrate(e.times, e.ENE, -1.0, own, lgdOwn, rule, &r.dvaByBucket);
    r.bcva = r.cva - r.dva;
    return r;
}
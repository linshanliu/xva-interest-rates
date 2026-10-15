#ifndef FLATHAZARDCURVE_HPP
#define FLATHAZARDCURVE_HPP

#include <cmath>
#include <stdexcept>
#include "CreditCurve.hpp"

// Constant hazard rate lambda (flat CDS curve):  S(t) = exp(-lambda * t)
// Rule of thumb: CDS spread ~ lambda * LGD.
class FlatHazardCurve : public CreditCurve {
public:
    explicit FlatHazardCurve(double hazardRate) : lambda_(hazardRate) {
        if (lambda_ < 0.0) throw std::invalid_argument("FlatHazardCurve: hazard rate must be non-negative");
    }

    double getSurvivalProbability(double t) const override {
        return std::exp(-lambda_ * t);
    }

    double hazardRate() const { return lambda_; }

private:
    double lambda_;
};

#endif // !FLATHAZARDCURVE_HPP
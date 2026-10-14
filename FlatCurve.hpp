#ifndef FLATCURVE_HPP
#define FLATCURVE_HPP

#include <cmath>
#include "DiscountCurve.hpp"

class FlatCurve : public DiscountCurve {
public:
    explicit FlatCurve(double rate) : rate_(rate) {}

    double getDiscountFactor(double time) const override {
        return std::exp(-rate_ * time);
    }

    double rate() const { return rate_; }

	double instForward(double t) const override {
		return rate_;
	}

private:
    double rate_;
};

#endif
#ifndef DISCOUNTCURVE_HPP
#define DISCOUNTCURVE_HPP

#include <cmath>

class DiscountCurve {

public:
	DiscountCurve() = default;
	virtual ~DiscountCurve() = default;
	virtual double getDiscountFactor(double time) const = 0;


	// return the instantaneous forward rate at time t, i.e. f(t) = -d/dt log(P(0,t))
    virtual double instForward(double t) const {
        constexpr double h = 1e-4;
        if (t >= h) {
            return -(std::log(getDiscountFactor(t + h)) -
                std::log(getDiscountFactor(t - h))) / (2.0 * h);
        }

		// t close to 0: P(0, t-h) is not defined, use second-order one-sided difference instead, error is also O(h^2)
        double l0 = std::log(getDiscountFactor(t));
        double l1 = std::log(getDiscountFactor(t + h));
        double l2 = std::log(getDiscountFactor(t + 2.0 * h));
        return -(-3.0 * l0 + 4.0 * l1 - l2) / (2.0 * h);
    }
};


#endif // !DISCOUNTCURVE_HPP

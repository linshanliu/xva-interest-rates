#include "HullWhite1F.hpp"



HullWhite1F::HullWhite1F(std::shared_ptr<const DiscountCurve> curve, double a, double sigma)
	: curve_(curve), a_(a), sigma_(sigma) {
}

std::size_t HullWhite1F::stateDimension() const {
	return 1;
}



std::vector<double> HullWhite1F::initialState() const {
	return { curve_->instForward(0.0) };
}


//     theta(t) = d/dt f^M(0,t) + a * f^M(0,t) + sigma^2 / (2a) * ( 1 - exp(-2 a t) )
double HullWhite1F::theta(double t) const {
	constexpr double h = 1e-4;
	double dfdt = (curve_->instForward(t + h) - curve_->instForward(t)) / h;
	return dfdt + a_ * curve_->instForward(t)
		+ sigma_ * sigma_ / (2.0 * a_) * (1.0 - std::exp(-2.0 * a_ * t));
}


//     dr(t) = ( theta(t) - a * r(t) ) dt + sigma * dW(t),     r(0) = f^M(0,0)
std::vector<double> HullWhite1F::step(
	double t,
	double dt,
	std::vector<double>& curr,
	std::vector<double>& normalDraws
)const {

	std::vector<double> nextState(1);
	nextState[0] = curr[0] + (theta(t) - a_*curr[0]) * dt + sigma_ * sqrt(dt) * normalDraws[0];
	return nextState;
}


double HullWhite1F::zcb(double t, double T, double r) const{
	// Calculate the zero-coupon bond price P(t,T) under the Hull-White model
	double B = (1.0 - std::exp(-a_ * (T - t))) / a_;
	double A = curve_->getDiscountFactor(T) / curve_->getDiscountFactor(t) *
		std::exp(B * curve_->instForward(t) -   sigma_ * sigma_ * B * B * (1.0 - std::exp(-2*a_*t))/ (4*a_));
	return A * std::exp(-B * r);
}

double HullWhite1F::a() const {
	return a_;
}

double HullWhite1F::sigma() const {
	return sigma_;
}

const DiscountCurve& HullWhite1F::curve() const {
	return *curve_;
}

double HullWhite1F::shortRate(double /*t*/, const std::vector<double>& state) const {
	return state[0];   // state is the short rate itself
}
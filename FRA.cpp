#include "FRA.hpp"

#include <stdexcept>

FRA::FRA(double T1, double T2, double K, int fixingId)
	: FRA(T1, T2, K, T2 - T1, fixingId) {
}

FRA::FRA(double T1, double T2, double K, double tau, int fixingId)
	: T1_(T1), T2_(T2), K_(K), tau_(tau), fixingId_(fixingId)
{
	if (T1_ < 0.0 || T2_ <= T1_) throw std::invalid_argument("FRA: require 0 <= T1 < T2");
	if (tau_ <= 0.0)             throw std::invalid_argument("FRA: tau must be positive");
}


double FRA::markToMarket(const MarketView& m) const{
	if (m.time() < T1_) return m.discount(T1_) - (1.0 + K_ * tau_) * m.discount(T2_);
	else if (m.time() < T2_) {
		double L = m.fixing(fixingId_);
		return (L - K_) * tau_ * m.discount(T2_);
	}
	else return 0.0;
}


double FRA::maturity() const {
	return T2_;
}


std::vector<double> FRA::eventDates() const {
	return { T1_, T2_ };
}

std::vector<FixingSpec> FRA::fixings() const {
	return { {fixingId_, T1_, T1_, T2_, tau_} };
}

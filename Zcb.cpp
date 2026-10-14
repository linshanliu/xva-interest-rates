#include "Zcb.hpp"
#include <stdexcept>

ZCB::ZCB(double T) : T_(T) {
	if (T_ < 0.0) throw std::invalid_argument("ZCB: maturity must be non-negative");
}

double ZCB::markToMarket(const MarketView& m) const {
	return m.time() < T_ ? m.discount(T_) : 0.0;
}
double ZCB::maturity() const {
	return T_;
}


std::vector<double> ZCB::eventDates() const {
	return { T_ };
}
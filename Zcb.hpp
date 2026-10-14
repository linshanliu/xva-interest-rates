#ifndef ZCB_HPP
#define ZCB_HPP

#include "IInstrument.hpp"


class ZCB : public Instrument {
public:
    explicit ZCB(double T);

    double markToMarket(const MarketView& m) const override;
    double maturity() const override;
    std::vector<double> eventDates() const override;

	// override? no, already implemented the default version in the base class
    // std::vector<FixingSpec> fixings() const { return {}; }

private:
	double T_; // maturity
};

#endif // !ZCB_HPP

#ifndef IINSTRUMENT_HPP
#define IINSTRUMENT_HPP

#include "MarketView.hpp"
#include <vector>


// Mainly for FRA, IRS, cap/floor the instrument needs to know 
// the fixing dates and the corresponding fixing ids.
struct FixingSpec {
    int    id;
    double fixDate;
    double start;
    double end;
    double tau;
};


class Instrument {
public:
    virtual ~Instrument() = default;

	// price at current time t, given market view m
    virtual double markToMarket(const MarketView& m) const = 0;


    virtual double maturity() const = 0;
    virtual std::vector<double> eventDates() const = 0;
    virtual std::vector<FixingSpec> fixings() const { return {}; }
};

#endif // !IINSTRUMENT_HPP

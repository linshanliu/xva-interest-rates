#ifndef IRS_HPP
#define IRS_HPP

#include "IInstrument.hpp"


struct Schedule {
	std::vector<double> dates;   // size n + 1
	std::vector<double> taus;    // size n
};


class IRS : public Instrument {
public:
    IRS(double K, Schedule fixedLeg, Schedule floatLeg, int firstFixingId);

    double markToMarket(const MarketView& m) const override;
    double maturity() const override;
    std::vector<double> eventDates() const override;
    std::vector<FixingSpec> fixings() const override;

private:
    double   K_;               // fixed rate
    Schedule fixed_;           // S_0 < ... < S_m
    Schedule float_;           // T_0 < ... < T_n
    int      firstFixingId_;   // floating period j uses firstFixingId_ + (j - 1)
};



#endif // !IRS_HPP

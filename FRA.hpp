#ifndef FRA_HPP
#define FRA_HPP

#include <vector>
#include "IInstrument.hpp"


class FRA : public Instrument {
private:
	double T1_;   // start date of the FRA
	double T2_;   // end date of the FRA
	double K_;	  // fixed rate of the FRA
	double tau_;  // year fraction between start and end
	int fixingId_;// unique id for the fixing of the floating rate
public:
	FRA(double T1, double T2, double K, int fixingID);
	FRA(double T1, double T2, double K, double tau, int fixingID);

	double markToMarket(const MarketView& m) const override;
	double maturity() const override;
	std::vector<double> eventDates() const override;
	std::vector<FixingSpec> fixings() const override;
};


#endif // !FRA_HPP

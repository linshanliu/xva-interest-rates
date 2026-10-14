#include "IRS.hpp"
#include <algorithm>

#include <stdexcept>
#include <utility>

namespace {
	void validateSchedule(const Schedule& s, const char* name) {
		if (s.dates.size() < 2 || s.taus.size() != s.dates.size() - 1)
			throw std::invalid_argument(std::string("IRS: inconsistent ") + name + " schedule");
		for (std::size_t i = 1; i < s.dates.size(); ++i) {
			if (s.dates[i] <= s.dates[i - 1])
				throw std::invalid_argument(std::string("IRS: ") + name + " dates must be increasing");
			if (s.taus[i - 1] <= 0.0)
				throw std::invalid_argument(std::string("IRS: ") + name + " taus must be positive");
		}
	}
} // namespace

IRS::IRS(double K, Schedule fixedLeg, Schedule floatLeg, int firstFixingId)
	: K_(K), fixed_(std::move(fixedLeg)), float_(std::move(floatLeg)), firstFixingId_(firstFixingId)
{
	validateSchedule(fixed_, "fixed");
	validateSchedule(float_, "float");
	if (fixed_.dates.front() != float_.dates.front() || fixed_.dates.back() != float_.dates.back())
		throw std::invalid_argument("IRS: legs must share start and end dates");
	if (float_.dates.front() < 0.0)
		throw std::invalid_argument("IRS: start date must be non-negative");
}

//IRS::IRS(double K, Schedule fixedLeg, Schedule floatLeg, int firstFixingId) : K_(K), fixed_(fixedLeg),
//	float_(floatLeg), firstFixingId_(firstFixingId){}


// 1 Structure 中的 data 用.去获得
// 2 如果需要用自己的 method 可以不用提 object

double IRS::markToMarket(const MarketView& m) const {
	const std::vector<double>&  S = fixed_.dates;  // S_0 < ... < S_m
	const std::vector<double>&  T = float_.dates;  // T_0 < ... < T_n
	const double t = m.time();

	if (t >= maturity()) { return 0; }

	// ------ fixed leg -------
	double PV_fixed = 0.0;
	for (std::size_t i = 1; i < S.size(); i++) {
		if (t < S[i]) {
			PV_fixed += K_ * fixed_.taus[i - 1] * m.discount(S[i]);
		}
	}

	// ------ floating leg -------
	double PV_float = 0.0;
	if (t < T.front()) {
		PV_float = m.discount(T.front()) - m.discount(T.back());
	}
	else if (t < T.back()) {
		// find j* such that T_{j*-1} <= t < T_{j*}
		// upper_bound returns an iterator to the first element in the range [first, last)
		// that is greater than value, or last if no such element is found.
		auto it = std::upper_bound(T.begin(), T.end(), t);

		// convert iterator to index
		const std::size_t jStar = static_cast<std::size_t>(it - T.begin());


		const double tau = float_.taus[jStar - 1];
		const double L = m.fixing(firstFixingId_ + static_cast<int>(jStar) - 1);
		PV_float = (1.0 + tau * L) * m.discount(T[jStar]) - m.discount(T.back());
	}

	return PV_float - PV_fixed;

}
double IRS::maturity() const {
	return std::max(fixed_.dates.back(), float_.dates.back());
}
std::vector<double> IRS::eventDates() const {
	// create a vector to hold event dates
	std::vector<double> dates;

	// reserve enough space for all dates from both legs
	dates.reserve(fixed_.dates.size() + float_.dates.size());

	// insert dates from both legs into the vector
	dates.insert(dates.end(), fixed_.dates.begin(), fixed_.dates.end());
	dates.insert(dates.end(), float_.dates.begin(), float_.dates.end());

	// sort the vector remove duplicates
	std::sort(dates.begin(), dates.end());

	//remove duplicates using std::unique and erase
	dates.erase(std::unique(dates.begin(), dates.end()), dates.end());
	return dates;
}


std::vector<FixingSpec> IRS::fixings() const {

	const std::vector<double>& T = float_.dates;
	const std::size_t n = float_.taus.size();

	// create a vector of FixingSpec for the floating leg
	std::vector<FixingSpec> specs;


	specs.reserve(n);
	for (std::size_t j = 1; j <= n; ++j) {
		specs.push_back({
			firstFixingId_ + static_cast<int>(j) - 1,   // id
			T[j - 1],                                   // fixDate
			T[j - 1],                                   // start
			T[j],                                       // end
			float_.taus[j - 1]                          // tau
			});
	}
	return specs;
}
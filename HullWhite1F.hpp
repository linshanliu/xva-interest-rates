#ifndef HULLWHITE1F_HPP
#define HULLWHITE1F_HPP

#include "IStochasticModel.hpp"
#include "DiscountCurve.hpp"
#include "IShortRateModel.hpp"
#include <memory>


class HullWhite1F : public IShortRateModel {
public:
	HullWhite1F() =delete;
	HullWhite1F(std::shared_ptr<const DiscountCurve> curve, double a, double sigma);

	~HullWhite1F() = default;

	std::size_t stateDimension() const override;
	std::vector<double> initialState() const override;

	double theta(double t) const;

	std::vector<double> step(
		double t,
		double dt,
		std::vector<double>& curr,
		std::vector<double>& normalDraws
	)const override;
	
	double zcb(double t, double T, double r) const;

	// Getters
	double a() const;
	double sigma() const;
	const DiscountCurve& curve() const;
	double shortRate(double t, const std::vector<double>& state) const override;
private:
	// DiscountCurve is a abstract class so we use a pointer.
	// We use a shared_ptr to manage the lifetime of the DiscountCurve object
	// because one curve can be used by multiple objects.
	std::shared_ptr<const DiscountCurve> curve_;
	double a_; // mean reversion speed
	double sigma_; // volatility

	
};




#endif // !HULLWHITE1F_HPP

#ifndef CREDITCURVE_HPP
#define CREDITCURVE_HPP


class CreditCurve {
public:
	CreditCurve() = default;
	virtual ~CreditCurve() = default;
	virtual double getSurvivalProbability(double time) const = 0;

};



#endif // !CREDITCURVE_HPP

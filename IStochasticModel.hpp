#ifndef ISTOCHASTICMODEL_HPP
#define ISTOCHASTICMODEL_HPP
#include <cstddef>
#include <vector>


class IStochasticModel {

public:
	IStochasticModel() = default;
	virtual ~IStochasticModel() = default;
	virtual std::size_t stateDimension() const = 0;
	virtual std::vector<double> initialState() const = 0;

	virtual std::vector<double> step(
		double t,
		double dt,
		std::vector<double>& curr,
		std::vector<double>& normaldraws
	)const = 0;



};


#endif // !ISTOCHASTICMODEL_HPP

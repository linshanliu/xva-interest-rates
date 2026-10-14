#ifndef ISHORTRATEMODEL_HPP
#define ISHORTRATEMODEL_HPP

#include <vector>
#include "IStochasticModel.hpp"

// A stochastic model whose state determines an instantaneous short rate r(t).
// PathGenerator needs r(t) along each path to build the stochastic discount factor
//     D(0,t) = exp( - int_0^t r(u) du ).
// Commodity models (Gibson-Schwartz) stay plain IStochasticModel; rates models derive from this.
class IShortRateModel : public IStochasticModel {
public:
    ~IShortRateModel() override = default;

    // Short rate at time t given the state vector
    virtual double shortRate(double t, const std::vector<double>& state) const = 0;
};

#endif // !ISHORTRATEMODEL_HPP
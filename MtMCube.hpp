#ifndef MTMCUBE_HPP
#define MTMCUBE_HPP

#include <cstddef>
#include <vector>

// Mark-to-market values of every trade on every path at every grid time.
//
// Layout (flattened, contiguous):
//   values[trade][path][k]  ->  values[(i * nPaths + p) * nSteps + k]
//
// Values are undiscounted MtM at t_k (unit notional, holder's perspective),
// excluding cash flows paid at or before t_k (ex-coupon).
struct MtMCube {
    std::vector<double> times;
    std::size_t nTrades = 0;
    std::size_t nPaths = 0;
    std::size_t nSteps = 0;
    std::vector<double> values;

    double value(std::size_t i, std::size_t p, std::size_t k) const {
        return values[(i * nPaths + p) * nSteps + k];
    }
    double& value(std::size_t i, std::size_t p, std::size_t k) {
        return values[(i * nPaths + p) * nSteps + k];
    }
};

#endif // !MTMCUBE_HPP
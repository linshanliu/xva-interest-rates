#ifndef SCENARIOSET_HPP
#define SCENARIOSET_HPP

#include <cstddef>
#include <vector>

// Output of the path generation step: simulated state paths and discount factors.
//
// Layout (flattened, contiguous):
//   states  [path][k][d]  ->  states[(p * nSteps + k) * dim + d]
//   discount[path][k]     ->  discount[p * nSteps + k]     D(0, t_k) = exp(-int_0^{t_k} r du)
//
// times[0] is the valuation date (0.0); nSteps == times.size().
struct ScenarioSet {
    std::vector<double> times;
    std::size_t nPaths = 0;
    std::size_t nSteps = 0;
    std::size_t dim = 0;
    std::vector<double> states;
    std::vector<double> discount;

    double state(std::size_t p, std::size_t k, std::size_t d = 0) const {
        return states[(p * nSteps + k) * dim + d];
    }

    // Copy of the full state vector at (p, k), e.g. to pass to a model method
    std::vector<double> stateVector(std::size_t p, std::size_t k) const {
        const double* first = &states[(p * nSteps + k) * dim];
        return std::vector<double>(first, first + dim);
    }

    double df(std::size_t p, std::size_t k) const {
        return discount[p * nSteps + k];
    }
};

#endif // !SCENARIOSET_HPP
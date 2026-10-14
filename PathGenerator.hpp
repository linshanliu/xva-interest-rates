#ifndef PATHGENERATOR_HPP
#define PATHGENERATOR_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "IInstrument.hpp"
#include "IShortRateModel.hpp"
#include "ScenarioSet.hpp"

// Builds the simulation grid: {0} U baseGrid U all instrument event dates,
// sorted, de-duplicated (tolerance tol), and truncated at the horizon.
// Event dates guarantee that payment and fixing dates lie exactly on the grid.
std::vector<double> buildTimeGrid(const std::vector<double>& baseGrid,
    const std::vector<std::shared_ptr<const Instrument>>& trades,
    double horizon,
    double tol = 1e-10);

// Generates state paths for a short-rate model on a given grid.
// Knows nothing about instruments: it only evolves the model and accumulates
// the discount factor D(0,t) = exp(-int r du) with the trapezoid rule.
class PathGenerator {
public:
    PathGenerator(std::shared_ptr<const IShortRateModel> model, std::uint64_t seed);

    // grid must start at 0 and be strictly increasing
    ScenarioSet generate(const std::vector<double>& grid, std::size_t nPaths) const;

private:
    std::shared_ptr<const IShortRateModel> model_;
    std::uint64_t seed_;
};

#endif // !PATHGENERATOR_HPP
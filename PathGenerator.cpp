#include "PathGenerator.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <utility>

std::vector<double> buildTimeGrid(const std::vector<double>& baseGrid,
    const std::vector<std::shared_ptr<const Instrument>>& trades,
    double horizon,
    double tol)
{
    // (time, isEventDate). Event dates win over base-grid points that coincide within tol,
    // so payment/fixing dates appear on the grid exactly as the instruments define them.
    std::vector<std::pair<double, bool>> pts;
    pts.emplace_back(0.0, true);
    for (double t : baseGrid) pts.emplace_back(t, false);
    for (const auto& trade : trades)
        for (double t : trade->eventDates()) pts.emplace_back(t, true);

    // keep only 0 <= t <= horizon
    pts.erase(std::remove_if(pts.begin(), pts.end(),
        [&](const std::pair<double, bool>& x) { return x.first < 0.0 || x.first > horizon + tol; }),
        pts.end());

    std::sort(pts.begin(), pts.end());

    // de-duplicate with tolerance, preferring event dates
    std::vector<std::pair<double, bool>> merged;
    merged.reserve(pts.size());
    for (const auto& x : pts) {
        if (!merged.empty() && x.first - merged.back().first <= tol) {
            if (x.second && !merged.back().second) merged.back() = x;
            continue;
        }
        merged.push_back(x);
    }

    std::vector<double> out;
    out.reserve(merged.size());
    for (const auto& x : merged) out.push_back(x.first);
    return out;
}

PathGenerator::PathGenerator(std::shared_ptr<const IShortRateModel> model, std::uint64_t seed)
    : model_(std::move(model)), seed_(seed)
{
    if (!model_) throw std::invalid_argument("PathGenerator: model is null");
}

ScenarioSet PathGenerator::generate(const std::vector<double>& grid, std::size_t nPaths) const {
    if (grid.empty() || grid.front() != 0.0)
        throw std::invalid_argument("PathGenerator: grid must start at 0");
    for (std::size_t k = 1; k < grid.size(); ++k)
        if (grid[k] <= grid[k - 1])
            throw std::invalid_argument("PathGenerator: grid must be strictly increasing");

    const std::size_t dim = model_->stateDimension();
    const std::size_t nSteps = grid.size();

    ScenarioSet scen;
    scen.times = grid;
    scen.nPaths = nPaths;
    scen.nSteps = nSteps;
    scen.dim = dim;
    scen.states.resize(nPaths * nSteps * dim);
    scen.discount.resize(nPaths * nSteps);

    std::mt19937_64 rng(seed_);
    std::normal_distribution<double> N01(0.0, 1.0);
    std::vector<double> z(dim);   // one normal per state dimension per step

    for (std::size_t p = 0; p < nPaths; ++p) {
        std::vector<double> s = model_->initialState();
        double integral = 0.0;

        for (std::size_t k = 0; k < nSteps; ++k) {
            // store state and discount factor at t_k
            std::copy(s.begin(), s.end(), scen.states.begin() + (p * nSteps + k) * dim);
            scen.discount[p * nSteps + k] = std::exp(-integral);

            // evolve to t_{k+1}
            if (k + 1 < nSteps) {
                const double t = grid[k];
                const double dt = grid[k + 1] - grid[k];
                for (double& zi : z) zi = N01(rng);

                std::vector<double> next = model_->step(t, dt, s, z);
                integral += 0.5 * (model_->shortRate(t, s) + model_->shortRate(t + dt, next)) * dt;
                s = std::move(next);
            }
        }
    }
    return scen;
}
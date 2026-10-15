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
	// in pts, first = time, second = isEventDate
    std::vector<std::pair<double, bool>> pts;


    // add 0.0, true in the end of the pts
    pts.emplace_back(0.0, true);

	// for baseGrid add t and false, for trades add all eventDates and true
    for (double t : baseGrid) pts.emplace_back(t, false);
    for (const auto& trade : trades)
        for (double t : trade->eventDates()) pts.emplace_back(t, true);

    // keep only 0 <= t <= horizon
	// std::remove_if(first, last, pred) removes all elements for which pred returns true to the end, and returns an iterator to the begin of remove part.
	// vec.erase(from, to) removes all elements in the range [from, to) from the vector.
    pts.erase(std::remove_if(pts.begin(), pts.end(),
        [&](const std::pair<double, bool>& x) { return x.first < 0.0 || x.first > horizon + tol; }),
        pts.end());

    std::sort(pts.begin(), pts.end());

    // de-duplicate with tolerance, preferring event dates
    std::vector<std::pair<double, bool>> merged;
    merged.reserve(pts.size());
    for (const auto& x : pts) {

		// if the merged vector is not empty and the difference between the current point and the last point in merged is less than or equal to tol
        if (!merged.empty() && x.first - merged.back().first <= tol) {
			// if the current point is an event date and the last point in merged is not, replace the last point in merged with the current point
            if (x.second && !merged.back().second) merged.back() = x;
            continue;
        }
        merged.push_back(x);
    }


	// extract the first element of each pair in merged to create the output vector
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


// given a grid of times and a number of paths, generate the state paths and discount factors
ScenarioSet PathGenerator::generate(const std::vector<double>& grid, std::size_t nPaths) const {
	// if the grid is empty or the first element of the grid is not 0.0, throw an invalid_argument exception (time grid must start by 0)
    if (grid.empty() || grid.front() != 0.0)
        throw std::invalid_argument("PathGenerator: grid must start at 0");

	// check that the grid of time is strictly increasing
    for (std::size_t k = 1; k < grid.size(); ++k)
        if (grid[k] <= grid[k - 1])
            throw std::invalid_argument("PathGenerator: grid must be strictly increasing");

    const std::size_t dim = model_->stateDimension();
    const std::size_t nSteps = grid.size();


	// define a ScenarioSet object to hold the data
    ScenarioSet scen;
    scen.times = grid;
    scen.nPaths = nPaths;
    scen.nSteps = nSteps;
    scen.dim = dim;
    scen.states.resize(nPaths * nSteps * dim);
    scen.discount.resize(nPaths * nSteps);


	// initialize a random number generator with the given seed 
    // and a normal distribution with mean 0 and standard deviation 1
    std::mt19937_64 rng(seed_);
    std::normal_distribution<double> N01(0.0, 1.0);
    std::vector<double> z(dim);


	// loop over each path
    for (std::size_t p = 0; p < nPaths; ++p) {

		// initialize the state vector to the initial state of the model and the integral of the short rate to 0
        std::vector<double> s = model_->initialState();
        double integral = 0.0;


		// loop over each time step
        for (std::size_t k = 0; k < nSteps; ++k) {
			// copy the current state vector to the ScenarioSet states
            std::copy(s.begin(), s.end(), scen.states.begin() + (p * nSteps + k) * dim);
			// compute the discount factor for the current time step using the trapezoidal rule and store it in the ScenarioSet discount
            scen.discount[p * nSteps + k] = std::exp(-integral);

			// if the current time step is not the last one, update t and dt,
            // generate a new set of standard normal random variables, and step the model forward to the next state
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
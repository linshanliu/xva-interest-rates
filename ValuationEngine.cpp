#include "ValuationEngine.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "FixingStore.hpp"
#include "HW1FView.hpp"

namespace {

    // Index k with |times[k] - t| <= tol, or -1 if none (times strictly increasing)
    long findGridIndex(const std::vector<double>& times, double t, double tol) {
        auto it = std::lower_bound(times.begin(), times.end(), t - tol);
        if (it != times.end() && std::abs(*it - t) <= tol) return static_cast<long>(it - times.begin());
        return -1;
    }

    bool sameFixing(const FixingSpec& x, const FixingSpec& y, double tol) {
        return std::abs(x.fixDate - y.fixDate) <= tol && std::abs(x.start - y.start) <= tol
            && std::abs(x.end - y.end) <= tol && std::abs(x.tau - y.tau) <= tol;
    }

} // namespace

ValuationEngine::ValuationEngine(std::shared_ptr<const HullWhite1F> model, double tol)
    : model_(std::move(model)), tol_(tol)
{
    if (!model_) throw std::invalid_argument("ValuationEngine: model is null");
}

MtMCube ValuationEngine::run(const ScenarioSet& scen,
    const std::vector<std::shared_ptr<const Instrument>>& trades) const
{
    const std::size_t nSteps = scen.nSteps;
    const std::size_t nPaths = scen.nPaths;
    const std::size_t nTrades = trades.size();
    if (nSteps == 0 || scen.times.size() != nSteps)
        throw std::invalid_argument("ValuationEngine: empty or inconsistent scenario set");
    const double horizon = scen.times.back();

    // ---- 1. Collect fixings, check ids, map each fixing to its grid step ----
    std::vector<FixingSpec> byId;                         // byId[id] = spec (fixDate < 0 means unused)
    std::vector<std::vector<FixingSpec>> fixingsAtStep(nSteps);

    for (const auto& trade : trades) {
        for (const FixingSpec& f : trade->fixings()) {
            if (f.id < 0) throw std::invalid_argument("ValuationEngine: negative fixing id");
            const std::size_t id = static_cast<std::size_t>(f.id);
            if (id >= byId.size()) byId.resize(id + 1, FixingSpec{ -1, -1.0, 0.0, 0.0, 0.0 });

            if (byId[id].fixDate >= 0.0) {                // id already registered
                if (!sameFixing(byId[id], f, tol_))
                    throw std::invalid_argument("ValuationEngine: fixing id " + std::to_string(f.id)
                        + " used for two different fixings");
                continue;                                 // same fixing shared: record once
            }
            byId[id] = f;

            if (f.fixDate > horizon + tol_) continue;     // beyond the simulation horizon
            const long k = findGridIndex(scen.times, f.fixDate, tol_);
            if (k < 0)
                throw std::invalid_argument("ValuationEngine: fixing date " + std::to_string(f.fixDate)
                    + " is not on the simulation grid");
            fixingsAtStep[static_cast<std::size_t>(k)].push_back(f);
        }
    }

    // ---- 2. Allocate output ----
    MtMCube cube;
    cube.times = scen.times;
    cube.nTrades = nTrades;
    cube.nPaths = nPaths;
    cube.nSteps = nSteps;
    cube.values.assign(nTrades * nPaths * nSteps, 0.0);

    // ---- 3. Path loop ----
    FixingStore fixings(byId.size());
    HW1FView view(*model_, fixings);

    for (std::size_t p = 0; p < nPaths; ++p) {
        fixings.clear();
        for (std::size_t k = 0; k < nSteps; ++k) {
            const double r = scen.state(p, k, 0);

            // a) fixings due at t_k: L = (1/P(T_s, T_e) - 1) / tau, evaluated at the exact fixing date
            for (const FixingSpec& f : fixingsAtStep[k]) {
                view.reset(f.fixDate, r);
                const double P = view.discount(f.end);
                fixings.set(f.id, (1.0 / P - 1.0) / f.tau);
            }

            // b) market view at t_k
            view.reset(scen.times[k], r);

            // c) value all trades
            for (std::size_t i = 0; i < nTrades; ++i)
                cube.value(i, p, k) = trades[i]->markToMarket(view);
        }
    }
    return cube;
}
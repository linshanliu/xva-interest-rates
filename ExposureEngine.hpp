#ifndef EXPOSUREENGINE_HPP
#define EXPOSUREENGINE_HPP

#include <memory>
#include <vector>

#include "ExposureCube.hpp"
#include "HullWhite1F.hpp"
#include "IInstrument.hpp"
#include "ScenarioSet.hpp"

// Values a set of trades on pre-generated scenarios.
//
// For each path and each grid point t_k, in this order:
//   a) record the floating-rate fixings due at t_k into the path's FixingStore
//   b) point the market view at (t_k, state)
//   c) value every trade
//
// Requirements (checked, std::invalid_argument otherwise):
//   - every fixing date within the grid horizon lies on the grid (within tol)
//   - fixing ids are non-negative; an id shared by two specs must describe the same fixing
class ExposureEngine {
public:
    explicit ExposureEngine(std::shared_ptr<const HullWhite1F> model, double tol = 1e-10);

    ExposureCube run(const ScenarioSet& scen,
        const std::vector<std::shared_ptr<const Instrument>>& trades) const;

private:
    std::shared_ptr<const HullWhite1F> model_;
    double tol_;
};

#endif // !EXPOSUREENGINE_HPP
#include "Exposure.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

    void validate(const MtMCube& cube, const NettingSet& ns) {
        if (ns.tradeIdx.size() != ns.quantity.size())
            throw std::invalid_argument("Exposure: netting set tradeIdx and quantity sizes differ");
        for (std::size_t idx : ns.tradeIdx)
            if (idx >= cube.nTrades)
                throw std::invalid_argument("Exposure: trade index " + std::to_string(idx) + " not in cube");
    }

} // namespace

std::vector<double> nettedValues(const MtMCube& cube, const NettingSet& ns) {
    validate(cube, ns);
    const std::size_t nPaths = cube.nPaths, nSteps = cube.nSteps;
    std::vector<double> v(nPaths * nSteps, 0.0);

    for (std::size_t j = 0; j < ns.tradeIdx.size(); ++j) {
        const std::size_t i = ns.tradeIdx[j];
        const double q = ns.quantity[j];
        for (std::size_t p = 0; p < nPaths; ++p)
            for (std::size_t k = 0; k < nSteps; ++k)
                v[p * nSteps + k] += q * cube.value(i, p, k);
    }
    return v;
}

ExposureProfile computeExposure(const MtMCube& cube,
    const ScenarioSet& scen,
    const NettingSet& ns,
    double pfeQuantile)
{
    if (!(pfeQuantile > 0.0 && pfeQuantile < 1.0))
        throw std::invalid_argument("Exposure: PFE quantile must be in (0, 1)");
    if (cube.nPaths != scen.nPaths || cube.nSteps != scen.nSteps)
        throw std::invalid_argument("Exposure: cube and scenario set have different shapes");
    if (cube.nPaths == 0)
        throw std::invalid_argument("Exposure: no paths");

    const std::size_t nPaths = cube.nPaths, nSteps = cube.nSteps;
    const std::vector<double> v = nettedValues(cube, ns);

    ExposureProfile prof;
    prof.times = cube.times;
    prof.pfeQuantile = pfeQuantile;
    prof.EE.assign(nSteps, 0.0);
    prof.ENE.assign(nSteps, 0.0);
    prof.EEu.assign(nSteps, 0.0);
    prof.PFE.assign(nSteps, 0.0);

    // index of the q-quantile in a sorted sample of size nPaths
    const std::size_t qIdx = std::min(nPaths - 1,
        static_cast<std::size_t>(std::ceil(pfeQuantile * static_cast<double>(nPaths))) - 1);

    std::vector<double> pos(nPaths);
    for (std::size_t k = 0; k < nSteps; ++k) {
        double ee = 0.0, ene = 0.0, eeu = 0.0;
        for (std::size_t p = 0; p < nPaths; ++p) {
            const double val = v[p * nSteps + k];
            const double d = scen.df(p, k);
            const double plus = std::max(val, 0.0);
            const double minus = std::min(val, 0.0);
            ee += d * plus;
            ene += d * minus;
            eeu += plus;
            pos[p] = plus;
        }
        prof.EE[k] = ee / static_cast<double>(nPaths);
        prof.ENE[k] = ene / static_cast<double>(nPaths);
        prof.EEu[k] = eeu / static_cast<double>(nPaths);

        std::nth_element(pos.begin(), pos.begin() + static_cast<std::ptrdiff_t>(qIdx), pos.end());
        prof.PFE[k] = pos[qIdx];
    }
    return prof;
}
// Tests for ExposureEngine / ExposureCube
//
// Assumes the classes from earlier steps: FlatCurve, HullWhite1F (derives IShortRateModel),
// FixingStore, HW1FView, ZCB, FRA(T1, T2, K, fixingId), IRS(K, fixed, float, firstFixingId),
// PathGenerator / buildTimeGrid / ScenarioSet.
//
// Note: exclude the other test files with their own main() from the build.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "FlatCurve.hpp"
#include "HullWhite1F.hpp"
#include "FixingStore.hpp"
#include "HW1FView.hpp"
#include "Zcb.hpp"
#include "FRA.hpp"
#include "IRS.hpp"
#include "PathGenerator.hpp"
#include "ExposureEngine.hpp"

namespace {

    int g_failures = 0;
    int g_tests = 0;

    void check(const std::string& name, double value, double expected, double tol) {
        ++g_tests;
        double err = std::abs(value - expected);
        bool ok = err <= tol;
        if (!ok) ++g_failures;
        std::cout << (ok ? "[PASS] " : "[FAIL] ") << name
            << "\n       value = " << std::setprecision(12) << value
            << ", expected = " << expected
            << ", |err| = " << std::setprecision(4) << err
            << ", tol = " << tol << "\n";
    }

    void checkTrue(const std::string& name, bool cond) {
        ++g_tests;
        if (!cond) ++g_failures;
        std::cout << (cond ? "[PASS] " : "[FAIL] ") << name << "\n";
    }

    void section(const std::string& title) {
        std::cout << "\n=== " << title << " ===\n";
    }

    Schedule makeTestSchedule(double start, double end, int paymentsPerYear) {
        const int n = static_cast<int>(std::round((end - start) * paymentsPerYear));
        Schedule s;
        for (int i = 0; i <= n; ++i) s.dates.push_back(start + i * (end - start) / n);
        for (int i = 1; i <= n; ++i) s.taus.push_back(s.dates[i] - s.dates[i - 1]);
        return s;
    }

    double simpleRate(double P, double tau) { return (1.0 / P - 1.0) / tau; }

    using Trades = std::vector<std::shared_ptr<const Instrument>>;

} // namespace

int main() {
    const double r0 = 0.03;
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    auto model = std::make_shared<HullWhite1F>(curve, a, sigma);
    ExposureEngine engine(model);

    // Portfolio used in most sections
    const double start = 1.0, end = 5.0;
    const Schedule fixedLeg = makeTestSchedule(start, end, 1);
    const Schedule floatLeg = makeTestSchedule(start, end, 4);
    double annuity0 = 0.0;
    for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i)
        annuity0 += fixedLeg.taus[i - 1] * curve->getDiscountFactor(fixedLeg.dates[i]);
    const double Kpar = (curve->getDiscountFactor(start) - curve->getDiscountFactor(end)) / annuity0;
    const double Kswap = Kpar + 0.005;
    const int swapFirstId = 5;

    Trades trades;
    trades.push_back(std::make_shared<ZCB>(3.0));                                          // trade 0
    trades.push_back(std::make_shared<FRA>(2.0, 2.5, 0.03, 0));                            // trade 1, id 0
    trades.push_back(std::make_shared<IRS>(Kswap, fixedLeg, floatLeg, swapFirstId));       // trade 2, ids 5..20

    // Weekly base grid 0 .. 5 merged with all event dates (event dates kept exactly)
    std::vector<double> weekly;
    for (int k = 1; k <= 52 * 5; ++k) weekly.push_back(k / 52.0);
    const std::vector<double> grid = buildTimeGrid(weekly, trades, 5.0);
    const int nGrid = static_cast<int>(grid.size()) - 1;

    // Grid index of time t (must be on the grid)
    auto idx = [&](double t) {
        for (std::size_t k = 0; k < grid.size(); ++k) if (std::abs(grid[k] - t) < 1e-10) return k;
        throw std::runtime_error("time not on grid: " + std::to_string(t));
        };

    // ------------------------------------------------------------------
    section("1. Input validation");
    // ------------------------------------------------------------------
    {
        ScenarioSet small = PathGenerator(model, 1).generate({ 0.0, 0.5, 1.0, 2.0, 3.0 }, 2);

        bool threwOffGrid = false;
        try {
            Trades t = { std::make_shared<FRA>(0.7, 1.2, 0.03, 0) };       // fixing at 0.7: not on grid
            engine.run(small, t);
        }
        catch (const std::invalid_argument&) { threwOffGrid = true; }
        checkTrue("fixing date off the grid throws", threwOffGrid);

        bool threwConflict = false;
        try {
            Trades t = { std::make_shared<FRA>(1.0, 1.5, 0.03, 0),
                         std::make_shared<FRA>(2.0, 2.5, 0.03, 0) };      // same id, different fixing
            engine.run(small, t);
        }
        catch (const std::invalid_argument&) { threwConflict = true; }
        checkTrue("same id for two different fixings throws", threwConflict);

        bool shared = true;
        try {
            Trades t = { std::make_shared<FRA>(1.0, 1.5, 0.03, 0),
                         std::make_shared<FRA>(1.0, 1.5, 0.04, 0) };      // same fixing, same id: allowed
            engine.run(small, t);
        }
        catch (const std::invalid_argument&) { shared = false; }
        checkTrue("same id for the same fixing is allowed", shared);

        bool beyond = true;
        try {
            Trades t = { std::make_shared<FRA>(4.0, 4.5, 0.03, 0) };       // fixing beyond horizon 3.0
            engine.run(small, t);
        }
        catch (const std::invalid_argument&) { beyond = false; }
        checkTrue("fixing beyond the horizon is ignored (no throw)", beyond);
    }

    // ------------------------------------------------------------------
    const std::size_t nPaths = 20000;
    ScenarioSet scen = PathGenerator(model, 42).generate(grid, nPaths);
    ExposureCube cube = engine.run(scen, trades);

    // ------------------------------------------------------------------
    section("2. Cube shape");
    // ------------------------------------------------------------------
    checkTrue("nTrades / nPaths / nSteps",
        cube.nTrades == 3 && cube.nPaths == nPaths && cube.nSteps == grid.size());
    checkTrue("values size", cube.values.size() == 3 * nPaths * grid.size());
    checkTrue("times copied from scenarios", cube.times == scen.times);

    // ------------------------------------------------------------------
    section("3. t = 0: every path equals the analytic value");
    // ------------------------------------------------------------------
    {
        FixingStore fx(32);
        HW1FView view(*model, fx);
        view.reset(0.0, model->initialState()[0]);
        for (std::size_t i = 0; i < trades.size(); ++i) {
            const double V0 = trades[i]->markToMarket(view);
            double maxErr = 0.0;
            for (std::size_t p = 0; p < nPaths; ++p) maxErr = std::max(maxErr, std::abs(cube.value(i, p, 0) - V0));
            check("trade " + std::to_string(i) + ": max_p |V(p,0) - V(0)|", maxErr, 0.0, 1e-15);
        }
    }

    // ------------------------------------------------------------------
    section("4. Same numbers as a hand-written loop (first 200 paths, all steps)");
    // ------------------------------------------------------------------
    {
        const std::size_t nCheck = 200;
        std::mt19937_64 rng(42);
        std::normal_distribution<double> N01(0.0, 1.0);
        FixingStore fx(32);
        HW1FView view(*model, fx);
        std::vector<double> z(1);
        double maxErr = 0.0;

        for (std::size_t p = 0; p < nCheck; ++p) {
            fx.clear();
            std::vector<double> s = model->initialState();
            for (int k = 0; k <= nGrid; ++k) {
                const double t = grid[k];
                for (const auto& tr : trades)
                    for (const FixingSpec& f : tr->fixings())
                        if (std::abs(f.fixDate - t) < 1e-10) {
                            view.reset(f.fixDate, s[0]);
                            fx.set(f.id, simpleRate(view.discount(f.end), f.tau));
                        }
                view.reset(t, s[0]);
                for (std::size_t i = 0; i < trades.size(); ++i)
                    maxErr = std::max(maxErr, std::abs(trades[i]->markToMarket(view) - cube.value(i, p, k)));
                if (k < nGrid) {
                    z[0] = N01(rng);
                    s = model->step(t, grid[k + 1] - grid[k], s, z);
                }
            }
        }
        check("max |V(manual) - V(engine)| over trades, paths, steps", maxErr, 0.0, 1e-14);
    }

    // ------------------------------------------------------------------
    section("5. Martingale checks from the cube");
    //   ZCB (single payment at 3): E[D(0,t) V(t)] = V(0) for every t < 3
    //   IRS: same identity before its first payment at 1.25
    // ------------------------------------------------------------------
    {
        auto meanDV = [&](std::size_t i, std::size_t k, double& se) {
            double s1 = 0.0, s2 = 0.0;
            for (std::size_t p = 0; p < nPaths; ++p) {
                double x = scen.df(p, k) * cube.value(i, p, k);
                s1 += x; s2 += x * x;
            }
            double m = s1 / nPaths;
            se = std::sqrt((s2 / nPaths - m * m) / nPaths);
            return m;
            };

        for (double t : { 1.0, 2.0, 2.75 }) {
            double se, m = meanDV(0, idx(t), se);
            check("ZCB: E[D V] at t = " + std::to_string(t), m, cube.value(0, 0, 0), 4.0 * se + 1e-4);
        }
        for (double t : { 0.5, 1.0, 1.25 - 1.0 / 52.0 }) {
            const std::size_t k = (t == 0.5 || t == 1.0) ? idx(t) : idx(1.25) - 1;   // last point before 1.25
            double se, m = meanDV(2, k, se);
            check("IRS: E[D V] at t = " + std::to_string(grid[k]), m, cube.value(2, 0, 0), 4.0 * se + 1e-4);
        }
    }

    // ------------------------------------------------------------------
    section("6. Event-date behaviour visible in the cube");
    // ------------------------------------------------------------------
    {
        auto meanV = [&](std::size_t i, std::size_t k) {
            double s = 0.0;
            for (std::size_t p = 0; p < nPaths; ++p) s += cube.value(i, p, k);
            return s / nPaths;
            };

        // ZCB(3) and FRA(2, 2.5) are worth 0 at and after maturity on every path
        bool zcbZero = true, fraZero = true, zcbAlive = true;
        for (std::size_t p = 0; p < nPaths; ++p) {
            zcbZero = zcbZero && cube.value(0, p, idx(3.0)) == 0.0 && cube.value(0, p, nGrid) == 0.0;
            zcbAlive = zcbAlive && cube.value(0, p, idx(3.0) - 1) > 0.0;
            fraZero = fraZero && cube.value(1, p, idx(2.5)) == 0.0;
        }
        checkTrue("ZCB worth > 0 one step before t = 3, all paths", zcbAlive);
        checkTrue("ZCB worth 0 from t = 3 on, all paths", zcbZero);
        checkTrue("FRA worth 0 from t = 2.5 on, all paths", fraZero);

        // Sawtooth: at t = 2 the payer swap pays fixed alpha*K and receives tau*L.
        // Mean MtM jumps up by about E[alpha K - tau L] (a one-week drift is small).
        const std::size_t k2 = idx(2.0);
        const double jump = meanV(2, k2) - meanV(2, k2 - 1);
        const double approx = 1.0 * Kswap - 0.25 * Kpar;
        std::cout << "  mean IRS MtM at t = 2 - dt: " << meanV(2, k2 - 1)
            << ", at t = 2: " << meanV(2, k2) << "\n";
        check("IRS mean MtM jumps up at fixed payment date t = 2", jump, approx, 3e-3);
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
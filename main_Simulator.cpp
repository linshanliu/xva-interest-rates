// Tests for buildTimeGrid, PathGenerator and ScenarioSet
//
// Assumes:
//   - HullWhite1F derives from IShortRateModel and implements
//       shortRate(double t, const std::vector<double>& s) const { return s[0]; }
//   - IStochasticModel has stateDimension(), initialState(),
//       step(double t, double dt, const std::vector<double>& curr,
//            const std::vector<double>& normalDraws) const
//   - ZCB(double T), IRS(double K, Schedule fixed, Schedule floating, int firstFixingId)
//
// Note: exclude the other test files with their own main() from the build.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "FlatCurve.hpp"
#include "HullWhite1F.hpp"
#include "PathGenerator.hpp"
#include "ScenarioSet.hpp"
#include "Zcb.hpp"
#include "IRS.hpp"

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

    std::vector<double> uniformGrid(double horizon, int stepsPerYear) {
        const int n = static_cast<int>(std::round(horizon * stepsPerYear));
        std::vector<double> g;
        for (int i = 1; i <= n; ++i) g.push_back(static_cast<double>(i) / stepsPerYear);
        return g;
    }

    bool onGrid(const std::vector<double>& grid, double t, double tol = 1e-10) {
        for (double g : grid) if (std::abs(g - t) <= tol) return true;
        return false;
    }

} // namespace

int main() {
    const double r0 = 0.03;
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    auto model = std::make_shared<HullWhite1F>(curve, a, sigma);

    // ------------------------------------------------------------------
    section("1. buildTimeGrid");
    // ------------------------------------------------------------------
    {
        std::vector<std::shared_ptr<const Instrument>> trades;
        trades.push_back(std::make_shared<ZCB>(2.3));                       // off the monthly grid
        trades.push_back(std::make_shared<IRS>(0.03, makeTestSchedule(1.0, 5.0, 1),
            makeTestSchedule(1.0, 5.0, 4), 0));
        trades.push_back(std::make_shared<ZCB>(8.0));                       // beyond the horizon

        const double horizon = 5.0;
        auto grid = buildTimeGrid(uniformGrid(horizon, 12), trades, horizon);

        checkTrue("starts at 0", !grid.empty() && grid.front() == 0.0);
        bool inc = true;
        for (std::size_t k = 1; k < grid.size(); ++k) inc = inc && (grid[k] > grid[k - 1] + 1e-10);
        checkTrue("strictly increasing, no near-duplicates", inc);
        checkTrue("last point <= horizon", grid.back() <= horizon + 1e-10);

        bool allEvents = true;
        for (const auto& tr : trades)
            for (double t : tr->eventDates())
                if (t <= horizon) allEvents = allEvents && onGrid(grid, t);
        checkTrue("every event date <= horizon is on the grid", allEvents);
        checkTrue("ZCB(2.3) maturity added", onGrid(grid, 2.3));
        checkTrue("ZCB(8.0) maturity dropped (beyond horizon)", !onGrid(grid, 8.0));
        check("grid size: 60 monthly + {0} + {2.3}", static_cast<double>(grid.size()), 62.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("2. ScenarioSet shape and initial values");
    // ------------------------------------------------------------------
    {
        PathGenerator gen(model, 7);
        std::vector<double> grid = { 0.0, 0.5, 1.0, 2.0 };
        ScenarioSet scen = gen.generate(grid, 100);

        checkTrue("nPaths / nSteps / dim", scen.nPaths == 100 && scen.nSteps == 4 && scen.dim == 1);
        checkTrue("states size == nPaths * nSteps * dim", scen.states.size() == 400);
        checkTrue("discount size == nPaths * nSteps", scen.discount.size() == 400);

        bool init = true;
        for (std::size_t p = 0; p < scen.nPaths; ++p)
            init = init && std::abs(scen.state(p, 0) - r0) < 1e-12 && scen.df(p, 0) == 1.0;
        checkTrue("every path starts at r0 with D(0,0) = 1", init);
    }

    // ------------------------------------------------------------------
    section("3. Invalid grids are rejected");
    // ------------------------------------------------------------------
    {
        PathGenerator gen(model, 7);
        bool threw1 = false, threw2 = false;
        try { gen.generate({ 0.5, 1.0 }, 10); }
        catch (const std::invalid_argument&) { threw1 = true; }
        try { gen.generate({ 0.0, 1.0, 1.0 }, 10); }
        catch (const std::invalid_argument&) { threw2 = true; }
        checkTrue("grid not starting at 0 throws", threw1);
        checkTrue("non-increasing grid throws", threw2);
    }

    // ------------------------------------------------------------------
    section("4. Deterministic limit (sigma -> 0)");
    // ------------------------------------------------------------------
    {
        auto det = std::make_shared<HullWhite1F>(curve, a, 1e-12);
        PathGenerator gen(det, 7);
        std::vector<double> grid = { 0.0 };
        auto g = uniformGrid(10.0, 52);
        grid.insert(grid.end(), g.begin(), g.end());
        ScenarioSet scen = gen.generate(grid, 3);

        const std::size_t kLast = scen.nSteps - 1;
        check("r(10) == r0", scen.state(1, kLast), r0, 1e-8);
        check("D(0,10) == exp(-r0 * 10)", scen.df(1, kLast), std::exp(-r0 * 10.0), 1e-8);
    }

    // ------------------------------------------------------------------
    section("5. Reproducibility");
    // ------------------------------------------------------------------
    {
        std::vector<double> grid = { 0.0, 0.25, 0.5, 1.0, 3.0 };
        ScenarioSet s1 = PathGenerator(model, 123).generate(grid, 50);
        ScenarioSet s2 = PathGenerator(model, 123).generate(grid, 50);
        ScenarioSet s3 = PathGenerator(model, 456).generate(grid, 50);
        checkTrue("same seed -> identical states", s1.states == s2.states);
        checkTrue("same seed -> identical discount factors", s1.discount == s2.discount);
        checkTrue("different seed -> different states", s1.states != s3.states);
    }

    // ------------------------------------------------------------------
    section("6. Same paths as a hand-written loop (seed 42, weekly grid)");
    // ------------------------------------------------------------------
    {
        const double dt = 1.0 / 52.0;
        const int nSteps = 52 * 3;
        const std::size_t nPaths = 5;

        std::vector<double> grid;
        for (int k = 0; k <= nSteps; ++k) grid.push_back(k * dt);
        ScenarioSet scen = PathGenerator(model, 42).generate(grid, nPaths);

        std::mt19937_64 rng(42);
        std::normal_distribution<double> N01(0.0, 1.0);
        std::vector<double> z(1);
        double maxErrR = 0.0, maxErrD = 0.0;
        for (std::size_t p = 0; p < nPaths; ++p) {
            std::vector<double> s = model->initialState();
            double integral = 0.0;
            for (int k = 0; k < nSteps; ++k) {
                z[0] = N01(rng);
                std::vector<double> next = model->step(grid[k], grid[k + 1] - grid[k], s, z);
                integral += 0.5 * (s[0] + next[0]) * (grid[k + 1] - grid[k]);
                s = std::move(next);
            }
            maxErrR = std::max(maxErrR, std::abs(s[0] - scen.state(p, nSteps)));
            maxErrD = std::max(maxErrD, std::abs(std::exp(-integral) - scen.df(p, nSteps)));
        }
        check("max |r_T(manual) - r_T(generator)|", maxErrR, 0.0, 1e-15);
        check("max |D(manual) - D(generator)|", maxErrD, 0.0, 1e-15);
    }

    // ------------------------------------------------------------------
    section("7. Monte Carlo checks on the generated scenarios");
    //   E[D(0,t_k)]                  = P^M(0,t_k)
    //   E[D(0,t) P(t,T | r_t)]       = P^M(0,T)     (using stored states)
    // ------------------------------------------------------------------
    {
        std::vector<double> grid = { 0.0 };
        auto g = uniformGrid(10.0, 52);
        grid.insert(grid.end(), g.begin(), g.end());
        const std::size_t nPaths = 20000;
        ScenarioSet scen = PathGenerator(model, 42).generate(grid, nPaths);

        for (double T : { 1.0, 5.0, 10.0 }) {
            const std::size_t k = static_cast<std::size_t>(std::round(T * 52));
            double s1 = 0.0, s2 = 0.0;
            for (std::size_t p = 0; p < nPaths; ++p) { double d = scen.df(p, k); s1 += d; s2 += d * d; }
            double mean = s1 / nPaths, se = std::sqrt((s2 / nPaths - mean * mean) / nPaths);
            check("E[D(0," + std::to_string(T) + ")] == P^M", mean, curve->getDiscountFactor(T), 4.0 * se + 1e-4);
        }

        {
            const double T = 10.0;
            const std::size_t k = static_cast<std::size_t>(std::round(5.0 * 52));
            const double t = scen.times[k];
            double s1 = 0.0, s2 = 0.0;
            for (std::size_t p = 0; p < nPaths; ++p) {
                double x = scen.df(p, k) * model->zcb(t, T, scen.state(p, k));
                s1 += x; s2 += x * x;
            }
            double mean = s1 / nPaths, se = std::sqrt((s2 / nPaths - mean * mean) / nPaths);
            check("E[D(0,5) P(5,10)] == P^M(0,10)", mean, curve->getDiscountFactor(T), 4.0 * se + 1e-4);
        }
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
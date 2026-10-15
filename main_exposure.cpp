// Tests for netting-set exposure (NettingSet, nettedValues, computeExposure)
//
// Builds on: FlatCurve, HullWhite1F, ZCB, IRS, buildTimeGrid / PathGenerator / ScenarioSet,
//            ValuationEngine / MtMCube.
// Also writes exposure_irs.csv (t, EE, ENE, EEu, PFE) for a par payer swap, for plotting.
//
// Note: exclude the other test files with their own main() from the build.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "FlatCurve.hpp"
#include "HullWhite1F.hpp"
#include "Zcb.hpp"
#include "IRS.hpp"
#include "PathGenerator.hpp"
#include "ValuationEngine.hpp"
#include "Exposure.hpp"

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

    double maxAbsDiff(const std::vector<double>& x, const std::vector<double>& y) {
        double m = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) m = std::max(m, std::abs(x[i] - y[i]));
        return m;
    }

    using Trades = std::vector<std::shared_ptr<const Instrument>>;

} // namespace

int main() {
    const double r0 = 0.03;
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    auto model = std::make_shared<HullWhite1F>(curve, a, sigma);

    // Spot-starting 5y swaps: fixed annual, float quarterly
    const Schedule fixedLeg = makeTestSchedule(0.0, 5.0, 1);
    const Schedule floatLeg = makeTestSchedule(0.0, 5.0, 4);
    double annuity0 = 0.0;
    for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i)
        annuity0 += fixedLeg.taus[i - 1] * curve->getDiscountFactor(fixedLeg.dates[i]);
    const double Kpar = (1.0 - curve->getDiscountFactor(5.0)) / annuity0;

    // Cube rows
    Trades trades;
    trades.push_back(std::make_shared<ZCB>(3.0));                                          // 0: ZCB
    trades.push_back(std::make_shared<IRS>(Kpar, fixedLeg, floatLeg, 0));                  // 1: par payer swap, ids 0..19
    trades.push_back(std::make_shared<IRS>(Kpar - 0.003, fixedLeg, floatLeg, 0));          // 2: same schedule -> shares fixings

    std::vector<double> weekly;
    for (int k = 1; k <= 52 * 5; ++k) weekly.push_back(k / 52.0);
    const std::vector<double> grid = buildTimeGrid(weekly, trades, 5.0);
    auto idx = [&](double t) {
        for (std::size_t k = 0; k < grid.size(); ++k) if (std::abs(grid[k] - t) < 1e-10) return k;
        throw std::runtime_error("time not on grid");
        };

    const std::size_t nPaths = 10000;
    ScenarioSet scen = PathGenerator(model, 42).generate(grid, nPaths);
    MtMCube cube = ValuationEngine(model).run(scen, trades);
    const std::size_t nSteps = cube.nSteps;

    auto singleNS = [](std::size_t i, double q) { NettingSet ns; ns.counterparty = "C"; ns.add(i, q); return ns; };

    // ------------------------------------------------------------------
    section("1. Input validation");
    // ------------------------------------------------------------------
    {
        bool t1 = false, t2 = false, t3 = false;
        NettingSet bad; bad.add(7, 1.0);
        try { computeExposure(cube, scen, bad); }
        catch (const std::invalid_argument&) { t1 = true; }
        NettingSet mism; mism.tradeIdx = { 0, 1 }; mism.quantity = { 1.0 };
        try { computeExposure(cube, scen, mism); }
        catch (const std::invalid_argument&) { t2 = true; }
        try { computeExposure(cube, scen, singleNS(0, 1.0), 1.0); }
        catch (const std::invalid_argument&) { t3 = true; }
        checkTrue("trade index not in cube throws", t1);
        checkTrue("tradeIdx / quantity size mismatch throws", t2);
        checkTrue("PFE quantile outside (0,1) throws", t3);
    }

    // ------------------------------------------------------------------
    section("2. Single long ZCB: exposure is the whole value");
    // ------------------------------------------------------------------
    {
        ExposureProfile e = computeExposure(cube, scen, singleNS(0, 1.0));
        const double P3 = curve->getDiscountFactor(3.0);
        check("EE(0) == P^M(0,3)", e.EE[0], P3, 1e-12);   // averaging 10k equal values: rounding ~1e-13
        check("PFE(0) == P^M(0,3)", e.PFE[0], P3, 1e-15);
        double maxENE = 0.0;
        for (double x : e.ENE) maxENE = std::max(maxENE, std::abs(x));
        check("ENE == 0 everywhere", maxENE, 0.0, 0.0);
        // E[D V] = P(0,3) before maturity (martingale); se ~ 1e-4 at 10k paths
        check("EE(2) ~= P^M(0,3)", e.EE[idx(2.0)], P3, 2e-3);
        check("EE(3) == 0 (matured)", e.EE[idx(3.0)], 0.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("3. Identities");
    // ------------------------------------------------------------------
    {
        NettingSet ns = singleNS(1, 1.0);
        ExposureProfile e = computeExposure(cube, scen, ns);

        // EE + ENE = E[D V]
        std::vector<double> v = nettedValues(cube, ns);
        double maxErr = 0.0;
        for (std::size_t k = 0; k < nSteps; ++k) {
            double m = 0.0;
            for (std::size_t p = 0; p < nPaths; ++p) m += scen.df(p, k) * v[p * nSteps + k];
            m /= static_cast<double>(nPaths);
            maxErr = std::max(maxErr, std::abs(e.EE[k] + e.ENE[k] - m));
        }
        check("max_k |EE + ENE - E[D V]|", maxErr, 0.0, 1e-15);

        // Par swap at t = 0: V0 = 0 -> EE(0) = ENE(0) = PFE(0) = 0
        check("par swap EE(0) == 0", e.EE[0], 0.0, 1e-15);
        check("par swap PFE(0) == 0", e.PFE[0], 0.0, 1e-15);

        // Ordering: 0 <= EEu, EE; PFE >= median-ish positive part  (PFE95 >= EEu holds here)
        bool ordered = true;
        for (std::size_t k = 1; k + 1 < nSteps; ++k)
            ordered = ordered && e.EE[k] >= 0.0 && e.ENE[k] <= 0.0 && e.PFE[k] >= e.EEu[k];
        checkTrue("EE >= 0, ENE <= 0, PFE95 >= EEu on the grid", ordered);

        // After the last payment nothing is left
        check("EE(5) == 0", e.EE.back(), 0.0, 0.0);
        check("ENE(5) == 0", e.ENE.back(), 0.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("4. Scaling with quantity");
    // ------------------------------------------------------------------
    {
        ExposureProfile e1 = computeExposure(cube, scen, singleNS(1, 1.0));
        ExposureProfile e2 = computeExposure(cube, scen, singleNS(1, 2.0));
        std::vector<double> twice(e1.EE.size());

        std::transform(e1.EE.begin(), e1.EE.end(), twice.begin(), [](double x) { return 2.0 * x; });
        check("EE(2N) == 2 EE(N)", maxAbsDiff(e2.EE, twice), 0.0, 1e-15);
        std::transform(e1.PFE.begin(), e1.PFE.end(), twice.begin(), [](double x) { return 2.0 * x; });
        check("PFE(2N) == 2 PFE(N)", maxAbsDiff(e2.PFE, twice), 0.0, 1e-15);

        // A short position swaps EE and ENE: EE(-1) = -ENE(+1)
        ExposureProfile es = computeExposure(cube, scen, singleNS(1, -1.0));
        std::vector<double> negENE(e1.ENE.size());
        std::transform(e1.ENE.begin(), e1.ENE.end(), negENE.begin(), [](double x) { return -x; });
        check("receiver EE == - payer ENE", maxAbsDiff(es.EE, negENE), 0.0, 1e-15);
    }

    // ------------------------------------------------------------------
    section("5. Netting");
    // ------------------------------------------------------------------
    {
        // Payer and offsetting receiver in one netting set: exposure vanishes
        NettingSet hedged; hedged.counterparty = "C"; hedged.add(1, 1.0); hedged.add(1, -1.0);
        ExposureProfile eh = computeExposure(cube, scen, hedged);
        double maxEE = 0.0;
        for (double x : eh.EE) maxEE = std::max(maxEE, x);
        check("payer + receiver, netted: max_k EE == 0", maxEE, 0.0, 0.0);

        // Same two trades without netting: exposures add up and are clearly positive
        ExposureProfile ep = computeExposure(cube, scen, singleNS(1, 1.0));
        ExposureProfile er = computeExposure(cube, scen, singleNS(1, -1.0));
        const std::size_t k = idx(2.0);
        std::cout << "  t = 2: EE payer = " << ep.EE[k] << ", EE receiver = " << er.EE[k] << "\n";
        checkTrue("payer + receiver, not netted: EE(2) > 0", ep.EE[k] + er.EE[k] > 1e-3);

        // Subadditivity for two different swaps: EE_NS <= EE_1 + EE_2, strict somewhere
        NettingSet both; both.counterparty = "C"; both.add(1, 1.0); both.add(2, -1.0);
        ExposureProfile eb = computeExposure(cube, scen, both);
        ExposureProfile e2 = computeExposure(cube, scen, singleNS(2, -1.0));
        bool sub = true; double maxBenefit = 0.0;
        for (std::size_t j = 0; j < nSteps; ++j) {
            double standalone = ep.EE[j] + e2.EE[j];
            sub = sub && eb.EE[j] <= standalone + 1e-15;
            maxBenefit = std::max(maxBenefit, standalone - eb.EE[j]);
        }
        checkTrue("EE_netted <= sum of standalone EE at every t", sub);
        checkTrue("netting benefit is strictly positive somewhere", maxBenefit > 1e-4);
    }

    // ------------------------------------------------------------------
    section("6. Profile of a par payer swap (written to exposure_irs.csv)");
    // ------------------------------------------------------------------
    {
        ExposureProfile e = computeExposure(cube, scen, singleNS(1, 1.0));
        std::ofstream out("exposure_irs.csv");
        out << "t,EE,ENE,EEu,PFE95\n";
        for (std::size_t k = 0; k < nSteps; ++k)
            out << std::setprecision(10) << e.times[k] << ',' << e.EE[k] << ',' << e.ENE[k] << ','
            << e.EEu[k] << ',' << e.PFE[k] << '\n';

        auto peak = std::max_element(e.EE.begin(), e.EE.end());
        std::cout << "  peak EE = " << *peak << " at t = " << e.times[peak - e.EE.begin()] << "\n";
        for (double t : { 0.25, 1.0, 2.0, 3.0, 4.0, 4.75 })
            std::cout << "  t = " << std::setw(5) << t << "  EE = " << std::setw(12) << e.EE[idx(t)]
            << "  ENE = " << std::setw(12) << e.ENE[idx(t)]
            << "  PFE95 = " << std::setw(12) << e.PFE[idx(t)] << "\n";

        // Hump shape: peak strictly inside the life of the swap
        const double tPeak = e.times[peak - e.EE.begin()];
        checkTrue("EE peaks strictly inside (0, 5)", tPeak > 0.5 && tPeak < 4.5);

        // Sawtooth: at a fixed payment date the payer pays K and receives the float coupon;
        // with K at par the net payment is usually negative for the payer, so EE jumps up.
        const std::size_t k2 = idx(2.0);
        std::cout << "  EE just before t = 2: " << e.EE[k2 - 1] << ", at t = 2: " << e.EE[k2] << "\n";
        checkTrue("EE jumps at the fixed payment date t = 2", std::abs(e.EE[k2] - e.EE[k2 - 1]) > 1e-3);
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
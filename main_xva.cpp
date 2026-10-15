// Tests for CVA / DVA (XVA.hpp) on top of the netting-set exposure profile
//
// Builds on: FlatCurve, HullWhite1F, ZCB, IRS, buildTimeGrid / PathGenerator / ScenarioSet,
//            ValuationEngine / MtMCube, NettingSet / Exposure, CreditCurve / FlatHazardCurve.
//
// Note: exclude the other test files with their own main() from the build.


// ============================================================================
// Pipeline overview
// ============================================================================
//
//                  +----------------- zcb(t, T, r) ------------------+
//                  |                                                  |
//                  |                                                  v
//  DiscountCurve --+--> HullWhite1F --> PathGenerator --> ScenarioSet  |
//                                                           |     |   |
//                                                    r(p,t) |     | D(p,t)
//                                                           v     |   |
//  Instruments (ZCB / FRA / IRS) ----------------> ValuationEngine <--+
//                                                        |        |
//                                                        v        |
//                                                     MtMCube     |
//                                                        |        |
//                                                        v        |
//  NettingSet ------------------------------------> computeExposure <-+
//                                                        |
//                                                        v
//                                                  EE / ENE / PFE
//                                                        |
//                                                        v
//  CreditCurve + LGD -----------------------------> computeXVA
//                                                        |
//                                                        v
//                                                CVA / DVA / BCVA
//
// Stages (each output can be stored and reused):
//
//   1. Scenarios    DiscountCurve, HullWhite1F, time grid -> ScenarioSet
//                   states r(p, t_k) and discount factors D(p, t_k) = exp(-int_0^t r du).
//                   Knows nothing about trades.
//
//   2. Valuation    ScenarioSet + Instruments -> MtMCube  V_i(p, t_k)
//                   per path and grid point: record due fixings, point the market view
//                   at (t_k, r), value every trade (unit notional, undiscounted, ex-coupon).
//                   The model is used here too, through zcb(t, T, r).
//
//   3. Exposure     MtMCube + NettingSet + D(p, t_k) -> EE / ENE / PFE per t_k
//                   net across trades on each path first, then take positive / negative
//                   parts, discount and average across paths.
//
//   4. XVA          EE / ENE + CreditCurve + LGD -> CVA / DVA / BCVA
//                   sum over time buckets of exposure x default probability x LGD.
//
// What to re-run when inputs change:
//   credit spreads / LGD        -> stage 4 only
//   netting, incremental trade  -> stages 3-4 (plus stage 2 for the new trade's row)
//   rates curve / model params  -> everything
// ============================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "FlatCurve.hpp"
#include "FlatHazardCurve.hpp"
#include "HullWhite1F.hpp"
#include "Zcb.hpp"
#include "IRS.hpp"
#include "PathGenerator.hpp"
#include "ValuationEngine.hpp"
#include "Exposure.hpp"
#include "XVA.hpp"

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

    // Synthetic profile on a uniform grid 0..T with n steps: EE(t) = f(t), ENE(t) = g(t)
    template <class F, class G>
    ExposureProfile syntheticProfile(double T, int n, F f, G g) {
        ExposureProfile e;
        for (int k = 0; k <= n; ++k) {
            double t = T * k / n;
            e.times.push_back(t);
            e.EE.push_back(f(t));
            e.ENE.push_back(g(t));
            e.EEu.push_back(f(t));
            e.PFE.push_back(f(t));
        }
        return e;
    }

    Schedule makeTestSchedule(double start, double end, int paymentsPerYear) {
        const int n = static_cast<int>(std::round((end - start) * paymentsPerYear));
        Schedule s;
        for (int i = 0; i <= n; ++i) s.dates.push_back(start + i * (end - start) / n);
        for (int i = 1; i <= n; ++i) s.taus.push_back(s.dates[i] - s.dates[i - 1]);
        return s;
    }

} // namespace

int main() {
    const double lgd = 0.6;

    // ------------------------------------------------------------------
    section("1. Input validation");
    // ------------------------------------------------------------------
    {
        FlatHazardCurve c(0.02);
        ExposureProfile ok = syntheticProfile(1.0, 4, [](double) { return 1.0; }, [](double) { return 0.0; });
        bool t1 = false, t2 = false, t3 = false, t4 = false;
        try { computeCVA(ok, c, 1.5); }
        catch (const std::invalid_argument&) { t1 = true; }
        ExposureProfile bad = ok; bad.times[2] = bad.times[1];
        try { computeCVA(bad, c, lgd); }
        catch (const std::invalid_argument&) { t2 = true; }
        ExposureProfile shortEE = ok; shortEE.EE.pop_back();
        try { computeCVA(shortEE, c, lgd); }
        catch (const std::invalid_argument&) { t3 = true; }
        try { FlatHazardCurve neg(-0.01); }
        catch (const std::invalid_argument&) { t4 = true; }
        checkTrue("LGD outside [0,1] throws", t1);
        checkTrue("non-increasing times throws", t2);
        checkTrue("inconsistent vector sizes throws", t3);
        checkTrue("negative hazard rate throws", t4);
    }

    // ------------------------------------------------------------------
    section("2. Analytic cases on synthetic profiles");
    // ------------------------------------------------------------------
    {
        const double lambda = 0.02, T = 5.0;
        FlatHazardCurve c(lambda);

        // Constant EE: every bucket rule gives LGD * c * (1 - S(T)) exactly
        ExposureProfile flat = syntheticProfile(T, 20, [](double) { return 0.01; }, [](double) { return -0.004; });
        const double expected = lgd * 0.01 * (1.0 - std::exp(-lambda * T));
        for (BucketRule r : { BucketRule::Left, BucketRule::Right, BucketRule::Average })
            check("constant EE: CVA == LGD c (1 - S(T))", computeCVA(flat, c, lgd, r), expected, 1e-15);
        check("constant ENE: DVA == LGD c (1 - S(T))", computeDVA(flat, c, lgd),
            lgd * 0.004 * (1.0 - std::exp(-lambda * T)), 1e-15);

        // Zero hazard: no default, no CVA
        FlatHazardCurve none(0.0);
        check("lambda = 0: CVA == 0", computeCVA(flat, none, lgd), 0.0, 0.0);

        // Linear EE(t) = t: continuous limit LGD * int_0^T t lambda e^{-lambda t} dt
        const double exact = lgd * (1.0 - std::exp(-lambda * T) * (1.0 + lambda * T)) / lambda;
        double errPrev = 0.0;
        for (int n : { 10, 40, 160 }) {
            ExposureProfile lin = syntheticProfile(T, n, [](double t) { return t; }, [](double) { return 0.0; });
            double err = std::abs(computeCVA(lin, c, lgd) - exact);
            std::cout << "  n = " << std::setw(4) << n << "  |CVA - exact| = " << err << "\n";
            if (n > 10) checkTrue("trapezoid error shrinks ~16x when grid is 4x finer (n = " + std::to_string(n) + ")",
                err < errPrev / 10.0);
            errPrev = err;
        }

        // Left <= Average <= Right for an increasing EE
        ExposureProfile lin = syntheticProfile(T, 20, [](double t) { return t; }, [](double) { return 0.0; });
        double l = computeCVA(lin, c, lgd, BucketRule::Left);
        double m = computeCVA(lin, c, lgd, BucketRule::Average);
        double rr = computeCVA(lin, c, lgd, BucketRule::Right);
        checkTrue("increasing EE: Left < Average < Right", l < m && m < rr);
        check("Average == (Left + Right) / 2", m, 0.5 * (l + rr), 1e-15);
    }

    // ------------------------------------------------------------------
    // Simulated portfolio
    // ------------------------------------------------------------------
    const double r0 = 0.03, a = 0.05, sigma = 0.01;
    auto curve = std::make_shared<FlatCurve>(r0);
    auto model = std::make_shared<HullWhite1F>(curve, a, sigma);

    const Schedule fixedLeg = makeTestSchedule(0.0, 5.0, 1);
    const Schedule floatLeg = makeTestSchedule(0.0, 5.0, 4);
    double annuity0 = 0.0;
    for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i)
        annuity0 += fixedLeg.taus[i - 1] * curve->getDiscountFactor(fixedLeg.dates[i]);
    const double Kpar = (1.0 - curve->getDiscountFactor(5.0)) / annuity0;

    std::vector<std::shared_ptr<const Instrument>> trades;
    trades.push_back(std::make_shared<ZCB>(3.0));                                 // 0
    trades.push_back(std::make_shared<IRS>(Kpar, fixedLeg, floatLeg, 0));         // 1: par payer
    trades.push_back(std::make_shared<IRS>(Kpar - 0.003, fixedLeg, floatLeg, 0)); // 2

    std::vector<double> weekly;
    for (int k = 1; k <= 52 * 5; ++k) weekly.push_back(k / 52.0);
    const std::vector<double> grid = buildTimeGrid(weekly, trades, 5.0);

    ScenarioSet scen = PathGenerator(model, 42).generate(grid, 10000);
    MtMCube cube = ValuationEngine(model).run(scen, trades);

    auto single = [](std::size_t i, double q) { NettingSet ns; ns.add(i, q); return ns; };

    FlatHazardCurve cpty(0.02);   // ~120bp spread at 60% LGD
    FlatHazardCurve own(0.01);    // ~60bp

    // ------------------------------------------------------------------
    section("3. ZCB: CVA only on the long side, DVA only on the short side");
    // ------------------------------------------------------------------
    {
        ExposureProfile longZ = computeExposure(cube, scen, single(0, 1.0));
        ExposureProfile shortZ = computeExposure(cube, scen, single(0, -1.0));
        check("long ZCB: DVA == 0", computeDVA(longZ, own, lgd), 0.0, 0.0);
        check("short ZCB: CVA == 0", computeCVA(shortZ, cpty, lgd), 0.0, 0.0);
        // Mirror: my DVA on the short side == counterparty's CVA on the long side with my curve
        check("short DVA(own) == long CVA(own curve)",
            computeDVA(shortZ, own, lgd), computeCVA(longZ, own, lgd), 1e-15);
        // Long ZCB: EE ~ P(0,3) until 3, so CVA ~ LGD P(0,3) (1 - S(3))
        const double approx = lgd * curve->getDiscountFactor(3.0) * (1.0 - std::exp(-0.02 * 3.0));
        check("long ZCB CVA ~ LGD P(0,3) (1 - S(3))", computeCVA(longZ, cpty, lgd), approx, 2e-4);
    }

    // ------------------------------------------------------------------
    section("4. Linearity and monotonicity");
    // ------------------------------------------------------------------
    {
        ExposureProfile e1 = computeExposure(cube, scen, single(1, 1.0));
        ExposureProfile e2 = computeExposure(cube, scen, single(1, 2.0));
        const double cva1 = computeCVA(e1, cpty, lgd);
        check("notional x2 -> CVA x2", computeCVA(e2, cpty, lgd), 2.0 * cva1, 1e-15);
        check("LGD x0.5 -> CVA x0.5", computeCVA(e1, cpty, 0.3), 0.5 * cva1, 1e-15);
        checkTrue("higher hazard rate -> higher CVA",
            computeCVA(e1, FlatHazardCurve(0.04), lgd) > cva1 && cva1 > computeCVA(e1, FlatHazardCurve(0.01), lgd));
    }

    // ------------------------------------------------------------------
    section("5. Netting reduces CVA");
    // ------------------------------------------------------------------
    {
        NettingSet both; both.add(1, 1.0); both.add(2, -1.0);
        const double netted = computeCVA(computeExposure(cube, scen, both), cpty, lgd);
        const double standalone = computeCVA(computeExposure(cube, scen, single(1, 1.0)), cpty, lgd)
            + computeCVA(computeExposure(cube, scen, single(2, -1.0)), cpty, lgd);
        std::cout << "  netted CVA = " << netted << ", sum of standalone CVA = " << standalone << "\n";
        checkTrue("CVA(netted) < sum of standalone CVA", netted < standalone);

        NettingSet hedged; hedged.add(1, 1.0); hedged.add(1, -1.0);
        check("perfect hedge in one netting set: CVA == 0",
            computeCVA(computeExposure(cube, scen, hedged), cpty, lgd), 0.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("6. Par payer swap: CVA, DVA, bucket contributions");
    // ------------------------------------------------------------------
    {
        ExposureProfile e = computeExposure(cube, scen, single(1, 1.0));
        XVAResult x = computeXVA(e, cpty, lgd, own, lgd);

        double sumC = 0.0, sumD = 0.0;
        for (double v : x.cvaByBucket) sumC += v;
        for (double v : x.dvaByBucket) sumD += v;
        check("sum of CVA bucket contributions == CVA", sumC, x.cva, 1e-15);
        check("sum of DVA bucket contributions == DVA", sumD, x.dva, 1e-15);
        check("BCVA == CVA - DVA", x.bcva, x.cva - x.dva, 0.0);
        check("computeXVA CVA == computeCVA", x.cva, computeCVA(e, cpty, lgd), 0.0);

        // Running-spread view: CVA spread = CVA / annuity, should be << counterparty spread (120bp)
        const double cvaSpreadBp = x.cva / annuity0 * 1e4;
        std::cout << std::setprecision(6)
            << "  notional 1, 5y par payer, lambda_C = 2%, lambda_B = 1%, LGD = 60%\n"
            << "  CVA  = " << x.cva * 1e4 << " bp of notional (" << cvaSpreadBp << " bp running)\n"
            << "  DVA  = " << x.dva * 1e4 << " bp of notional\n"
            << "  BCVA = " << x.bcva * 1e4 << " bp of notional\n";
        checkTrue("CVA and DVA are positive", x.cva > 0.0 && x.dva > 0.0);
        checkTrue("CVA running spread well below the counterparty spread", cvaSpreadBp > 0.0 && cvaSpreadBp < 120.0);

        // Bucket rules differ only by discretisation error on a weekly grid
        const double cl = computeCVA(e, cpty, lgd, BucketRule::Left);
        const double cr = computeCVA(e, cpty, lgd, BucketRule::Right);
        std::cout << "  CVA Left / Average / Right = " << cl << " / " << x.cva << " / " << cr << "\n";
        checkTrue("bucket rules agree within 2%", std::abs(cl - cr) < 0.02 * x.cva);
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
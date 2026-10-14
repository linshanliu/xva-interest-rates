// Tests for the FRA instrument (via HW1FView + FixingStore)
//
// Assumes:
//   - HullWhite1F(std::shared_ptr<const DiscountCurve> curve, double a, double sigma)
//     with initialState(), step(t, dt, curr, z) const, zcb(t, T, r) const
//   - FixingStore(std::size_t n) with set(id, v), get(id), clear(); unset entries are NaN
//   - HW1FView(const HullWhite1F& model, const FixingStore& fixings), reset(t, r)
//   - FRA(double T1, double T2, double K, int fixingId)          (tau = T2 - T1)
//     unit notional, receive floating / pay fixed: tau * (L - K) paid at T2
//
// Simulation convention tested here: at each grid point, any fixing due at that
// time is recorded BEFORE instruments are valued.
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
#include "FixingStore.hpp"
#include "HW1FView.hpp"
#include "FRA.hpp"

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

    // Simple rate for [T1, T2] implied by a bond price P(T1, T2)
    double simpleRate(double P, double tau) {
        return (1.0 / P - 1.0) / tau;
    }

} // namespace

int main() {
    const double r0 = 0.03;
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    HullWhite1F model(curve, a, sigma);
    const double rInit = model.initialState()[0];

    // FRA on [T1, T2] = [2, 2.5]
    const double T1 = 2.0;
    const double T2 = 2.5;
    const double tau = T2 - T1;
    const int    id = 3;                 // deliberately not 0, to catch hard-coded ids
    const std::size_t nFixings = 5;

    const double P1 = curve->getDiscountFactor(T1);
    const double P2 = curve->getDiscountFactor(T2);
    const double Kpar = simpleRate(P2 / P1, tau);   // (P1/P2 - 1) / tau

    FixingStore fixings(nFixings);
    HW1FView view(model, fixings);

    // ------------------------------------------------------------------
    section("1. Contract data");
    // ------------------------------------------------------------------
    {
        FRA fra(T1, T2, 0.03, id);
        check("maturity() == T2", fra.maturity(), T2, 0.0);

        auto ev = fra.eventDates();
        checkTrue("eventDates() == {T1, T2}", ev.size() == 2 && ev[0] == T1 && ev[1] == T2);

        auto fx = fra.fixings();
        checkTrue("fixings() has one entry", fx.size() == 1);
        if (fx.size() == 1) {
            checkTrue("fixing spec: id", fx[0].id == id);
            checkTrue("fixing spec: fixDate", fx[0].fixDate == T1);
            checkTrue("fixing spec: [start, end] == [T1, T2]", fx[0].start == T1 && fx[0].end == T2);
            check("fixing spec: tau", fx[0].tau, tau, 1e-15);
        }
    }

    // ------------------------------------------------------------------
    section("2. Value at t = 0");
    // ------------------------------------------------------------------
    view.reset(0.0, rInit);
    {
        FRA par(T1, T2, Kpar, id);
        check("par FRA: V(0) == 0", par.markToMarket(view), 0.0, 1e-14);

        const double K = 0.04;
        FRA fra(T1, T2, K, id);
        check("V(0) == P(0,T1) - (1 + tau K) P(0,T2)",
            fra.markToMarket(view), P1 - (1.0 + tau * K) * P2, 1e-14);
    }

    // ------------------------------------------------------------------
    section("3. Reading a fixing before it is set gives NaN");
    // ------------------------------------------------------------------
    {
        fixings.clear();
        FRA fra(T1, T2, 0.03, id);
        view.reset(2.2, 0.03);
        checkTrue("MtM in (T1, T2) with no fixing recorded is NaN",
            std::isnan(fra.markToMarket(view)));
    }

    // ------------------------------------------------------------------
    section("4. Continuity across the fixing date T1");
    //   Same state r at T1: value just before fixing vs value just after.
    // ------------------------------------------------------------------
    {
        const double K = 0.035;
        FRA fra(T1, T2, K, id);

        for (double r : {0.01, 0.03, 0.06}) {
            fixings.clear();
            double L = simpleRate(model.zcb(T1, T2, r), tau);
            fixings.set(id, L);

            view.reset(T1, r);
            double vAfter = fra.markToMarket(view);              // fixed branch

            double vFormula = 1.0 - (1.0 + tau * K) * model.zcb(T1, T2, r);   // pre-fixing formula at T1
            check("V(T1) fixed branch == pre-fixing formula (r=" + std::to_string(r) + ")",
                vAfter, vFormula, 1e-14);

            view.reset(T1 - 1e-8, r);
            double vBefore = fra.markToMarket(view);             // unfixed branch
            check("V(T1-eps) ~= V(T1) (r=" + std::to_string(r) + ")", vBefore, vAfter, 1e-8);
        }
    }

    // ------------------------------------------------------------------
    section("5. Fixed branch and after maturity");
    // ------------------------------------------------------------------
    {
        const double K = 0.03;
        FRA fra(T1, T2, K, id);
        fixings.clear();
        const double L = 0.045;
        fixings.set(id, L);

        for (double t : {2.1, 2.4}) {
            const double r = 0.04;
            view.reset(t, r);
            check("V(" + std::to_string(t) + ") == tau (L - K) P(t,T2)",
                fra.markToMarket(view), tau * (L - K) * model.zcb(t, T2, r), 1e-14);
        }
        for (double t : {T2, 3.0}) {
            view.reset(t, 0.04);
            check("V(" + std::to_string(t) + ") == 0 after T2", fra.markToMarket(view), 0.0, 0.0);
        }
    }

    // ------------------------------------------------------------------
    section("6. Monte Carlo with fixings recorded on the path");
    //   E[ D(0,t) V(t) ] = V(0)  for t before and after the fixing date
    // ------------------------------------------------------------------
    {
        const double K = Kpar + 0.005;          // off-market, so V(0) != 0
        FRA fra(T1, T2, K, id);
        view.reset(0.0, rInit);
        const double V0 = fra.markToMarket(view);

        const double dt = 1.0 / 52.0;
        const int kA = static_cast<int>(std::round(1.0 / dt));   // t = 1.0  (before fixing)
        const int kFix = static_cast<int>(std::round(T1 / dt));   // t = 2.0  (fixing date)
        const int kB = static_cast<int>(std::round(2.25 / dt));   // t = 2.25 (after fixing)
        const int kEnd = kB;
        const int nPaths = 20000;

        std::mt19937_64 rng(42);
        std::normal_distribution<double> N01(0.0, 1.0);

        double sumA = 0.0, sumA2 = 0.0, sumB = 0.0, sumB2 = 0.0;
        double sumL = 0.0;
        std::vector<double> z(1);

        for (int p = 0; p < nPaths; ++p) {
            fixings.clear();
            std::vector<double> s = model.initialState();
            double integral = 0.0;

            for (int k = 0; k <= kEnd; ++k) {
                const double t = k * dt;
                view.reset(t, s[0]);

                // 1) record fixings due at this grid point
                if (k == kFix) {
                    double L = simpleRate(view.discount(T2), tau);   // P(T1, T2 | r_T1)
                    fixings.set(id, L);
                    sumL += L;
                }

                // 2) value
                if (k == kA) {
                    double x = std::exp(-integral) * fra.markToMarket(view);
                    sumA += x; sumA2 += x * x;
                }
                if (k == kB) {
                    double x = std::exp(-integral) * fra.markToMarket(view);
                    sumB += x; sumB2 += x * x;
                }

                // 3) evolve
                if (k < kEnd) {
                    z[0] = N01(rng);
                    std::vector<double> next = model.step(t, dt, s, z);
                    integral += 0.5 * (s[0] + next[0]) * dt;
                    s = std::move(next);
                }
            }
        }

        auto stats = [&](double s1, double s2, double& mean, double& se) {
            mean = s1 / nPaths;
            double var = s2 / nPaths - mean * mean;
            se = std::sqrt(var / nPaths);
            };
        double mA, seA, mB, seB;
        stats(sumA, sumA2, mA, seA);
        stats(sumB, sumB2, mB, seB);

        std::cout << "  V(0) = " << std::setprecision(8) << V0
            << ", mean fixing E[L] = " << sumL / nPaths
            << " (par rate " << Kpar << ")\n";

        check("before fixing: E[D(0,1) V(1)] == V(0)", mA, V0, 4.0 * seA + 5e-5);
        check("after fixing:  E[D(0,2.25) V(2.25)] == V(0)", mB, V0, 4.0 * seB + 5e-5);
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
// Tests for the ZCB instrument (via HW1FView)
//
// Assumes:
//   - HullWhite1F(std::shared_ptr<const DiscountCurve> curve, double a, double sigma)
//     with initialState(), step(t, dt, curr, z) const, zcb(t, T, r) const
//   - FixingStore with an explicit constructor FixingStore(std::size_t n)
//   - HW1FView(const HullWhite1F& model, const FixingStore& fixings), reset(t, r)
//   - ZCB(double T): unit-notional zero-coupon bond paying 1 at T
//
// Note: if this file is in the same project as the model test main.cpp, exclude one of
// them from the build (two main() functions cannot be linked together).

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
#include "Zcb.hpp"

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

    void section(const std::string& title) {
        std::cout << "\n=== " << title << " ===\n";
    }

} // namespace

int main() {
    const double r0 = 0.03;
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    HullWhite1F model(curve, a, sigma);
    const double rInit = model.initialState()[0];

    FixingStore fixings(0);            // ZCB has no fixings
    HW1FView view(model, fixings);

    // ------------------------------------------------------------------
    section("1. Contract data");
    // ------------------------------------------------------------------
    {
        ZCB bond(10.0);
        check("maturity() == 10", bond.maturity(), 10.0, 0.0);
        auto ev = bond.eventDates();
        check("eventDates().size() == 1", static_cast<double>(ev.size()), 1.0, 0.0);
        check("eventDates()[0] == 10", ev.empty() ? -1.0 : ev[0], 10.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("2. t = 0: MtM == P^M(0,T)");
    // ------------------------------------------------------------------
    view.reset(0.0, rInit);
    for (double T : {0.5, 1.0, 5.0, 10.0, 30.0}) {
        ZCB bond(T);
        check("MtM(0) for T = " + std::to_string(T),
            bond.markToMarket(view), curve->getDiscountFactor(T), 1e-12);
    }

    // ------------------------------------------------------------------
    section("3. After maturity: MtM == 0");
    // ------------------------------------------------------------------
    {
        ZCB bond(5.0);
        for (double t : {5.0, 6.0, 20.0}) {
            view.reset(t, 0.04);
            check("MtM(" + std::to_string(t) + ") for T = 5", bond.markToMarket(view), 0.0, 0.0);
        }
    }

    // ------------------------------------------------------------------
    section("4. Consistency with model.zcb(t, T, r)");
    // ------------------------------------------------------------------
    {
        ZCB bond(10.0);
        for (double t : {1.0, 4.5, 9.9}) {
            for (double r : {-0.01, 0.02, 0.06}) {
                view.reset(t, r);
                check("MtM(t=" + std::to_string(t) + ", r=" + std::to_string(r) + ")",
                    bond.markToMarket(view), model.zcb(t, 10.0, r), 1e-14);
            }
        }
    }

    // ------------------------------------------------------------------
    section("5. Polymorphism: call through base-class pointer");
    // ------------------------------------------------------------------
    {
        std::shared_ptr<const Instrument> inst = std::make_shared<ZCB>(7.0);
        view.reset(2.0, 0.035);
        check("Instrument* -> ZCB::markToMarket",
            inst->markToMarket(view), model.zcb(2.0, 7.0, 0.035), 1e-14);
        check("Instrument* -> fixings() is empty (base default)",
            static_cast<double>(inst->fixings().size()), 0.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("6. Monte Carlo martingale through the instrument");
    //   E[ D(0,t) * MtM(t) ] = P^M(0,T),   t = 5, T = 10
    // ------------------------------------------------------------------
    {
        const double T = 10.0;
        const double tm = 5.0;
        const double dt = 1.0 / 52.0;
        const int kMid = static_cast<int>(std::round(tm / dt));
        const int nPaths = 20000;

        ZCB bond(T);
        std::mt19937_64 rng(42);
        std::normal_distribution<double> N01(0.0, 1.0);

        double sum = 0.0, sum2 = 0.0;
        std::vector<double> z(1);

        for (int p = 0; p < nPaths; ++p) {
            std::vector<double> s = model.initialState();
            double integral = 0.0;
            for (int k = 0; k < kMid; ++k) {
                z[0] = N01(rng);
                std::vector<double> next = model.step(k * dt, dt, s, z);
                integral += 0.5 * (s[0] + next[0]) * dt;   // trapezoid rule
                s = std::move(next);
            }
            view.reset(kMid * dt, s[0]);
            double x = std::exp(-integral) * bond.markToMarket(view);
            sum += x;
            sum2 += x * x;
        }

        double mean = sum / nPaths;
        double var = sum2 / nPaths - mean * mean;
        double se = std::sqrt(var / nPaths);
        check("E[D(0,5) * MtM(5)] == P^M(0,10)", mean, curve->getDiscountFactor(T), 4.0 * se + 1e-4);
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
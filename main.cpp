// Tests for HullWhite1F (Euler step, state = short rate r_t)
//
// Assumes:
//   - HullWhite1F(std::shared_ptr<const DiscountCurve> curve, double a, double sigma)
//   - stateDimension(), initialState()
//   - step(double t, double dt, const std::vector<double>& curr,
//          const std::vector<double>& normalDraws) const
//   - theta(double t) const          (public; if private, delete Test 2)
//   - zcb(double t, double T, double r) const
//   - FlatCurve(double rate) derives from DiscountCurve and implements
//     getDiscountFactor / instForward
//
// Model:  dr = (theta(t) - a r) dt + sigma dW,   r(0) = f^M(0,0)
// ZCB:    P(t,T) = A(t,T) exp(-B(t,T) r),   B(t,T) = (1 - exp(-a(T-t))) / a

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "FlatCurve.hpp"
#include "HullWhite1F.hpp"

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

    double Bfun(double a, double t, double T) {
        return (1.0 - std::exp(-a * (T - t))) / a;
    }

} // namespace

int main() {
    const double r0 = 0.03;   // flat curve level
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    HullWhite1F model(curve, a, sigma);
    const double rInit = model.initialState()[0];

    // ------------------------------------------------------------------
    section("1. Initial state");
    // ------------------------------------------------------------------
    check("stateDimension == 1", static_cast<double>(model.stateDimension()), 1.0, 0.0);
    check("initialState == f^M(0,0)", rInit, r0, 1e-8);

    // ------------------------------------------------------------------
    section("2. theta(t) on a flat curve");
    //   theta(t) = a r + sigma^2/(2a) (1 - exp(-2 a t))
    // ------------------------------------------------------------------
    for (double t : {0.0, 1.0, 5.0, 20.0}) {
        double expected = a * r0 + sigma * sigma / (2.0 * a) * (1.0 - std::exp(-2.0 * a * t));
        check("theta(" + std::to_string(t) + ")", model.theta(t), expected, 1e-6);
    }

    // ------------------------------------------------------------------
    section("3. Deterministic path (sigma -> 0) stays at r0");
    // ------------------------------------------------------------------
    {
        HullWhite1F detModel(curve, a, 1e-12);
        const double dt = 1.0 / 52.0;
        const int nSteps = 30 * 52;
        std::vector<double> s = detModel.initialState();
        std::vector<double> z(1, 0.0);
        for (int k = 0; k < nSteps; ++k) s = detModel.step(k * dt, dt, s, z);
        check("r(30y) == r0", s[0], r0, 1e-6);
    }

    // ------------------------------------------------------------------
    section("4. zcb: analytic properties");
    // ------------------------------------------------------------------
    // 4a. t = 0 reproduces today's curve
    for (double T : {0.5, 1.0, 5.0, 10.0, 30.0}) {
        check("zcb(0, " + std::to_string(T) + ", r0) == P^M(0,T)",
            model.zcb(0.0, T, rInit), curve->getDiscountFactor(T), 1e-12);
    }

    // 4b. Maturity: P(T,T) = 1 for any r
    for (double r : {-0.01, 0.0, 0.03, 0.10}) {
        check("zcb(5, 5, " + std::to_string(r) + ") == 1",
            model.zcb(5.0, 5.0, r), 1.0, 1e-12);
    }

    // 4c. sigma -> 0: P(t,T | r0) = exp(-r0 (T - t))
    {
        HullWhite1F detModel(curve, a, 1e-12);
        for (double t : {1.0, 5.0}) {
            double T = t + 5.0;
            check("sigma->0: zcb(" + std::to_string(t) + ", " + std::to_string(T) + ", r0)",
                detModel.zcb(t, T, r0), std::exp(-r0 * (T - t)), 1e-10);
        }
    }

    // 4d. Sensitivity: d ln P / dr = -B(t,T)   (central finite difference)
    {
        const double t = 2.0, T = 7.0, r = 0.04, h = 1e-6;
        double dlnP = (std::log(model.zcb(t, T, r + h)) - std::log(model.zcb(t, T, r - h))) / (2.0 * h);
        check("d ln P(2,7) / dr == -B(2,7)", dlnP, -Bfun(a, t, T), 1e-6);
    }

    // 4e. Monotonicity: higher rate -> lower bond price
    {
        bool ok = model.zcb(2.0, 7.0, 0.05) < model.zcb(2.0, 7.0, 0.03);
        check("zcb decreasing in r", ok ? 1.0 : 0.0, 1.0, 0.0);
    }

    // ------------------------------------------------------------------
    section("5. Monte Carlo: distribution, curve fit, martingale");
    //   E[r_T]                         = f^M(0,T) + sigma^2/(2a^2) (1 - exp(-aT))^2
    //   Var[r_T]                       = sigma^2/(2a) (1 - exp(-2aT))
    //   E[exp(-int_0^T r)]             = P^M(0,T)
    //   E[exp(-int_0^t r) P(t,T|r_t)]  = P^M(0,T)      (martingale test)
    // ------------------------------------------------------------------
    {
        const double T = 10.0;
        const double tm = 5.0;                         // intermediate date
        const double dt = 1.0 / 52.0;
        const int nSteps = static_cast<int>(std::round(T / dt));
        const int kMid = static_cast<int>(std::round(tm / dt));
        const int nPaths = 20000;

        std::mt19937_64 rng(42);
        std::normal_distribution<double> N01(0.0, 1.0);

        double sumR = 0.0, sumR2 = 0.0;
        double sumD = 0.0, sumD2 = 0.0;
        double sumM = 0.0, sumM2 = 0.0;
        std::vector<double> z(1);

        for (int p = 0; p < nPaths; ++p) {
            std::vector<double> s = model.initialState();
            double integral = 0.0;
            double mart = 0.0;

            for (int k = 0; k < nSteps; ++k) {
                if (k == kMid) {
                    // D(0,tm) * P(tm, T | r_tm)
                    mart = std::exp(-integral) * model.zcb(tm, T, s[0]);
                }
                z[0] = N01(rng);
                std::vector<double> next = model.step(k * dt, dt, s, z);
                integral += 0.5 * (s[0] + next[0]) * dt;   // trapezoid rule
                s = std::move(next);
            }

            double rT = s[0];
            double D = std::exp(-integral);
            sumR += rT;   sumR2 += rT * rT;
            sumD += D;    sumD2 += D * D;
            sumM += mart; sumM2 += mart * mart;
        }

        auto mean = [&](double s1) { return s1 / nPaths; };
        auto var = [&](double s1, double s2) { double m = s1 / nPaths; return s2 / nPaths - m * m; };

        double meanR = mean(sumR), varR = var(sumR, sumR2);
        double meanD = mean(sumD), varD = var(sumD, sumD2);
        double meanM = mean(sumM), varM = var(sumM, sumM2);

        double g = Bfun(a, 0.0, T);
        double expMeanR = r0 + 0.5 * sigma * sigma * g * g;
        double expVarR = sigma * sigma / (2.0 * a) * (1.0 - std::exp(-2.0 * a * T));
        double PT = curve->getDiscountFactor(T);

        double seMeanR = std::sqrt(expVarR / nPaths);
        double seVarR = expVarR * std::sqrt(2.0 / (nPaths - 1));
        double seD = std::sqrt(varD / nPaths);
        double seM = std::sqrt(varM / nPaths);

        // Tolerance: 4 standard errors, plus a small allowance for Euler / trapezoid bias
        check("E[r_T]", meanR, expMeanR, 4.0 * seMeanR + 1e-4);
        check("Var[r_T]", varR, expVarR, 4.0 * seVarR);
        check("E[exp(-int_0^T r)] == P^M(0,T)", meanD, PT, 4.0 * seD + 1e-4);
        check("Martingale: E[D(0,5) P(5,10)] == P^M(0,10)",
            meanM, PT, 4.0 * seM + 1e-4);
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
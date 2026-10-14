// Tests for the IRS instrument (via HW1FView + FixingStore)
//
// Assumes:
//   - HullWhite1F(std::shared_ptr<const DiscountCurve> curve, double a, double sigma)
//     with initialState(), step(t, dt, curr, z) const, zcb(t, T, r) const
//   - FixingStore(std::size_t n) with set(id, v), get(id), clear(); unset entries are NaN
//   - HW1FView(const HullWhite1F& model, const FixingStore& fixings), reset(t, r)
//   - Schedule { std::vector<double> dates; std::vector<double> taus; }   (in IRS.hpp)
//   - IRS(double K, Schedule fixedLeg, Schedule floatLeg, int firstFixingId)
//     unit notional, payer: pay fixed, receive floating
//
// Swap tested: forward-starting, start 1y, end 5y, fixed annual, float quarterly.
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

    // Regular schedule: n = round((end - start) * freq) periods, tau = date difference
    Schedule makeTestSchedule(double start, double end, int paymentsPerYear) {
        const int n = static_cast<int>(std::round((end - start) * paymentsPerYear));
        Schedule s;
        s.dates.reserve(n + 1);
        s.taus.reserve(n);
        for (int i = 0; i <= n; ++i) s.dates.push_back(start + i * (end - start) / n);
        for (int i = 1; i <= n; ++i) s.taus.push_back(s.dates[i] - s.dates[i - 1]);
        return s;
    }

    double simpleRate(double P, double tau) { return (1.0 / P - 1.0) / tau; }

} // namespace

int main() {
    const double r0 = 0.03;
    const double a = 0.05;
    const double sigma = 0.01;

    auto curve = std::make_shared<FlatCurve>(r0);
    HullWhite1F model(curve, a, sigma);
    const double rInit = model.initialState()[0];

    const double start = 1.0, end = 5.0;
    const Schedule fixedLeg = makeTestSchedule(start, end, 1);   // annual
    const Schedule floatLeg = makeTestSchedule(start, end, 4);   // quarterly
    const int firstId = 5;                                       // deliberately not 0
    const std::size_t nFloat = floatLeg.taus.size();             // 16
    const std::size_t storeSize = firstId + nFloat + 4;

    // Par rate at t = 0
    double annuity0 = 0.0;
    for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i)
        annuity0 += fixedLeg.taus[i - 1] * curve->getDiscountFactor(fixedLeg.dates[i]);
    const double Kpar = (curve->getDiscountFactor(start) - curve->getDiscountFactor(end)) / annuity0;

    FixingStore fixings(storeSize);
    HW1FView view(model, fixings);

    // ------------------------------------------------------------------
    section("1. Contract data");
    // ------------------------------------------------------------------
    {
        IRS swap(0.03, fixedLeg, floatLeg, firstId);
        check("maturity() == 5", swap.maturity(), end, 0.0);

        auto ev = swap.eventDates();
        checkTrue("eventDates(): 17 unique dates (quarterly grid 1..5)", ev.size() == 17);
        bool sorted = true;
        for (std::size_t i = 1; i < ev.size(); ++i) sorted = sorted && (ev[i] > ev[i - 1]);
        checkTrue("eventDates(): strictly increasing", sorted);

        auto fx = swap.fixings();
        checkTrue("fixings(): one per floating period (16)", fx.size() == nFloat);
        bool idsOk = true, datesOk = true;
        for (std::size_t j = 1; j <= fx.size(); ++j) {
            const FixingSpec& f = fx[j - 1];
            idsOk = idsOk && (f.id == firstId + static_cast<int>(j) - 1);
            datesOk = datesOk && (f.fixDate == floatLeg.dates[j - 1]) && (f.start == floatLeg.dates[j - 1])
                && (f.end == floatLeg.dates[j]) && (f.tau == floatLeg.taus[j - 1]);
        }
        checkTrue("fixings(): ids are firstId, firstId+1, ...", idsOk);
        checkTrue("fixings(): fixDate/start/end/tau match the float schedule", datesOk);
    }

    // ------------------------------------------------------------------
    section("2. Value at t = 0");
    // ------------------------------------------------------------------
    view.reset(0.0, rInit);
    {
        IRS par(Kpar, fixedLeg, floatLeg, firstId);
        check("par swap: V(0) == 0", par.markToMarket(view), 0.0, 1e-14);

        const double K = 0.04;
        IRS swap(K, fixedLeg, floatLeg, firstId);
        double expected = curve->getDiscountFactor(start) - curve->getDiscountFactor(end) - K * annuity0;
        check("V(0) == P(0,T0) - P(0,Tn) - K * A(0)", swap.markToMarket(view), expected, 1e-14);
    }

    // ------------------------------------------------------------------
    section("3. Reading a fixing before it is set gives NaN");
    // ------------------------------------------------------------------
    {
        fixings.clear();
        IRS swap(0.03, fixedLeg, floatLeg, firstId);
        view.reset(1.1, 0.03);
        checkTrue("MtM after T0 with no fixing recorded is NaN", std::isnan(swap.markToMarket(view)));
    }

    // ------------------------------------------------------------------
    section("4. Telescoped formula vs period-by-period sum");
    //   Floating leg at t in (T_{j*-1}, T_{j*}):
    //   tau L P(t,T_{j*}) + sum_{j > j*} [P(t,T_{j-1}) - P(t,T_j)]
    // ------------------------------------------------------------------
    {
        const double K = 0.035;
        IRS swap(K, fixedLeg, floatLeg, firstId);
        const double t = 2.6, r = 0.04;                 // period [2.5, 2.75] is current

        fixings.clear();
        for (std::size_t j = 1; j <= nFloat; ++j)       // set every fixing to a distinct value
            fixings.set(firstId + static_cast<int>(j) - 1, 0.02 + 0.001 * static_cast<double>(j));

        // brute force
        const auto& T = floatLeg.dates;
        double pvFloat = 0.0;
        for (std::size_t j = 1; j <= nFloat; ++j) {
            if (T[j] <= t) continue;                                   // already paid
            if (T[j - 1] <= t) {                                       // current period: fixed
                double L = fixings.get(firstId + static_cast<int>(j) - 1);
                pvFloat += floatLeg.taus[j - 1] * L * model.zcb(t, T[j], r);
            }
            else {                                                   // future period
                pvFloat += model.zcb(t, T[j - 1], r) - model.zcb(t, T[j], r);
            }
        }
        double pvFixed = 0.0;
        for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i)
            if (fixedLeg.dates[i] > t) pvFixed += K * fixedLeg.taus[i - 1] * model.zcb(t, fixedLeg.dates[i], r);

        view.reset(t, r);
        check("V(2.6) == brute-force sum over periods", swap.markToMarket(view), pvFloat - pvFixed, 1e-13);
    }

    // ------------------------------------------------------------------
    section("5. Behaviour at event dates");
    // ------------------------------------------------------------------
    {
        const double K = 0.035;
        IRS swap(K, fixedLeg, floatLeg, firstId);
        const double eps = 1e-8;
        const double r = 0.04;
        const auto& T = floatLeg.dates;

        // 5a. Start date T0 (fixing only, no payment): continuous
        {
            fixings.clear();
            fixings.set(firstId, simpleRate(model.zcb(T[0], T[1], r), floatLeg.taus[0]));
            view.reset(T[0] - eps, r);
            double vBefore = swap.markToMarket(view);
            view.reset(T[0], r);
            double vAfter = swap.markToMarket(view);
            check("T0 = 1.0: continuous across first fixing", vBefore, vAfter, 1e-7);
        }

        // 5b. Payment dates: V(T_j - eps) - V(T_j) == net payment at T_j
        //     float pays tau_j L_j; fixed pays alpha K if T_j is also a fixed payment date
        for (std::size_t j : {2u, 4u, 7u}) {             // T_2 = 1.5, T_4 = 2.0 (also fixed), T_7 = 2.75
            const double Tj = T[j];
            const double Lj = 0.045;                     // fixing of the period ending at T_j
            fixings.clear();
            fixings.set(firstId + static_cast<int>(j) - 1, Lj);
            fixings.set(firstId + static_cast<int>(j),       // next period, fixed at T_j from state r
                simpleRate(model.zcb(Tj, T[j + 1], r), floatLeg.taus[j]));

            view.reset(Tj - eps, r);
            double vBefore = swap.markToMarket(view);
            view.reset(Tj, r);
            double vAfter = swap.markToMarket(view);

            double netPayment = floatLeg.taus[j - 1] * Lj;
            for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i)
                if (fixedLeg.dates[i] == Tj) netPayment -= K * fixedLeg.taus[i - 1];

            check("jump at T = " + std::to_string(Tj) + " equals net payment",
                vBefore - vAfter, netPayment, 1e-7);
        }

        // 5c. After maturity
        for (double t : {end, 6.0}) {
            view.reset(t, r);
            check("V(" + std::to_string(t) + ") == 0 after maturity", swap.markToMarket(view), 0.0, 0.0);
        }
    }

    // ------------------------------------------------------------------
    section("6. Monte Carlo with fixings recorded on the path");
    //   Martingale for the discounted GAINS process (value + cash flows already paid):
    //   E[ D(0,t) V(t) + sum_{payment dates <= t} D(0,T_j) CF_j ] = V(0)
    //   V(t) excludes paid cash flows (ex-coupon), so E[D(0,t) V(t)] alone != V(0)
    //   once payments have occurred.
    // ------------------------------------------------------------------
    {
        const double K = Kpar + 0.005;                   // off-market, so V(0) != 0
        IRS swap(K, fixedLeg, floatLeg, firstId);
        view.reset(0.0, rInit);
        const double V0 = swap.markToMarket(view);

        const double dt = 1.0 / 52.0;                    // quarterly dates land on the grid (13 steps)
        const std::vector<int> kObs = { 26, 135, 254 };  // t ~ 0.5 (before start), 2.60, 4.88
        const int kEnd = kObs.back();
        const int nPaths = 20000;

        // Map each fixing to the grid step at which it must be recorded
        std::vector<std::vector<FixingSpec>> fixingsAtStep(kEnd + 1);
        for (const FixingSpec& f : swap.fixings()) {
            int k = static_cast<int>(std::round(f.fixDate / dt));
            if (std::abs(k * dt - f.fixDate) > 1e-9) {
                std::cout << "Fixing date " << f.fixDate << " is not on the grid\n";
                return 1;
            }
            if (k <= kEnd) fixingsAtStep[k].push_back(f);
        }

        // Map each payment to the grid step at which it is paid.
        // Floating payment at T_j: + tau_j * L_j (needs fixing id of period j)
        // Fixed payment at S_i:    - alpha_i * K
        struct FloatPay { int fixingId; double tau; };
        std::vector<std::vector<FloatPay>> floatPayAtStep(kEnd + 1);
        std::vector<double> fixedPayAtStep(kEnd + 1, 0.0);
        for (std::size_t j = 1; j < floatLeg.dates.size(); ++j) {
            int k = static_cast<int>(std::round(floatLeg.dates[j] / dt));
            if (k <= kEnd)
                floatPayAtStep[k].push_back({ firstId + static_cast<int>(j) - 1, floatLeg.taus[j - 1] });
        }
        for (std::size_t i = 1; i < fixedLeg.dates.size(); ++i) {
            int k = static_cast<int>(std::round(fixedLeg.dates[i] / dt));
            if (k <= kEnd) fixedPayAtStep[k] += K * fixedLeg.taus[i - 1];
        }

        std::mt19937_64 rng(42);
        std::normal_distribution<double> N01(0.0, 1.0);

        std::vector<double> sum(kObs.size(), 0.0), sum2(kObs.size(), 0.0);
        std::vector<double> z(1);

        for (int p = 0; p < nPaths; ++p) {
            fixings.clear();
            std::vector<double> s = model.initialState();
            double integral = 0.0;
            double paidPV = 0.0;   // sum of D(0, T_pay) * CF over payments already made

            for (int k = 0; k <= kEnd; ++k) {
                const double t = k * dt;
                const double D = std::exp(-integral);   // D(0, t_k)

                // 0) cash flows paid at this step (uses fixings recorded at earlier steps)
                double cf = -fixedPayAtStep[k];
                for (const FloatPay& fp : floatPayAtStep[k]) cf += fp.tau * fixings.get(fp.fixingId);
                paidPV += D * cf;

                // 1) record fixings due at this step (use the exact fixing date)
                for (const FixingSpec& f : fixingsAtStep[k]) {
                    view.reset(f.fixDate, s[0]);
                    fixings.set(f.id, simpleRate(view.discount(f.end), f.tau));
                }

                // 2) value at observation dates: discounted value + discounted cash flows already paid
                for (std::size_t o = 0; o < kObs.size(); ++o) {
                    if (k == kObs[o]) {
                        view.reset(t, s[0]);
                        double x = D * swap.markToMarket(view) + paidPV;
                        sum[o] += x; sum2[o] += x * x;
                    }
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

        std::cout << "  V(0) = " << std::setprecision(8) << V0 << ", par rate = " << Kpar << "\n";
        for (std::size_t o = 0; o < kObs.size(); ++o) {
            double mean = sum[o] / nPaths;
            double var = sum2[o] / nPaths - mean * mean;
            double se = std::sqrt(var / nPaths);
            check("E[D(0,t) V(t) + paid CF] == V(0) at t = " + std::to_string(kObs[o] * dt),
                mean, V0, 4.0 * se + 1e-4);
        }
    }

    std::cout << "\n" << (g_tests - g_failures) << " / " << g_tests << " tests passed.\n";
    return g_failures == 0 ? 0 : 1;
}
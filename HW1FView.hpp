#ifndef HW1FVIEW_HPP
#define HW1FVIEW_HPP
#include "MarketView.hpp"
#include "HullWhite1F.hpp"
#include "FixingStore.hpp"


class HW1FView : public MarketView {
private:
    const HullWhite1F* model_;
    const FixingStore* fix_;
    double t_ = 0.0;
    double r_ = 0.0;
public:
    HW1FView(const HullWhite1F& m, const FixingStore& f) : model_(&m), fix_(&f) {}
    void reset(double t, double r) { t_ = t; r_ = r; }

    double time() const override { return t_; }
    double discount(double T) const override { return T <= t_ ? 1.0 : model_->zcb(t_, T, r_); }
    double fixing(int id) const override { return fix_->get(id); }
};

#endif // !HW1FVIEW_HPP

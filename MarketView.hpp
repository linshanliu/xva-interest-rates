#ifndef MARKETVIEW_HPP
#define MARKETVIEW_HPP

class MarketView {
public:
    virtual ~MarketView() = default;
    virtual double time() const = 0;                 // 当前时刻 t
    virtual double discount(double T) const = 0;     // P(t,T)，T <= t 时返回 1
    virtual double fixing(int fixingId) const = 0;   // 已发生的定盘值 L_j
};


#endif // !MARKETVIEW_HPP

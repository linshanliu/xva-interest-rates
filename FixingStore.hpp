#ifndef FIXINGSTORE_HPP
#define FIXINGSTORE_HPP

#include <limits>
#include <vector>


// FixingStore only records the floating rates L_j fixed on fixing dates
// that have already passed on the current path.

class FixingStore {
public:
    explicit FixingStore(std::size_t n)
        : values_(n, std::numeric_limits<double>::quiet_NaN()) {
    }

    void set(int id, double value) { values_[id] = value; }
    double get(int id) const { return values_[id]; }

    void clear() {
        std::fill(values_.begin(), values_.end(),
            std::numeric_limits<double>::quiet_NaN());
    }

private:
    std::vector<double> values_;
};

#endif // !FIXINGSTORE_HPP

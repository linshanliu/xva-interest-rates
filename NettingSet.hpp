#ifndef NETTINGSET_HPP
#define NETTINGSET_HPP

#include <cstddef>
#include <string>
#include <vector>

// Trades with one counterparty that are settled net on default.
//
// Each entry refers to a row of the MtMCube (trade index) and carries a signed quantity:
// the cube holds unit-notional values, so the position value is quantity * V_i.
//   quantity > 0 : long the instrument as defined (e.g. payer swap, long bond)
//   quantity < 0 : short (e.g. receiver swap built from a payer IRS)
// The same cube row may appear in several netting sets (or twice in one, with opposite signs).
struct NettingSet {
    std::string counterparty;
    std::vector<std::size_t> tradeIdx;   // rows of the MtMCube
    std::vector<double> quantity;        // same length as tradeIdx

    void add(std::size_t idx, double q) {
        tradeIdx.push_back(idx);
        quantity.push_back(q);
    }
};

#endif // !NETTINGSET_HPP
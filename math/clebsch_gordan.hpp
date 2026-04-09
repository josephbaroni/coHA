#pragma once
#include <vector>
#include <cassert>
#include <map>
#include <tuple>
#include <cmath>

namespace math {

    // Single CG coefficient entry
struct CGEntry {
    int m1;
    int m2;
    int m;
    float value;
};

// CG table for (l1, l2) -> l
struct ClebschGordanTable {
    int l1;
    int l2;
    int l;

    std::vector<CGEntry> entries;

    ClebschGordanTable(int l1_, int l2_, int l_) : l1(l1_), l2(l2_), l(l_) {
        assert(std::abs(l1 - l2) <= l && l <= l1 + l2);
    }

    [[nodiscard]] int dim1() const {
        return 2 * l1 + 1;
    }

    [[nodiscard]] int dim2() const {
        return 2 * l2 + 1;
    }

    [[nodiscard]] int dim() const {
        return 2 * l + 1;
    }

};

// Container for all CGs up to some l_max
class ClebschGordan {
public:
    explicit ClebschGordan(int l_max);

    const ClebschGordanTable&
    table(int l1, int l2, int l) const;

private:
    int l_max_;
    std::map<std::tuple<int,int,int>, ClebschGordanTable> tables_;
};

}
#pragma once
#include <vector>
#include <cassert>
#include <cstdint>

namespace math {

struct Irrep {
    int l; // degree

    explicit Irrep(int l_) : l(l_) {
        assert(l >= 0);
    }

    [[nodiscard]] int dim() const {
        return 2 * l + 1;
    }

    bool operator==(const Irrep& other) const {
        return l == other.l;
    }

};

struct IrrepMultiplicity {
    Irrep irrep;
    int multiplicity;

    IrrepMultiplicity(Irrep irrep_, int mult_) : irrep(irrep_), multiplicity(mult_) {
        assert(mult_ > 0);        
    }
};

struct IrrepCollection {
    std::vector<IrrepMultiplicity> blocks;

    [[nodiscard]] int total_dim() const {
        int sum = 0;
        for (const auto& b : blocks) {
            sum += b.multiplicity * b.irrep.dim();
        }
        return sum;
    }
};

}


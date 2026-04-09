#pragma once
#include <vector>
#include <cassert>

namespace math {

    // Real Wigner D matrix for a single l
    struct WignerDReal {
        int l;
        std::vector<float> data;        // size (2l+1)^2

        WignerDReal(int l_) : l(l_), data((2*l_ + 1)*(2*l_+1), 0.f) {}

        float& operator()(int m, int mp) {
            int dim = 2 * l + 1;
            return data[(m + l) * dim + (mp + l)];
        }

        float operator()(int m, int mp) const {
            int dim = 2 * l + 1;
            return data[(m + l) * dim + (mp + l)];
        }
    };

    //Compute real Wigner D using complex-to-real change of basis
    WignerDReal real_wigner_d(int l, float alpha, float beta, float gamma);

}
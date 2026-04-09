#pragma once
#include <vector>
#include <complex>
#include <cassert>

namespace math {

// Wigner D matrix for a single l
struct WignerD {
    int l;
    std::vector<std::complex<float>> data;    // size (2l+1)^2

    explicit WignerD(int l_) : l(l_), data((2*l_+1)*(2*l_+1), {0.f, 0.f}) {}

    std::complex<float>& operator()(int m, int mp) {
        int dim = 2 * l + 1;
        return data[(m+l) * dim + (mp + l)];        
    }

    std::complex<float> operator()(int m, int mp) const {
        int dim = 2 * l + 1;
        return data[(m + l) * dim + (mp + l)];
    }
};

// Compute complex Wigner D using ZYZ Euler angles
WignerD wigner_d_matrix(int l, float alpha, float beta, float gamma);

}
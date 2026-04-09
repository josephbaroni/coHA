#include "real_basis.hpp"
#include <cmath>

namespace math {
    
    std::vector<std::vector<std::complex<float>>> real_basis_change(int l) {

        int dim = 2*l +1;
        std::vector<std::vector<std::complex<float>>> U(
            dim, std::vector<std::complex<float>>(dim, {0.f, 0.f})
        );

        auto ic = [&](int m) { return m + l; };

        // m = 0
        U[0][ic(0)] = {1.f, 0.f};

        for (int m = 1; m <=l; ++m) {
            int row_c = 2*m - 1;
            int row_s = 2*m;

            float s = (m % 2 == 0) ? 1.f : -1.f;

            // cosine part
            U[row_c][ic( m)] = { 1.f/std::sqrt(2.f), 0.f };
            U[row_c][ic(-m)] = { s/std::sqrt(2.f), 0.f };

            // sine part
            U[row_s][ic( m)] = { 0.f, -1.f/std::sqrt(2.f) };
            U[row_s][ic(-m)] = { 0.f, s/std::sqrt(2.f) };
        }

        return U;
    }
}

// This matrix is unitary, real-valued after application,
#include "wigner_real.hpp"
#include "wigner.hpp"       // for complex D
#include "real_basis.hpp"   // for U_l
#include <complex>

namespace math {

    WignerDReal real_wigner_d(int l, float alpha, float beta, float gamma) {

        auto Dc = wigner_d_matrix(l, alpha, beta, gamma);       // complex D
        auto U = real_basis_change(l);                          // complex U_l

        int dim = 2*l + 1;
        WignerDReal Dr(l);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                std::complex<float> sum = 0.f;
                for (int a = 0; a < dim; ++a) {
                    for (int b = 0; b < dim; ++b) {
                        sum += U[i][a] * Dc(a-l, b-l) * std::conj(U[j][b]);
                    }
                }
                Dr(i, j) = sum.real();
            }
        }
        return Dr;
    }

}
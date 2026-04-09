#include "wigner.hpp"
#include <cmath>

namespace math {

namespace {

inline double factorial(int n) {
    static double fact[128] = {0.0};
    static bool init = false;
    if (!init) {
        fact[0] = 1.0;
        for (int i = 1; i < 128; ++i)
            fact[i] = i * fact[i-1];
        init = true;
    }
    return fact[n];
}

// Wigner small-d element d^l_{m,m'}(beta)
double wigner_small_d(int l, int m, int mp, double beta) {
    double sum = 0.0;

    int k_min = std::max(0, m - mp);
    int k_max = std::min(l + m, l - mp);

    double prefactor =
        std::sqrt(
            factorial(l+m) * factorial(l-m) *
            factorial(l+mp) * factorial(l-mp)
        );

    for (int k = k_min; k <= k_max; ++k) {
        double denom =
            factorial(l + m - k) *
            factorial(k) *
            factorial(mp - m + k) *
            factorial(l - mp - k);

        double term = prefactor / denom;
        term *= std::pow(-1.0, k + mp - m);
        term *= std::pow(std::cos(beta/2), 2*l + m - mp - 2*k);
        term *= std::pow(std::sin(beta/2), mp - m + 2*k);

        sum += term;
    }
    return sum;
}

}

WignerD wigner_d_matrix(
    int l,
    float alpha,
    float beta,
    float gamma) {

    WignerD D(l);

    for (int m = -l; m <= l; ++m) {
        for (int mp = -l; mp <= l; ++mp) {
            double d = wigner_small_d(l, m, mp, beta);

            std::complex<float> phase = 
                std::exp(std::complex<float>(0.f, -m * alpha)) * 
                std::exp(std::complex<float>(0.f, -mp * gamma));
                
            D(m, mp) = static_cast<float>(d) * phase;
        }
    }
    return D;
}

}

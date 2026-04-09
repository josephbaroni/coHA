#include "spherical_harmonics.hpp"
#include <cmath>
#include <cassert>

namespace math {

namespace {

// factorial cache (small l only)
inline float factorial(int n) {
    static float fact[64] = {0};
    static bool init = false;
    if (!init) {
        fact[0] = 1.f;
        for (int i = 1; i < 64; ++i)
            fact[i] = i * fact[i-1];
        init = true;
    }
    return fact[n];
}

// normalization constant
inline float K(int l, int m) {
    float num = (2.f*l + 1.f) * factorial(l - m);
    float den = 4.f * float(M_PI) * factorial(l + m);
    return std::sqrt(num / den);
}

} // anonymous namespace

SphericalHarmonics::SphericalHarmonics(int l_max)
    : l_max_(l_max) {
    assert(l_max_ >= 0);
}

int SphericalHarmonics::output_dim() const {
    int dim = 0;
    for (int l = 0; l <= l_max_; ++l)
        dim += 2*l + 1;
    return dim;
}

// Associated Legendre P_l^m(x), m >= 0
float SphericalHarmonics::legendre(int l, int m, float x) {
    assert(m >= 0 && m <= l);

    float pmm = 1.f;
    if (m > 0) {
        float somx2 = std::sqrt((1.f - x) * (1.f + x));
        float fact = 1.f;
        for (int i = 1; i <= m; ++i) {
            pmm *= -fact * somx2;
            fact += 2.f;
        }
    }

    if (l == m)
        return pmm;

    float pmmp1 = x * (2*m + 1) * pmm;
    if (l == m + 1)
        return pmmp1;

    float pll = 0.f;
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = ((2*ll - 1)*x*pmmp1 - (ll + m - 1)*pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

void SphericalHarmonics::evaluate(
    float x, float y, float z,
    std::vector<float>& out) const {

    float r = std::sqrt(x*x + y*y + z*z);
    assert(r > 0.f);
    x /= r; y /= r; z /= r;

    float theta = std::acos(z);
    float phi = std::atan2(y, x);

    out.resize(output_dim());

    int offset = 0;
    for (int l = 0; l <= l_max_; ++l) {
        for (int m = -l; m <= l; ++m) {
            float value;
            if (m == 0) {
                value = K(l,0) * legendre(l,0,std::cos(theta));
            } else if (m > 0) {
                value = std::sqrt(2.f) * K(l,m)
                      * std::cos(m*phi)
                      * legendre(l,m,std::cos(theta));
            } else {
                value = std::sqrt(2.f) * K(l,-m)
                      * std::sin(-m*phi)
                      * legendre(l,-m,std::cos(theta));
            }
            out[offset++] = value;
        }
    }
}

}

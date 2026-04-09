#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <complex>

#include "../math/clebsch_gordan.hpp"
#include "../math/wigner.hpp"
#include "../math/rotation.hpp"
#include "../math/wigner_real.hpp"

using namespace math;

// ---------------- utilities ----------------

std::vector<std::complex<float>> random_vector(int dim) {
    static std::mt19937 gen(1234);
    std::normal_distribution<float> dist(0.f, 1.f);

    std::vector<std::complex<float>> v(dim);
    for (auto& x : v) {
        x = {dist(gen), 0.f};       // real-valued initial state
    }
    return v;
}

float l2_norm(const std::vector<std::complex<float>>& v) {
    float s = 0.f;
    for (auto x : v) {
        s += std::norm(x);      // |x|^2
    }
    return std::sqrt(s);
}

std::vector<std::complex<float>> apply_wigner(
    const math::WignerD& D,
    const std::vector<std::complex<float>>& v) {

    int l = D.l;
    int dim = 2*l + 1;
    std::cerr << "apply_wigner: D.l=" << D.l << ", v.size()=" << v.size() << std::endl;
    assert((int)v.size() == dim);

    std::vector<std::complex<float>> out(dim, {0.f, 0.f});

    for (int m = -l; m <= l; ++m) {
        for (int mp = -l; mp <= l; ++mp) {
            out[m + l] += D(m, mp) * v[mp + l];
        }
    }
    return out;
}

// CG contraction z = x \otimes y -> l
std::vector<std::complex<float>> cg_contract(
    const std::vector<std::complex<float>>& x,
    const std::vector<std::complex<float>>& y,
    const ClebschGordanTable& table) {

    // output irrep l
    const int l  = table.l;
    const int l1 = table.l1;
    const int l2 = table.l2;

    assert((int)x.size() == 2*l1 + 1);
    assert((int)y.size() == 2*l2 + 1);

    std::vector<std::complex<float>> z(2*l + 1, {0.f, 0.f});

    for (const auto& e : table.entries) {
        // e.m in [-l, l]
        // e.m1 in [-l1, l1]
        // e.m2 in [-l2, l2]
        z[e.m + l] +=
            e.value *
            x[e.m1 + l1] *
            y[e.m2 + l2];
    }

    return z;
}


// ---------------- test ----------------

void test_equivariance(
    int l1, int l2, int l,
    const ClebschGordan& cg) {

    const auto& table = cg.table(l1, l2, l);

    auto x = random_vector(2*l1 + 1);
    auto y = random_vector(2*l2 + 1);

    // Random rotation
    float alpha = 0.7f;
    float beta  = 1.1f;
    float gamma = -0.4f;

    auto D1 = wigner_d_matrix(l1, alpha, beta, gamma);
    auto D2 = wigner_d_matrix(l2, alpha, beta, gamma);
    auto D  = wigner_d_matrix(l,  alpha, beta, gamma);

    // Left-hand side
    auto x_rot = apply_wigner(D1, x);
    auto y_rot = apply_wigner(D2, y);
    auto z_lhs = cg_contract(x_rot, y_rot, table);

    // Right-hand side
    auto z = cg_contract(x, y, table);
    auto z_rhs = apply_wigner(D, z);

    // Compare
    std::vector<std::complex<float>> diff(z.size());
    for (size_t i = 0; i < z.size(); ++i)
        diff[i] = z_lhs[i] - z_rhs[i];

    float err = l2_norm(diff);
    float ref = l2_norm(z_rhs);

    std::cout
        << "Equivariance test "
        << "(l1=" << l1
        << ", l2=" << l2
        << ", l=" << l
        << "): error = " << err
        << std::endl;

    assert(err < 1e-4f * (ref + 1e-6f));
}

int main() {
    int l_max = 4;
    ClebschGordan cg(l_max);

    for (int l1 = 0; l1 <= l_max; ++l1) {
        for (int l2 = 0; l2 <= l_max; ++l2) {
            for (int l = std::abs(l1 - l2); l <= l1 + l2; ++l) {
                test_equivariance(l1, l2, l, cg);
            }
        }
    }

    std::cout << "All equivariance tests passed" << std::endl;
    return 0;
}

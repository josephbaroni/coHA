/*
Conventions:
- Integer l only
- m, m1, m2 in Z
- real basis CGs (signs handled explicity)
- m=m1+m2 enforced
- triangle inequality enforced
*/
#include "clebsch_gordan.hpp"
#include <cmath>
#include <algorithm>

namespace math {

namespace {

// ---------- factorial ----------
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

// ---------- Wigner 3j ----------
// (l1 l2 l)
// (m1 m2 m3)
double wigner_3j(
    int l1, int l2, int l,
    int m1, int m2, int m3) {

    if (m1 + m2 + m3 != 0) return 0.0;
    if (std::abs(m1) > l1 || std::abs(m2) > l2 || std::abs(m3) > l)
        return 0.0;
    if (l < std::abs(l1 - l2) || l > l1 + l2)
        return 0.0;

    double prefactor =
        std::pow(-1.0, l1 - l2 - m3) *
        std::sqrt(
            factorial(l1 + l2 - l) *
            factorial(l1 - l2 + l) *
            factorial(-l1 + l2 + l) /
            factorial(l1 + l2 + l + 1)
        );

    prefactor *= std::sqrt(
        factorial(l1 + m1) * factorial(l1 - m1) *
        factorial(l2 + m2) * factorial(l2 - m2) *
        factorial(l  + m3) * factorial(l  - m3)
    );

    double sum = 0.0;

    int k_min = std::max({
        0,
        l2 - l - m1,
        l1 - l + m2
    });

    int k_max = std::min({
        l1 + l2 - l,
        l1 - m1,
        l2 + m2
    });

    for (int k = k_min; k <= k_max; ++k) {
        double denom =
            factorial(k) *
            factorial(l1 + l2 - l - k) *
            factorial(l1 - m1 - k) *
            factorial(l2 + m2 - k) *
            factorial(l - l2 + m1 + k) *
            factorial(l - l1 - m2 + k);

        sum += std::pow(-1.0, k) / denom;
    }

    return prefactor * sum;
}

// ---------- Clebsch–Gordan ----------
double clebsch_gordan(
    int l1, int m1,
    int l2, int m2,
    int l,  int m) {

    if (m != m1 + m2)
        return 0.0;

    return std::pow(-1.0, l1 - l2 + m)
         * std::sqrt(2.0*l + 1.0)
         * wigner_3j(l1, l2, l, m1, m2, -m);
}

} // anonymous namespace

// ---------- ClebschGordan implementation ----------

ClebschGordan::ClebschGordan(int l_max)
    : l_max_(l_max) {

    assert(l_max_ >= 0);

    for (int l1 = 0; l1 <= l_max_; ++l1) {
        for (int l2 = 0; l2 <= l_max_; ++l2) {
            for (int l = std::abs(l1 - l2); l <= l1 + l2; ++l) {

                ClebschGordanTable table(l1, l2, l);

                for (int m1 = -l1; m1 <= l1; ++m1) {
                    for (int m2 = -l2; m2 <= l2; ++m2) {
                        int m = m1 + m2;
                        if (std::abs(m) > l)
                            continue;

                        double c = clebsch_gordan(
                            l1, m1,
                            l2, m2,
                            l,  m
                        );

                        if (std::abs(c) > 1e-12) {
                            table.entries.push_back({
                                m1, m2, m,
                                static_cast<float>(c)
                            });
                        }
                    }
                }

                tables_.emplace(
                    std::make_tuple(l1, l2, l),
                    std::move(table)
                );
            }
        }
    }
}

const ClebschGordanTable&
ClebschGordan::table(int l1, int l2, int l) const {
    assert(l1 <= l_max_ && l2 <= l_max_);
    assert(std::abs(l1 - l2) <= l && l <= l1 + l2);

    auto it = tables_.find({l1, l2, l});
    assert(it != tables_.end());
    return it->second;
}

} // namespace math

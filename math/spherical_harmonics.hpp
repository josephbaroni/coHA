#pragma once
#include <vector>
#include <cassert>

namespace math {

    // evaluate all real spherical harmonics up to l_max
    // Input: unit vector (x, y, z)
    // Output: concatenated blocks [l=0][l=1]...[l=l_max]
class SphericalHarmonics {
public:
    explicit SphericalHarmonics(int l_max);

    void evaluate(float x, float y, float z, std::vector<float>& out) const;

    [[nodiscard]] int output_dim() const;

private:
    int l_max_;

    static float legendre(int l, int m, float x);

};

}
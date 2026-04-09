#pragma once
#include <array>

namespace math {

// Simple SO(3) rotation via matrix
struct Rotation {
    std::array<float, 9> R;

    static Rotation identity();
    static Rotation from_euler_zyz(float alpha, float beta, float gamma);

    std::array<float, 3> apply(const std::array<float, 3>& v) const;

};

}
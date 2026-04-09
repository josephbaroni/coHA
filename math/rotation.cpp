#include "rotation.hpp"
#include <cmath>

namespace math {

Rotation Rotation::identity() {
    return Rotation{{{
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f,
        0.f, 0.f, 1.f
    }}};
}

// R = Rz(alpha) * Ry(beta) * Rz(gamma)
Rotation Rotation::from_euler_zyz(float alpha, float beta, float gamma) {
    float ca = std::cos(alpha);
    float sa = std::sin(alpha);
    float cb = std::cos(beta);
    float sb = std::sin(beta);
    float cg = std::cos(gamma);
    float sg = std::sin(gamma);

    Rotation rot;
    rot.R = {{
        ca*cb*cg - sa*sg,   -ca*cb*sg - sa*cg,   ca*sb,
        sa*cb*cg + ca*sg,   -sa*cb*sg + ca*cg,   sa*sb,
        -sb*cg,              sb*sg,              cb
    }};
    return rot;
}

std::array<float,3> Rotation::apply(const std::array<float,3>& v) const {
    return {
        R[0]*v[0] + R[1]*v[1] + R[2]*v[2],
        R[3]*v[0] + R[4]*v[1] + R[5]*v[2],
        R[6]*v[0] + R[7]*v[1] + R[8]*v[2]
    };
}

}

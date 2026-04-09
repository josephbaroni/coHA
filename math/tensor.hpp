#pragma once
#include <vector>
#include <cassert>

namespace math {

struct Tensor {
    
    std::vector<float> data;
    std::vector<int> shape;

    Tensor() = default;

    Tensor(std::vector<int> shape_) : shape(std::move(shape_)) {
        int size = 1;
        for (int s : shape) {
            assert(s > 0);
            size *= s;
        }
        data.resize(size, 0.0f);
    }

    float *ptr() {
        return data.data();
    }

    const float *ptr() const {
        return data.data();
    }

    [[nodiscard]] int size() const {
        return static_cast<int>(data.size());
    }

};

}
// memory layout:
// [l block][channel][m]
// GPU friendly for later
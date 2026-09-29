#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "types.hpp"
#include "tolerances.hpp"

namespace sov {

// HVector: Sparse-aware dense work vector for simplex FTRAN, BTRAN, and pricing
class HVector {
public:
    Idx size = 0;
    std::vector<double> val;
    std::vector<Idx> idx;
    Idx count = 0;
    bool packed_valid = true;
    double ema_density = 0.05;

    HVector() = default;
    explicit HVector(Idx n) { resize(n); }

    void resize(Idx n) {
        size = n;
        val.assign(n, 0.0);
        idx.resize(n);
        count = 0;
        packed_valid = true;
    }

    void clear() {
        std::fill(val.begin(), val.end(), 0.0);
        count = 0;
        packed_valid = true;
    }

    void add(Idx i, double v) {
        if (v == 0.0) return;
        if (packed_valid) {
            if (val[i] == 0.0) {
                idx[count++] = i;
                val[i] = v;
            } else {
                val[i] += v;
                if (val[i] == 0.0) {
                    val[i] = kTinySentinel; // Keep membership marker until tighten()
                }
            }
        } else {
            val[i] += v;
        }
    }

    void set(Idx i, double v) {
        if (v == 0.0) {
            if (val[i] != 0.0) {
                val[i] = kTinySentinel;
            }
            return;
        }
        if (packed_valid && val[i] == 0.0) {
            idx[count++] = i;
        }
        val[i] = v;
    }

    void tighten(double drop = kDropTol) {
        if (packed_valid) {
            Idx new_count = 0;
            for (Idx k = 0; k < count; ++k) {
                Idx i = idx[k];
                if (std::abs(val[i]) > drop) {
                    idx[new_count++] = i;
                } else {
                    val[i] = 0.0;
                }
            }
            count = new_count;
        } else {
            count = 0;
            for (Idx i = 0; i < size; ++i) {
                if (std::abs(val[i]) > drop) {
                    idx[count++] = i;
                } else {
                    val[i] = 0.0;
                }
            }
            packed_valid = true;
        }
        if (size > 0) {
            double d = static_cast<double>(count) / static_cast<double>(size);
            ema_density = 0.9 * ema_density + 0.1 * d;
        }
    }

    double density() const {
        return size > 0 ? static_cast<double>(count) / static_cast<double>(size) : 0.0;
    }

    double norm2_sq() const {
        double s = 0.0;
        if (packed_valid) {
            for (Idx k = 0; k < count; ++k) {
                double v = val[idx[k]];
                s += v * v;
            }
        } else {
            for (Idx i = 0; i < size; ++i) s += val[i] * val[i];
        }
        return s;
    }

    void copy_from(const HVector& other) {
        if (size != other.size) resize(other.size);
        clear();
        if (other.packed_valid) {
            count = other.count;
            for (Idx k = 0; k < count; ++k) {
                Idx i = other.idx[k];
                idx[k] = i;
                val[i] = other.val[i];
            }
            packed_valid = true;
        } else {
            val = other.val;
            packed_valid = false;
            tighten();
        }
    }
};

} // namespace sov

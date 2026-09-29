#pragma once
#include <cstdint>
#include <algorithm>
#include <vector>

namespace sov {

inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

inline uint64_t hash_combine(uint64_t h, uint64_t v) {
    uint64_t x = h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
    return splitmix64(x);
}

// Deterministic xoshiro256** random number generator
class Rng {
public:
    explicit Rng(uint64_t seed = 1) {
        reseed(seed);
    }

    void reseed(uint64_t seed) {
        uint64_t x = seed;
        for (int i = 0; i < 4; ++i) {
            s_[i] = splitmix64(x);
        }
    }

    uint64_t next_u64() {
        const uint64_t result = rotl(s_[1] * 5, 7) * 9;
        const uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    // Uniform in [0, 1)
    double uniform01() {
        return (next_u64() >> 11) * 0x1.0p-53;
    }

    // Rejection-sampled uniform integer in [lo, hi] inclusive
    int64_t uniform_int(int64_t lo, int64_t hi) {
        if (lo >= hi) return lo;
        uint64_t range = static_cast<uint64_t>(hi - lo) + 1ULL;
        uint64_t threshold = (0ULL - range) % range;
        for (;;) {
            uint64_t r = next_u64();
            if (r >= threshold) {
                return lo + static_cast<int64_t>(r % range);
            }
        }
    }

    template <typename T>
    void shuffle(std::vector<T>& v) {
        for (size_t i = v.size(); i > 1; --i) {
            size_t j = static_cast<size_t>(uniform_int(0, static_cast<int64_t>(i - 1)));
            std::swap(v[i - 1], v[j]);
        }
    }

private:
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
    uint64_t s_[4];
};

} // namespace sov

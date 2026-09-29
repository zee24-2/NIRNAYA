#pragma once
#include <cmath>
#include <vector>

namespace sov {

// Knuth two-sum: s + err = a + b exactly
inline void two_sum(double a, double b, double& s, double& err) {
    s = a + b;
    double bb = s - a;
    err = (a - (s - bb)) + (b - bb);
}

// Fast two-sum when |a| >= |b|
inline void quick_two_sum(double a, double b, double& s, double& err) {
    s = a + b;
    err = b - (s - a);
}

// Two-product using hardware/standard fused multiply-add: p + err = a * b exactly
inline void two_prod(double a, double b, double& p, double& err) {
    p = a * b;
    err = std::fma(a, b, -p);
}

// Double-double (~32 decimal digits of precision) for residual & checker calculations
struct dd_real {
    double hi = 0.0;
    double lo = 0.0;

    constexpr dd_real() = default;
    constexpr dd_real(double h, double l = 0.0) : hi(h), lo(l) {}

    double to_double() const { return hi + lo; }

    dd_real operator-() const { return dd_real(-hi, -lo); }

    dd_real& operator+=(double b) {
        double s1, s2;
        two_sum(hi, b, s1, s2);
        s2 += lo;
        quick_two_sum(s1, s2, hi, lo);
        return *this;
    }

    dd_real& operator+=(const dd_real& b) {
        double s1, s2, t1, t2;
        two_sum(hi, b.hi, s1, s2);
        two_sum(lo, b.lo, t1, t2);
        s2 += t1;
        quick_two_sum(s1, s2, s1, s2);
        s2 += t2;
        quick_two_sum(s1, s2, hi, lo);
        return *this;
    }

    dd_real& operator-=(const dd_real& b) {
        return (*this) += (-b);
    }

    dd_real& operator*=(double b) {
        double p1, p2;
        two_prod(hi, b, p1, p2);
        p2 += lo * b;
        quick_two_sum(p1, p2, hi, lo);
        return *this;
    }

    dd_real& operator*=(const dd_real& b) {
        double p1, p2;
        two_prod(hi, b.hi, p1, p2);
        p2 += (hi * b.lo + lo * b.hi);
        quick_two_sum(p1, p2, hi, lo);
        return *this;
    }

    void add_prod(double a, double b) {
        double p, e;
        two_prod(a, b, p, e);
        (*this) += dd_real(p, e);
    }
};

inline dd_real operator+(dd_real a, const dd_real& b) { a += b; return a; }
inline dd_real operator-(dd_real a, const dd_real& b) { a -= b; return a; }
inline dd_real operator*(dd_real a, const dd_real& b) { a *= b; return a; }

// Neumaier compensated sum for fast accurate dot products
struct NeumaierSum {
    double s = 0.0;
    double c = 0.0;

    void add(double x) {
        double t = s + x;
        if (std::abs(s) >= std::abs(x)) {
            c += (s - t) + x;
        } else {
            c += (x - t) + s;
        }
        s = t;
    }

    double result() const { return s + c; }
};

} // namespace sov

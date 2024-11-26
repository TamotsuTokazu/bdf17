#ifndef POLY_H
#define POLY_H

#include <cstdint>
#include <stdexcept>
#include <type_traits>

#include <iostream>

#include "ntt.h"

template <typename NTT>
class Poly {
public:

    static const NTT ntt;

    static constexpr uint64_t p = NTT::p;
    static constexpr size_t O = NTT::O;
    static constexpr size_t N = NTT::N;

    using ZZ = Zp<p>;

    bool is_coeff = false;
    uint64_t a[N];

    Poly(bool is_coeff_ = true) : is_coeff(is_coeff_), a{} {}

    void ToCoeff() {
        if (is_coeff) {
            return;
        }
        NTT::GetInstance().InverseNTT(a);
        is_coeff = true;
    }

    void ToNTT() {
        if (!is_coeff) {
            return;
        }
        NTT::GetInstance().ForwardNTT(a);
        is_coeff = false;
    }

    Poly operator+(const Poly &rhs) const {
        if (is_coeff != rhs.is_coeff) {
            throw std::runtime_error("Addition is not supported between coefficient and NTT domain");
        }

        Poly ret(is_coeff);
        ret.is_coeff = is_coeff;
        for (size_t i = 0; i < N; i++) {
            ret.a[i] = ZZ::Add(a[i], rhs.a[i]);
        }
        return ret;
    }

    Poly operator-(const Poly &rhs) const {
        if (is_coeff != rhs.is_coeff) {
            throw std::runtime_error("Addition is not supported between coefficient and NTT domain");
        }

        Poly ret(is_coeff);
        for (size_t i = 0; i < N; i++) {
            ret.a[i] = ZZ::Sub(a[i], rhs.a[i]);
        }
        return ret;
    }

    Poly operator*(const Poly &rhs) const {
        if (is_coeff || rhs.is_coeff) {
            throw std::runtime_error("Multiplication is not supported in coefficient domain");
        }

        Poly ret(is_coeff);
        for (size_t i = 0; i < N; i++) {
            ret.a[i] = ZZ::Mul(a[i], rhs.a[i]);
        }
        return ret;
    }

    Poly operator*(uint64_t rhs) const {
        Poly ret(is_coeff);
        for (size_t i = 0; i < N; i++) {
            ret.a[i] = ZZ::Mul(a[i], rhs);
        }
        return ret;
    }

    void print() {
        if (is_coeff) {
            std::cout << "Coeff: ";
        } else {
            std::cout << "NTT: ";
        }
        for (size_t i = 0; i < N; i++) {
            std::cout << a[i] << " ";
        }
        std::cout << std::endl;
    }
};

#endif // POLY_H
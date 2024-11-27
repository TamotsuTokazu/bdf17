#ifndef RLWE_H
#define RLWE_H

#include <vector>
#include <array>
#include <random>

#include "poly.h"

template <typename Poly, uint64_t B>
class SchemeImpl {
public:
    using RLWEKey = std::vector<Poly>;
    using RLWECiphertext = std::vector<Poly>;
    using RLWEGadgetCiphertext = std::vector<RLWECiphertext>;
    using RLWESwitchingKey = std::vector<RLWEGadgetCiphertext>;
    using RGSWCiphertext = std::pair<RLWEGadgetCiphertext, RLWEGadgetCiphertext>;

    constexpr static uint64_t Q = Poly::p;

    constexpr static size_t G = []() {
        size_t G = 0;
        for (__uint128_t t = 1; t <= Q; t *= B) {
            G++;
        }
        return G;
    }();

    constexpr static auto gadget = []() {
        std::array<uint64_t, G> gadget;
        uint64_t t = 1;
        for (size_t i = 0; i < G; i++) {
            gadget[i] = t;
            t *= B;
        }
        return gadget;
    }();

    RLWEKey sk;
    Poly skp;

    std::vector<RLWESwitchingKey> ksk_galois;

    std::mt19937_64 engine;
    std::uniform_int_distribution<uint64_t> distribution;

    SchemeImpl(bool keygen = true);
    void GaloisKeyGen();

    Poly GaloisConjugate(const Poly &x, const size_t &a);

    template <typename T>
    std::vector<T> GaloisConjugate(const std::vector<T> &x, const size_t &a);

    template <typename T1, typename T2>
    std::pair<T1, T2> GaloisConjugate(const std::pair<T1, T2> &x, const size_t &a);

    void ModSwitch(Poly &x, uint64_t q);
    void ModSwitch(RLWECiphertext &ct, uint64_t q);

    RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, uint64_t q_plain);
    RLWEGadgetCiphertext RLWEGadgetEncrypt(const Poly &m, const RLWEKey &sk, uint64_t q_plain);
    RGSWCiphertext RGSWEncrypt(const Poly &m, const RLWEKey &sk);
    Poly RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, uint64_t q_plain);

    RLWECiphertext Mult(Poly a, RLWEGadgetCiphertext ct);
    RLWECiphertext ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW);

    RLWESwitchingKey KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN);
    RLWECiphertext KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K);

    std::vector<RGSWCiphertext> BootstrappingKeyGen(std::vector<uint64_t> z);
    RLWECiphertext Process(const std::vector<RGSWCiphertext> &BK, std::vector<uint64_t> a, uint64_t b, uint64_t q_plain);
};

#include "rlwe-impl.h"

#endif // RLWE_H
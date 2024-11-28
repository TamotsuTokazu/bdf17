#ifndef RLWE_H
#define RLWE_H

#include <vector>
#include <array>
#include <random>

#include "poly.h"

template <size_t l>
class GaussianSampler {
public:
    static GaussianSampler &GetInstance() {
        static GaussianSampler instance;
        return instance;
    }

    std::vector<int64_t> SampleE(double var) {
        std::vector<int64_t> a(l, 0);
        for (size_t i = 0; i < var * l / 2; i++) {
            a[distribution(engine)]++;
            a[distribution(engine)]--;
        }
        return a;
    }

    std::vector<int64_t> SampleSk(double density) {
        std::vector<int64_t> a(l, 0);
        size_t count1 = density * l;
        while (count1 > 0) {
            size_t i = distribution(engine);
            if (a[i] == 0) {
                a[i] = 1;
                count1--;
            }
        }
        count1 = density * l;
        while (count1 > 0) {
            size_t i = distribution(engine);
            if (a[i] == 0) {
                a[i] = -1;
                count1--;
            }
        }
        return a;
    }
private:
    std::mt19937 engine;
    std::uniform_int_distribution<int64_t> distribution;

    GaussianSampler() : engine(std::random_device{}()), distribution(0, l - 1) {}
};

template <typename Poly_, uint64_t B>
class SchemeImpl {
public:
    using Poly = Poly_;

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

    SchemeImpl();

    SchemeImpl(std::vector<int64_t> skVec);

    void GaloisKeyGen();

    static Poly GaloisConjugate(const Poly &x, const size_t &a);

    template <typename T>
    static std::vector<T> GaloisConjugate(const std::vector<T> &x, const size_t &a);

    template <typename T1, typename T2>
    static std::pair<T1, T2> GaloisConjugate(const std::pair<T1, T2> &x, const size_t &a);

    static void ModSwitch(Poly &x, uint64_t q);

    template <typename S>
    static typename S::RLWECiphertext ModSwitch(const RLWECiphertext &ct);

    RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, uint64_t q_plain);
    RLWEGadgetCiphertext RLWEGadgetEncrypt(const Poly &m, const RLWEKey &sk, uint64_t q_plain);
    RGSWCiphertext RGSWEncrypt(const Poly &m, const RLWEKey &sk);
    static Poly RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, uint64_t q_plain);

    static RLWECiphertext Mult(Poly a, RLWEGadgetCiphertext ct);
    static RLWECiphertext ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW);

    RLWESwitchingKey KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN);
    static RLWECiphertext KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K);

    std::vector<RGSWCiphertext> BootstrappingKeyGen(std::vector<int64_t> z);
    RLWECiphertext Process(const std::vector<RGSWCiphertext> &BK, std::vector<int64_t> a, int64_t b, uint64_t q_plain);
};

#include "rlwe-impl.h"

#endif // RLWE_H
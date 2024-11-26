#ifndef RLWE_IMPL_H
#define RLWE_IMPL_H

template <typename Poly, uint64_t B>
SchemeImpl<Poly, B>::SchemeImpl(bool keygen) : sk(), skp() {
    if (keygen) {
        skp.a[0] = 1;
        skp.ToNTT();
        sk.push_back(skp);
    }
}

template <typename Poly, uint64_t B>
void SchemeImpl<Poly, B>::GaloisKeyGen() {
    ksk_galois.resize(Poly::O);
    skp.ToNTT();
    for (size_t a = 2; a < Poly::O; a++) {
        auto sk_a = GaloisConjugate(sk, a);
        auto ksk = KeySwitchGen(sk_a, sk);
        ksk_galois[a] = ksk;
    }
}

template <typename Poly, uint64_t B>
Poly SchemeImpl<Poly, B>::GaloisConjugate(const Poly &x, const size_t &a) {
    Poly ret(x.is_coeff);
    if (x.is_coeff) {
        for (size_t i = 1; i <= Poly::N; i++) {
            ret.a[i * a % Poly::O - 1] = x.a[i - 1];
        }
    } else {
        for (size_t i = 1; i <= Poly::N; i++) {
            ret.a[i - 1] = x.a[i * a % Poly::O - 1];
        }
    }
    return ret;
}

template <typename Poly, uint64_t B>
template <typename T>
std::vector<T> SchemeImpl<Poly, B>::GaloisConjugate(const std::vector<T> &x, const size_t &a) {
    std::vector<T> ret;
    ret.reserve(x.size());
    for (const auto &elem : x) {
        ret.push_back(GaloisConjugate(elem, a));
    }
    return ret;
}

template <typename Poly, uint64_t B>
template <typename T1, typename T2>
std::pair<T1, T2> SchemeImpl<Poly, B>::GaloisConjugate(const std::pair<T1, T2> &x, const size_t &a) {
    return {GaloisConjugate(x.first, a), GaloisConjugate(x.second, a)};
}

template <typename Poly, uint64_t B>
typename SchemeImpl<Poly, B>::RLWECiphertext SchemeImpl<Poly, B>::RLWEEncrypt(const Poly &m, const RLWEKey &sk, uint64_t q_plain) {
    size_t k = sk.size();
    RLWECiphertext ct;

    Poly result(false);
    for (size_t i = 0; i < k; i++) {
        Poly a(false);
        result = result + a * sk[i];
        ct.push_back(a);
    }
    ct.push_back(result + m * (Q / q_plain));
    return ct;
}

template <typename Poly, uint64_t B>
typename SchemeImpl<Poly, B>::RLWEGadgetCiphertext SchemeImpl<Poly, B>::RLWEGadgetEncrypt(const Poly &m, const RLWEKey &sk, uint64_t q_plain) {
    RLWEGadgetCiphertext ct(G);

    for (size_t i = 0; i < G; i++) {
        ct[i] = RLWEEncrypt(m * gadget[i], sk, q_plain);
    }
    return ct;
}

template <typename Poly, uint64_t B>
typename SchemeImpl<Poly, B>::RGSWCiphertext SchemeImpl<Poly, B>::RGSWEncrypt(const Poly &m, const typename SchemeImpl<Poly, B>::RLWEKey &sk) {
    if (sk.size() != 1) {
        throw std::runtime_error("RGSW encryption requires a secret key of size 1");
    }
    Poly s = sk[0];
    return std::make_pair(RLWEGadgetEncrypt(m * s, sk, Q), RLWEGadgetEncrypt(m, sk, Q));
}

template <typename Poly, uint64_t B>
Poly SchemeImpl<Poly, B>::RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, uint64_t q_plain) {
    size_t k = sk.size();
    Poly result = ct[k];
    for (size_t i = 0; i < k; i++) {
        result = result - ct[i] * sk[i];
    }
    result.ToCoeff();
    ModSwitch(result, q_plain);
    return result;
}

template <typename Poly, uint64_t B>
void SchemeImpl<Poly, B>::ModSwitch(Poly &x, uint64_t q) {
    size_t n = Poly::N;
    __uint128_t Q = Poly::p;
    __uint128_t halfQ = Q / 2;
    __uint128_t qq = q;
    for (size_t i = 0; i < n; ++i) {
        __uint128_t xi = x.a[i];
        xi = (xi * qq + halfQ) / Q;
        if (xi >= qq) {
            xi -= qq;
        }
        x.a[i] = (uint64_t)xi;
    }
}

template <typename Poly, uint64_t B>
void SchemeImpl<Poly, B>::ModSwitch(RLWECiphertext &ct, uint64_t q) {
    for (auto &elem : ct) {
        ModSwitch(elem, q);
    }
}

template <typename Poly, uint64_t B>
typename SchemeImpl<Poly, B>::RLWESwitchingKey SchemeImpl<Poly, B>::KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN) {
    RLWESwitchingKey result;
    for (size_t i = 0; i < sk.size(); i++) {
        result.push_back(RLWEGadgetEncrypt(sk[i], skN, Q));
    }
    return result;
}

template <typename Poly, uint64_t B>
typename SchemeImpl<Poly, B>::RLWECiphertext SchemeImpl<Poly, B>::KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K) {
    size_t k = ct.size() - 1;
    size_t kN = K[0][0].size() - 1;

    RLWECiphertext result(kN + 1);
    for (size_t i = 0; i <= kN; i++) {
        result[i].is_coeff = false;
    }
    result[kN] = ct[k];

    for (size_t i = 0; i < k; i++) {
        auto a = ct[i];
        a.toCoeff();
        for (size_t j = 0; j < G; j++) {
            uint64_t t = gadget[j];
            Poly a0(true);
            for (size_t l = 0; l < Poly::N; l++) {
                a0.a[l] = a.a[l] / t % B;
            }
            a0.toNTT();
            for (size_t l = 0; l <= kN; l++) {
                result[l] = result[l] - a0 * K[i][j][l];
            }
        }
    }
    return result;
}

#endif // RLWE_IMPL_H
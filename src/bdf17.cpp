#include <iostream>
#include <random>
#include <chrono>

#include "ntt.h"
#include "rlwe.h"

#define START_TIMER start = std::chrono::system_clock::now()
#define END_TIMER std::cout << "Time: " << (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()) << std::endl
auto start = std::chrono::system_clock::now();

using NTTp = CircNTT<72057421557668737LL, 5LL, 1153, 5>;
using NTTq = CircNTT<72057421557668737LL, 5LL, 1297, 10>;

using NTTpt = CircNTT<108533126017LL, 10LL, 1153, 5>;
using NTTqt = CircNTT<108533126017LL, 10LL, 1297, 10>;

using NTTpq = TensorNTTImpl<NTTpt, NTTqt>;

using PolyP = Poly<NTTp>;
using PolyQ = Poly<NTTq>;
using PolyPt = Poly<NTTpt>;
using PolyQt = Poly<NTTqt>;

using PolyPQ = Poly<NTTpq>;

using Z = NTTpq::Z;

constexpr size_t n = 600;
constexpr uint64_t Qplain = 64;
constexpr uint64_t Bks = 1 << 8;
constexpr size_t pq = PolyP::N * PolyQ::N;
constexpr size_t N_tests = 8;

using SchemeP = SchemeImpl<PolyP, Bks>;
using SchemeQ = SchemeImpl<PolyQ, Bks>;
using SchemePQ = SchemeImpl<PolyPQ, Bks>;
using SchemePt = SchemeImpl<PolyPt, Bks>;
using SchemeQt = SchemeImpl<PolyQt, Bks>;

std::mt19937 engine(std::random_device{}());
std::uniform_int_distribution<size_t> distributionp(0, pq - 1);

template <typename PolyPT, typename PolyQT>
Poly<TensorNTTImpl<typename PolyPT::NTT, typename PolyQT::NTT>> Tensor(const PolyPT &a, const PolyQT &b) {
    if (a.is_coeff || b.is_coeff) {
        throw std::runtime_error("Tensor product is not supported for coefficient domain");
    }

    Poly<TensorNTTImpl<typename PolyPT::NTT, typename PolyQT::NTT>> c(false);
    for (size_t i = 0; i < PolyPT::N; i++) {
        for (size_t j = 0; j < PolyQT::N; j++) {
            c.a[i * PolyQT::N + j] = Z::Mul(a.a[i], b.a[j]);
        }
    }
    return c;
}

template <typename PolyPT, typename PolyQT>
Poly<TensorNTTImpl<typename PolyPT::NTT, typename PolyQT::NTT>> GenKey(std::vector<int64_t> &sk) {
    PolyPT skp(true);
    PolyQT one(false);
    for (size_t i = 0; i < sk.size(); i++) {
        skp.a[i] = sk[i] < 0 ? sk[i] + PolyPT::p : sk[i];
    }
    for (size_t i = 0; i < PolyQT::N; i++) {
        one.a[i] = 1;
    }
    skp.ToNTT();
    return Tensor(skp, one);
}

typename SchemePQ::RLWEKey TensorKey(const SchemePt::RLWEKey &skp, const SchemeQt::RLWEKey &skq) {
    typename SchemePQ::RLWEKey skpq;
    PolyPt skp0 = skp[0];
    PolyQt skq0 = skq[0];
    PolyPt p1(false);
    PolyQt q1(false);
    for (size_t i = 0; i < PolyPt::N; i++) {
        p1.a[i] = 1;
    }
    for (size_t i = 0; i < PolyQt::N; i++) {
        q1.a[i] = 1;
    }
    PolyPQ skpq0 = Tensor(skp0, skq0);
    for (size_t i = 0; i < PolyPQ::N; i++) {
        skpq0.a[i] = Z::Sub(0, skpq0.a[i]);
    }
    skpq.push_back(skpq0);
    skpq.push_back(Tensor(skp0, q1));
    skpq.push_back(Tensor(p1, skq0));
    return skpq;
}

typename SchemePQ::RLWECiphertext TensorCt(const SchemePt::RLWECiphertext &ctp, const SchemeQt::RLWECiphertext &ctq) {
    typename SchemePQ::RLWECiphertext ct;
    uint64_t z = SchemePQ::Q - Qplain;
    ct.push_back(Tensor(ctp[0], ctq[0]) * z);
    ct.push_back(Tensor(ctp[0], ctq[1]) * z);
    ct.push_back(Tensor(ctp[1], ctq[0]) * z);
    ct.push_back(Tensor(ctp[1], ctq[1]) * z);
    return ct;
}

PolyPt TracePQtoP(const PolyPQ &a) {
    if (a.is_coeff) {
        PolyPt b(true);
        for (size_t i = 0; i < PolyPt::N; i++) {
            b.a[i] = a.a[i * PolyQt::N];
        }
        return b;
    } else {
        uint64_t z = Z::Pow(PolyQt::N, Z::p - 2);
        PolyPt b(false);
        for (size_t i = 0; i < PolyPt::N; i++) {
            size_t ii = i;
            for (size_t j = 0; j < PolyQt::N; j++) {
                b.a[ii] = Z::Add(b.a[ii], a.a[i * PolyQt::N + j]);
            }
            b.a[ii] = Z::Mul(b.a[ii], z);
        }
        return b;
    }
}

template <typename Poly>
uint64_t TracePtoZ(const Poly &a) {
    if (a.is_coeff) {
        return a.a[0];
    } else {
        uint64_t z = 0;
        for (size_t i = 0; i < Poly::N; i++) {
            z = Poly::Z::Add(z, a.a[i]);
        }
        return Poly::Z::Mul(z, Z::Pow(Poly::N, Z::p - 2));
    }
}

PolyPQ ConstructF(std::vector<size_t> &f) {
    PolyPQ a(true);
    for (size_t k = 0; k < f.size(); k++) {
        size_t i = k % PolyP::O;
        size_t j = k % PolyQ::O;
        a.a[i * PolyQ::N + j] = f[k];
    }
    return a;
}

std::uniform_int_distribution<size_t> distribution(0, Qplain - 1);

int main() {

    START_TIMER;

    std::vector<int64_t> sk = GaussianSampler<n>::GetInstance().SampleSk(0.33);

    std::vector<size_t> f_plain(Qplain, 0);

    // define f
    for (size_t i = 0; i < Qplain; i++) {
        f_plain[i] = i & 1;
    }

    std::vector<size_t> f_ct(pq, 0);
    for (size_t i = 0; i < pq; i++) {
        f_ct[i] = f_plain[(size_t)(0.5 + (double)Qplain * i / pq) % Qplain];
    }

    for (size_t i = 1, j = pq - 1; i < j; i++, j--) {
        std::swap(f_ct[i], f_ct[j]);
    }

    auto skp = GaussianSampler<PolyP::N>::GetInstance().SampleSk(0.3);
    auto skq = GaussianSampler<PolyQ::N>::GetInstance().SampleSk(0.3);

    SchemeP schemeP(skp);
    schemeP.GaloisKeyGen();
    auto BKp = schemeP.BootstrappingKeyGen(sk);
    SchemeQ schemeQ(skq);
    schemeQ.GaloisKeyGen();
    auto BKq = schemeQ.BootstrappingKeyGen(sk);

    SchemePt schemePt(skp);
    SchemeQt schemeQt(skq);

    SchemePQ schemePQ;
    auto skpq = TensorKey(schemePt.sk, schemeQt.sk);
    SchemePQ::RLWEKey skp0{GenKey<PolyPt, PolyQt>(sk)};
    auto tensorBK = schemePQ.KeySwitchGen(skpq, skp0);

    END_TIMER;

    for (size_t n_test = 0; n_test < N_tests; n_test++) {

        std::vector<int64_t> a(n);
        for (size_t i = 0; i < n; i++) {
            a[i] = distributionp(engine);
        }

        int64_t b0 = distribution(engine);

        int64_t b = b0 * pq / Qplain;
        for (size_t i = 0; i < n; i++) {
            b = (b + sk[i] * a[i]) % (int64_t)pq;
        }
        b = (b + pq) % pq;

        START_TIMER;

        auto ctp = SchemeP::ModSwitch<SchemePt>(schemeP.Process(BKp, a, b, Qplain));
        auto ctq = SchemeQ::ModSwitch<SchemeQt>(schemeQ.Process(BKq, a, b, Qplain));

        END_TIMER;

        START_TIMER;

        auto ctpq = TensorCt(ctp, ctq);

        auto tensor_ct = schemePQ.KeySwitch(ctpq, tensorBK);

        auto f = ConstructF(f_ct);
        f.ToNTT();

        tensor_ct[0] = f * tensor_ct[0];
        tensor_ct[1] = f * tensor_ct[1];

        SchemePt::RLWECiphertext ct_trace = {TracePQtoP(tensor_ct[0]), TracePQtoP(tensor_ct[1])};
        SchemePt::RLWEKey sk_trace = {TracePQtoP(skp0[0])};

        auto b_out = TracePtoZ(ct_trace[1]);

        std::vector<uint64_t> a_out(n);
        ct_trace[0].ToCoeff();
        a_out[0] = ct_trace[0].a[0];
        for (size_t i = 1; i < n; i++) {
            a_out[i] = ct_trace[0].a[PolyP::N - i];
        }

        END_TIMER;

        for (size_t i = 0; i < n; i++) {
            b_out = Z::Sub(b_out, Z::Mul(sk[i] + Z::p, a_out[i]));
        }
        std::cout << "result: " << (size_t)(0.5 + (double)Qplain * b_out / Z::p) % Qplain << std::endl;
        std::cout << "expected: " << f_plain[b0] << std::endl;

        if (f_plain[b0] != (size_t)(0.5 + (double)Qplain * b_out / Z::p) % Qplain) {
            std::cout << "Error" << std::endl;
            return 1;
        }
    }

    return 0;
}
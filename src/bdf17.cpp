#include <iostream>
#include <random>
#include <chrono>

#include "ntt.h"
#include "rlwe.h"

#define START_TIMER start = std::chrono::system_clock::now()
#define END_TIMER std::cout << "Time: " << (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()) << std::endl
auto start = std::chrono::system_clock::now();

using NTTp = CircNTT<1152920977604149249LL, 7LL, 1153, 5>;
using NTTq = CircNTT<1152920977604149249LL, 7LL, 1297, 10>;

// using NTTp = CircNTT<1152920959144329601LL, 22LL, 433, 5>;
// using NTTq = CircNTT<1152920959144329601LL, 22LL, 487, 3>;

using NTTpq = TensorNTTImpl<NTTp, NTTq>;

using PolyP = Poly<NTTp>;
using PolyQ = Poly<NTTq>;

using PolyPQ = Poly<NTTpq>;

using ZZ = NTTpq::ZZ;

const size_t n = 600;
const uint64_t Qplain = 64;
const uint64_t Bks = 1 << 8;
const size_t pq = PolyP::N * PolyQ::N;

using SchemeP = SchemeImpl<PolyP, Bks>;
using SchemeQ = SchemeImpl<PolyQ, Bks>;
using SchemePQ = SchemeImpl<PolyPQ, Bks>;

std::mt19937 engine(std::random_device{}());
std::uniform_int_distribution<size_t> distribution3(0, 2);
std::uniform_int_distribution<size_t> distributionp(0, pq - 1);

PolyPQ Tensor(const PolyP &a, const PolyQ &b) {
    if (a.is_coeff || b.is_coeff) {
        throw std::runtime_error("Tensor product is not supported for coefficient domain");
    }
    
    PolyPQ c(false);
    for (size_t i = 0; i < PolyP::N; i++) {
        for (size_t j = 0; j < PolyQ::N; j++) {
            c.a[i * PolyQ::N + j] = ZZ::Mul(a.a[i], b.a[j]);
        }
    }
    return c;
}

PolyPQ GenKey(std::vector<size_t> &sk) {
    PolyP skp(true);
    PolyQ one(false);
    for (size_t i = 0; i < sk.size(); i++) {
        skp.a[i] = sk[i];
    }
    for (size_t i = 0; i < PolyQ::N; i++) {
        one.a[i] = 1;
    }
    skp.ToNTT();
    return Tensor(skp, one);
}

typename SchemePQ::RLWEKey TensorKey(const SchemeP::RLWEKey &skp, const SchemeQ::RLWEKey &skq) {
    typename SchemePQ::RLWEKey skpq;
    PolyP skp0 = skp[0];
    PolyQ skq0 = skq[0];
    PolyP p1(false);
    PolyQ q1(false);
    for (size_t i = 0; i < PolyP::N; i++) {
        p1.a[i] = 1;
    }
    for (size_t i = 0; i < PolyQ::N; i++) {
        q1.a[i] = 1;
    }
    PolyPQ skpq0 = Tensor(skp0, skq0);
    for (size_t i = 0; i < PolyPQ::N; i++) {
        skpq0.a[i] = ZZ::Sub(0, skpq0.a[i]);
    }
    skpq.push_back(skpq0);
    skpq.push_back(Tensor(skp0, q1));
    skpq.push_back(Tensor(p1, skq0));
    return skpq;
}

typename SchemePQ::RLWECiphertext TensorCt(const SchemeP::RLWECiphertext &ctp, const SchemeQ::RLWECiphertext &ctq) {
    typename SchemePQ::RLWECiphertext ct;
    uint64_t z = SchemePQ::Q - Qplain;
    ct.push_back(Tensor(ctp[0], ctq[0]) * z);
    ct.push_back(Tensor(ctp[0], ctq[1]) * z);
    ct.push_back(Tensor(ctp[1], ctq[0]) * z);
    ct.push_back(Tensor(ctp[1], ctq[1]) * z);
    return ct;
}

PolyP TracePQtoP(const PolyPQ &a) {
    if (a.is_coeff) {
        PolyP b(true);
        for (size_t i = 0; i < PolyP::N; i++) {
            b.a[i] = a.a[i * PolyQ::N];
        }
        return b;
    } else {
        uint64_t z = ZZ::Pow(PolyQ::N, ZZ::p - 2);
        PolyP b(false);
        for (size_t i = 0; i < PolyP::N; i++) {
            size_t ii = i;
            for (size_t j = 0; j < PolyQ::N; j++) {
                b.a[ii] = ZZ::Add(b.a[ii], a.a[i * PolyQ::N + j]);
            }
            b.a[ii] = ZZ::Mul(b.a[ii], z);
        }
        return b;
    }
}

uint64_t TracePtoZ(const PolyP &a) {
    if (a.is_coeff) {
        return a.a[0];
    } else {
        uint64_t z = 0;
        for (size_t i = 0; i < PolyP::N; i++) {
            z = ZZ::Add(z, a.a[i]);
        }
        return ZZ::Mul(z, ZZ::Pow(PolyP::N, ZZ::p - 2));
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

int main() {

    std::vector<size_t> a(n);
    std::vector<size_t> sk(n);

    for (size_t i = 0; i < n; i++) {
        sk[i] = distribution3(engine);
        if (sk[i] == 2) {
            sk[i] = pq - 1;
        }
        a[i] = distributionp(engine);
    }

    uint64_t b0 = 18;
    uint64_t b = b0 * pq / Qplain;
    for (size_t i = 0; i < n; i++) {
        b = (b + sk[i] * a[i]) % pq;
    }

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

    SchemeP schemeP(true);
    schemeP.GaloisKeyGen();
    auto BKp = schemeP.BootstrappingKeyGen(sk);
    SchemeQ schemeQ(true);
    schemeQ.GaloisKeyGen();
    auto BKq = schemeQ.BootstrappingKeyGen(sk);

    SchemePQ schemePQ(false);
    auto skpq = TensorKey(schemeP.sk, schemeQ.sk);
    SchemePQ::RLWEKey skp0{GenKey(sk)};
    auto tensorBK = schemePQ.KeySwitchGen(skpq, skp0);

    START_TIMER;

    auto ctp = schemeP.Process(BKp, a, b, Qplain);
    auto ctq = schemeQ.Process(BKq, a, b, Qplain);

    END_TIMER;

    START_TIMER;

    auto ctpq = TensorCt(ctp, ctq);

    auto tensor_ct = schemePQ.KeySwitch(ctpq, tensorBK);

    auto f = ConstructF(f_ct);
    f.ToNTT();

    tensor_ct[0] = f * tensor_ct[0];
    tensor_ct[1] = f * tensor_ct[1];

    SchemeP::RLWECiphertext ct_trace = {TracePQtoP(tensor_ct[0]), TracePQtoP(tensor_ct[1])};
    SchemeP::RLWEKey sk_trace = {TracePQtoP(skp0[0])};

    auto b_out = TracePtoZ(ct_trace[1]);

    std::vector<uint64_t> a_out(n);
    ct_trace[0].ToCoeff();
    a_out[0] = ct_trace[0].a[0];
    for (size_t i = 1; i < n; i++) {
        a_out[i] = ct_trace[0].a[PolyP::N - i];
    }

    END_TIMER;

    for (size_t i = 0; i < n; i++) {
        b_out = ZZ::Sub(b_out, ZZ::Mul(sk[i], a_out[i]));
    }
    std::cout << "b: " << b_out << " " << (b_out * Qplain + ZZ::p / 2) / ZZ::p << std::endl;
    std::cout << "expected: " << f_plain[b0] << std::endl;

    return 0;
}
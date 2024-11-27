#include <iostream>
#include <random>
#include <chrono>

#include "ntt.h"
#include "rlwe.h"

using NTTp = NTT<1152920977604149249LL, 7LL, 1153, 5>;
using NTTq = NTT<1152920977604149249LL, 7LL, 1297, 10>;

// using NTTp = NTT<1152921504602622529LL, 14LL, 37, 2>;
// using NTTq = NTT<1152921504602622529LL, 14LL, 17, 3>;

using NTTpq = TensorNTTImpl<NTTp, NTTq>;

using PolyP = Poly<NTTp>;

SchemeImpl<PolyP, 1 << 8> schemeP;

std::mt19937 engine(std::random_device{}());
std::uniform_int_distribution<size_t> distribution3(0, 2);
std::uniform_int_distribution<size_t> distributionp(0, PolyP::O - 1);

int main() {
    size_t n = 600;
    const uint64_t Qplain = 64;

    std::vector<size_t> a(n);
    std::vector<size_t> sk(n);

    for (size_t i = 0; i < n; i++) {
        sk[i] = distribution3(engine);
        a[i] = distributionp(engine);
    }

    uint64_t b = 10;
    for (size_t i = 0; i < n; i++) {
        b = (b + sk[i] * a[i]) % PolyP::O;
    }

    schemeP.GaloisKeyGen();
    auto BK = schemeP.BootstrappingKeyGen(sk);

    auto start = std::chrono::high_resolution_clock::now();

    auto ct = schemeP.Process(BK, a, b, Qplain);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    auto m = schemeP.RLWEDecrypt(ct, schemeP.sk, Qplain);
    m.print();
    return 0;
}
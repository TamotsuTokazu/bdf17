#include <iostream>
#include "ntt.h"
#include "rlwe.h"

using NTTp = NTT<1152921504606588289LL, 7, 13, 2>;
using NTTq = NTT<1152921504606588289LL, 7, 17, 3>;

using NTTpq = TensorNTTImpl<NTTp, NTTq>;

using PolyP = Poly<NTTp>;

SchemeImpl<PolyP, 2> schemeP;

int main() {
    PolyP m;
    m.a[0] = 1;
    m.print();
    m.ToNTT();
    schemeP.GaloisKeyGen();
    auto ct = schemeP.RLWEEncrypt(m, schemeP.sk, 64);
    auto m2 = schemeP.RLWEDecrypt(ct, schemeP.sk, 64);
    m2.print();
    return 0;
}
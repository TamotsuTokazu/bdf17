def gen_prime(N, m, q=64):
    m = lcm(m, q)
    N = N // m * m + 1;
    while not is_prime(N):
        N -= m
    g = primitive_root(N)
    return N, g

p, q = 433, 487
m = p * q * (p - 1) * (q - 1)
N, g = gen_prime(2^60, m)
print(factor(N - 1))

print(f'{N}LL, {g}LL, {p}, {primitive_root(p)}')
print(f'{N}LL, {g}LL, {q}, {primitive_root(q)}')
def gen_prime(N, m):
    N = N // m * m + 1;
    while not is_prime(N):
        N -= m
    g = primitive_root(N)
    return N, g

p, q = 37, 17
m = p * q * (p - 1) * (q - 1)
N, g = gen_prime(2^60, m)

print(f'{N}LL, {g}LL, {p}, {primitive_root(p)}')
print(f'{N}LL, {g}LL, {q}, {primitive_root(q)}')
from polynomial import Polynomial
from polynomial_factoring import find_number_of_factors
from arithmetic import generate_primes


def build_pth_cyclotomic_polynomial(p: int, c: int):
    x = Polynomial([0, 1], c)
    return (x ** p - 1) // (x - 1)


for prime in generate_primes(10000):
    try:
        pol = build_pth_cyclotomic_polynomial(11, prime)
        print(f'p = {prime}, n_factors = {11 - 1 - find_number_of_factors(pol)}')
    except:
        print(f'p = {prime}, aaaaaaaaaaaaaaaaaa')

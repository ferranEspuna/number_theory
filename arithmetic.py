from tqdm import tqdm
import itertools

from polynomial import Polynomial
from polynomial_factoring import is_irreducible


def generate_primes(limit=float('inf')):
    yield 2
    bar = tqdm(total=limit)
    bar.update(2)
    primes = [2]
    n = 3
    while True:

        if n > limit:
            break

        bar.update(1)
        for p in primes:
            if n % p == 0:
                break
            if p * p >= n:
                primes.append(n)
                yield n
                break
        n += 1


def generate_primes_in_zpx(p, degree):
        if degree > 0:
            for x in range(1, p):
                yield Polynomial([0, x], p)

        a = [range(0, p)] * degree + [range(1, p)]
        for coefs in tqdm(itertools.product(*a), total=p ** degree):
            try:
                pol = Polynomial(list(coefs)[::-1], p)
                if is_irreducible(pol) == 1:
                    yield pol
            except:
                pass


if __name__ == '__main__':

    for p in generate_primes(100):
        print(p)

    for p in generate_primes_in_zpx(2, 19):
        print(p)

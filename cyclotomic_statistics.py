import multiprocessing
import time

from tqdm import tqdm

from polynomial import Polynomial
from polynomial_factoring import find_one_irreducible, find_splitting_degree
from arithmetic import generate_primes, generate_primes_in_zpx, generate_batched_primes
from collections import Counter
import multiprocessing as mp


def build_pth_cyclotomic_polynomial(p: int, c: int):
    x = Polynomial([0, 1], c)
    return (x ** p - 1) // (x - 1)


counter = Counter()
counter_total = Counter()
base_prime = 11


def worker(primes):
    c1 = Counter()
    c2 = Counter()
    print(f'my assigned primes are {primes[0]}:{primes[-1]}')

    for prime in primes:
        try:
            pol = build_pth_cyclotomic_polynomial(base_prime, prime)
            # degree = find_one_irreducible(pol).degree()
            degree = find_splitting_degree(pol)
            remainder = prime % base_prime
            n_factors, res = divmod((base_prime - 1), degree)
            assert res == 0
            c1[(n_factors, remainder)] += 1
            c2[n_factors] += 1

        except Exception as e:
            print(e)
            print(prime)
            assert prime == base_prime

    return c1, c2, primes[0], primes[-1]


print(multiprocessing.cpu_count())
with mp.Pool() as pool:
    t0 = time.time()

    for cfactors, cremainder, first_prime, last_prime in pool.imap(worker, generate_batched_primes(100, 200_000)):
        counter += cfactors
        counter_total += cremainder
        print(first_prime, last_prime)
        print(counter)
        print(counter_total)

    t1 = time.time()

print(t1 - t0)

print(counter)
print(counter_total)

"""for prime in generate_primes(10_000):

    try:
        pol = build_pth_cyclotomic_polynomial(base_prime, prime)
        #degree = find_one_irreducible(pol).degree()
        degree = find_splitting_degree(pol)
        remainder = prime % base_prime
        n_factors, res = divmod((base_prime - 1), degree)
        assert res == 0
        counter[(n_factors, remainder)] += 1
        counter_total[n_factors] += 1

    except Exception as e:
        print(e)
        print(prime)
        assert prime == base_prime



print(counter)
print(counter_total)"""

"""x = Polynomial([0, 1], 2)
mypol = x**2 + x + 1

mycount = Counter()

for p in generate_primes_in_zpx(2, 19):

    residue = p % mypol
    mycount[str(residue)] += 1

print(mycount)"""

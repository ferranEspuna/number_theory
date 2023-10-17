import numpy as np
from polynomial import Polynomial, gcd
from linalg import reduce, find_some_non_trivial_in_reduced

def build_T_minus_I_matrix(p: Polynomial):
    assert gcd(p, p.differentiate()).degree() == 0, 'Polynomial must be separable!'
    c = p.char
    n = p.degree()

    powers = [(Polynomial([0, 1], c) ** (c * i) -
               Polynomial([0, 1], c) ** i
               ) % p
              for i in range(n)]

    matrix = np.zeros((n, n), dtype=int)
    for i, q in enumerate(powers):
        d = q.degree()
        matrix[:d + 1, i] = q.coefficients

    return matrix


def find_nontrivial_factorization(f: Polynomial, verbose=False):
    c = f.char

    if verbose:
        print(f'Checking irreducibility of f(x) = {f} over F_{c}...\n')

    M = build_T_minus_I_matrix(f)

    if verbose:
        print(f'the relevant matrix is:\n{M}\n')

    rank = reduce(M, c)

    if verbose:
        print(f'The reduced matrix is:\n{M}\n')
        print(f'It has rank {rank}.\n')

    if rank == M.shape[0] - 1:
        if verbose:
            print('The matrix is (almost) full rank, so the polynomial is irreducible.\n')
        return None

    coef = find_some_non_trivial_in_reduced(M, c)


    assert coef is not None, 'Something went wrong. No non-trivial solution found.'

    h = Polynomial(coef, c)
    assert h.degree() > 0, 'Something went wrong. No non-trivial solution found.'

    must_be_zero = ((h ** c) - h) % f
    assert must_be_zero == 0, 'This is not an actual solution!'

    if verbose:
        print(f'A non-trivial solution is h(x) = {h}\n')

    for a in range(c):
        f1 = gcd(f, h - a).monic()
        if f1.degree() > 0:
            if verbose:
                print(f'Found a non-trivial factor: f1 = gcd(f, h - {a}) = {f1}')

            f2 = f // f1

            if verbose:
                print(f'And f2 = f / f1 = {f2}')

            assert f1.degree() > 0 and f2.degree() > 0
            assert f1 * f2 == f

            if verbose:
                print(f'Indeed, f1 * f2 = ({f1}) * ({f2}) = {f} = f')

            return f1, f2

    
def factor_into_irreducibles(f: Polynomial, verbose=False):

    if verbose:
        print(f'\nFactoring f(x) = {f} into irreducibles over F_{f.char}...')

    factors = []

    def factor_into_irreducibles_rec(f_: Polynomial, verbose=False):
        if f_.degree() == 1:
            factors.append(f_)
            return

        partial_frac = find_nontrivial_factorization(f_, verbose=verbose)

        if partial_frac is None:
            factors.append(f_)
            return
        else:
            f1, f2 = partial_frac

        factor_into_irreducibles_rec(f1, verbose=False)
        factor_into_irreducibles_rec(f2, verbose=False)

    factor_into_irreducibles_rec(f, verbose=verbose)
    factors.sort(reverse=True)

    if verbose:

        print(f'f(x) = {f} factors into irreducibles as',
              '\nf(x) = \n', " *\n".join([f"\t({str(f)})" for f in factors]))

    return factors


if __name__ == '__main__':

    char = input('Enter a prime number: ')
    char = int(char)
    _f = input('Enter a polynomial: ')
    _f = Polynomial(_f, char)

    factor_into_irreducibles(_f, verbose=True)


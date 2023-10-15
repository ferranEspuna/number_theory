from copy import copy
from typing import Any, List

import numpy as np


class Polynomial:

    @classmethod
    def get_term(cls, term: str):

        term = term.strip()

        if term == '':
            return 0, 0

        if 'x' not in term:
            return int(term.strip()), 0

        parts = f' {term} '.split('x')
        if len(parts) > 2:
            raise ValueError(f'Invalid term: {term}')

        coef, exp = parts

        coef = coef.strip()
        final_coef = 1
        for i in range(len(coef), 0, -1):
            try:
                final_coef = int(coef[:i])
                break
            except:
                pass

        exp = exp.strip()
        final_exp = 1
        for i in range(len(exp)):
            try:
                final_exp = int(exp[i:])
                break
            except:
                pass

        return final_coef, final_exp

    @classmethod
    def print_term(cls, coef: int, exp: int):

        if coef == 0:
            return ''

        out = ''

        if coef != 1:
            out += f'{coef}'

        if exp != 0:
            out += 'x'
        if exp > 1:
            out += f'^{exp}'
        if out == '':
            out = '1'

        return out

    def __init__(self, coefficients: Any, characteristic: int = 5):
        self.char = characteristic
        if isinstance(coefficients, np.ndarray):
            assert len(coefficients.shape) == 1
            assert coefficients.dtype == int, 'Polynomial coefficients must be integers'
            self.coefficients = list(coefficients)
        elif isinstance(coefficients, List):
            assert all([isinstance(c, int) for c in coefficients]), 'Polynomial coefficients must be integers'
            self.coefficients = coefficients
        elif isinstance(coefficients, str):
            terms_prev = coefficients.split('+')
            terms = []

            for t in terms_prev:
                sign = 1
                for final_term in t.split('-'):
                    c, e = self.get_term(final_term)
                    terms.append(((sign * c) % self.char, e))
                    sign = -1

            self.coefficients = [0] * (max([t[1] for t in terms]) + 1)
            for coef, exp in terms:
                self.coefficients[exp] += coef
                self.coefficients[exp] %= self.char

    def __call__(self, x):
        y = 0
        for c in self.coefficients[1:][::-1]:
            y += c
            y %= self.char
            y *= x
            y %= self.char
        y += self.coefficients[0]
        y %= self.char
        return y

    def __add__(self, other):

        if isinstance(other, int):
            other = Polynomial([other], self.char)

        assert self.char == other.char, 'Polynomials must have the same characteristic'

        if len(self.coefficients) < len(other.coefficients):
            return other + self

        result = copy(self.coefficients)
        for i, c in enumerate(other.coefficients):
            result[i] += c
            result[i] %= self.char

        degree = max(i for i, c in enumerate(result) if c != 0)
        result = result[:degree + 1]

        return Polynomial(result)

    def __mul__(self, other):

        if isinstance(other, int):
            other = Polynomial([other], self.char)

        assert self.char == other.char, 'Polynomials must have the same characteristic'

        result = [0] * (len(self.coefficients) + len(other.coefficients) - 1)
        for i, c1 in enumerate(self.coefficients):
            for j, c2 in enumerate(other.coefficients):
                result[i + j] += c1 * c2 % self.char
                result[i + j] %= self.char

        return Polynomial(result)

    def __str__(self):
        s = ''
        for i, c in enumerate(self.coefficients):

            t = self.print_term(c, i)
            if t:
                s = f'{self.print_term(c, i)} + ' + s

        if len(s) == 0:
            return '0'

        return s[:-3]


if __name__ == '__main__':
    char = 5
    p1 = Polynomial('2x^3 + 3x^2 + 7x + x^2 - 1', char)
    print('p1:', p1)
    p2 = Polynomial('3x^3 + 4x+ x^2 + 34', char)
    print('p2:', p2)

    p3 = p1 * p2
    p4 = p1 + p2

    print('p1*p2:', p3)
    print('p1+p2:', p4)
    print()

    for i in range(5):
        print(f'x = {i}')
        print(f'{p1(i)} * {p2(i)} == {p3(i)} mod {char}')
        assert (p3(i) == (p1(i) * p2(i)) % 5)
        print(f'{p1(i)} + {p2(i)} == {p4(i)} mod {char}')
        assert (p4(i) == (p1(i) + p2(i)) % 5)

        print()
    print('All tests passed!')


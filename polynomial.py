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
        for j in range(len(coef), 0, -1):
            try:
                final_coef = int(coef[:j])
                break
            except:
                pass

        exp = exp.strip()
        final_exp = 1
        for j in range(len(exp)):
            try:
                final_exp = int(exp[j:])
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

    def reduce(self):

        i = len(self.coefficients) - 1
        while i >= 0 and self.coefficients[i] == 0:
            i -= 1
        self.coefficients = self.coefficients[:i + 1]

    def __init__(self, coefficients: Any, characteristic: int) -> object:
        self.char = characteristic
        if isinstance(coefficients, np.ndarray):
            assert len(coefficients.shape) == 1
            self.coefficients = list(coefficients)
        elif isinstance(coefficients, List):
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

            self.reduce()

    def copy(self):
        return Polynomial(copy(self.coefficients), self.char)

    def degree(self):
        self.reduce()
        return len(self.coefficients) - 1

    def lead(self):
        self.reduce()
        return 0 if len(self.coefficients) == 0 else self.coefficients[-1]

    def __eq__(self, other):
        if isinstance(other, int):
            other = Polynomial([other], self.char)

        self.reduce()
        other.reduce()
        return self.coefficients == other.coefficients and self.char == other.char

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

        p = Polynomial(result, self.char)
        p.reduce()
        return p

    def __neg__(self):
        return Polynomial([0 if c == 0 else self.char - c for c in self.coefficients], self.char)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):

        if isinstance(other, int):
            other = Polynomial([other], self.char)

        assert self.char == other.char, 'Polynomials must have the same characteristic'

        # todo maybe use fft to speed up
        result = [0] * (len(self.coefficients) + len(other.coefficients) - 1)
        for i, c1 in enumerate(self.coefficients):
            for j, c2 in enumerate(other.coefficients):
                result[i + j] += c1 * c2 % self.char
                result[i + j] %= self.char

        return Polynomial(result, self.char)

    def __divmod__(self, other):
        if not isinstance(other, Polynomial):
            other = Polynomial([other], self.char)
        assert self.char == other.char, 'Polynomials must have the same characteristic'

        d = other.degree()

        if self.degree() < d:
            return Polynomial([0], self.char), self

        if other == 0:
            raise ZeroDivisionError('Polynomial division by zero')

        q = Polynomial([0], self.char)
        r = self.copy()

        l = other.lead()
        l_inv = pow(int(l), -1, self.char)

        while True:

            rd = r.degree()

            if rd < d:
                break

            c = (r.lead() * l_inv) % self.char
            e = rd - d
            t = Polynomial([0] * e + [c], self.char)
            q += t
            r -= t * other

        return q, r

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __pow__(self, power, modulo=None):
        if modulo is not None:
            raise NotImplementedError('Modular exponentiation not implemented')
        if power == 0:
            return Polynomial([1], self.char)
        if power == 1:
            return self.copy()

        q, r = divmod(power, self.char)
        if q > 0:
            ez_coef = [0] * (q * self.degree() * self.char + 1)
            for i, c in enumerate(self.coefficients):
                ez_coef[q * i * self.char] = c

            ez = Polynomial(ez_coef, self.char)
            return ez * pow(self, r)

        if power == 2:
            return self * self

        if power % 2 == 0:
            return pow(pow(self, power // 2), 2)
        else:
            return self * pow(self, power - 1)

    def differentiate(self):
        result = [0] * self.degree()
        for i, c in enumerate(self.coefficients[1:]):
            result[i] = ((i + 1) * c) % self.char

        return Polynomial(result, self.char)

    def monic(self):
        return self // self.lead()

    def __str__(self):
        s = ''
        for i, c in enumerate(self.coefficients):

            t = self.print_term(c, i)
            if t:
                s = f'{self.print_term(c, i)} + ' + s

        if len(s) == 0:
            return '0'

        return s[:-3]


def gcd(a: Polynomial, b: Polynomial) -> Polynomial:
    assert a.char == b.char, 'Polynomials must have the same characteristic'
    if b == 0:
        return a
    return gcd(b, a % b)


if __name__ == '__main__':
    char = 7
    p1 = Polynomial('2x^3 + 3x^2 + 7x + x^2 - 1', char)
    print('p1:', p1)
    p2 = Polynomial('3x^3 + 4x+ x^2 + 34', char)
    print('p2:', p2)

    p3 = p1 * p2
    p4 = p1 + p2

    print('p1*p2:', p3)
    print('p1+p2:', p4)
    print()

    for i in range(char):
        print(f'x = {i}')
        print(f'{p1(i)} * {p2(i)} == {p3(i)} mod {char}')
        assert (p3(i) == (p1(i) * p2(i)) % char)
        print(f'{p1(i)} + {p2(i)} == {p4(i)} mod {char}')
        assert (p4(i) == (p1(i) + p2(i)) % char)

        print()
    print('All tests passed!')

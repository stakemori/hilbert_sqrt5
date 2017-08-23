# -*- coding: utf-8; mode: sage -*-
from sage.all import QuadraticField, PolynomialRing, ZZ, flatten, gcd

K = QuadraticField(5)
R = PolynomialRing(K, names='g2, g5, g6')


def gens():
    return R.gens()


def to_pol(l):
    sqrt5 = K.gen()
    g2, g5, g6 = gens()
    return sum(g2**a * g5**b * g6**c * (ZZ(s) + ZZ(t) * sqrt5) / ZZ(2) for (a, b, c), s, t in l)


def to_pols_normalized(ll):
    pols = [R(to_pol(l)) for l in ll]
    l = flatten([a.coefficients() for a in pols])
    idl = K.ideal(l)
    a = idl.gens_reduced()[0]
    return [p / a for p in pols]


def syzygy(f, g):
    f = R(f)
    g = R(g)
    d = gcd(f.lt(), g.lt())
    return (f.lt() / d) * g - (g.lt() / d) * f

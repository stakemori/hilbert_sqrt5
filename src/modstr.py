# -*- coding: utf-8; mode: sage -*-
from sage.all import QuadraticField, PolynomialRing, ZZ, flatten, gcd

K = QuadraticField(5)
R = PolynomialRing(K, names='g2, g5, g6')
g2, g5, g6 = R.gens()


def gens():
    return R.gens()


def to_pol(l):
    sqrt5 = K.gen()
    return sum(g2**a * g5**b * g6**c * (ZZ(s) + ZZ(t) * sqrt5) / ZZ(2) for (a, b, c), s, t in l)


def to_pol1(tpl):
    sqrt5 = K.gen()
    l, dnm = tpl
    denm = (ZZ(dnm[0]) + ZZ(dnm[1]) * sqrt5)
    return sum(g2**a * g5**b * g6**c * (ZZ(s) + ZZ(t) * sqrt5) /
               ZZ(2) for (a, b, c), (s, t) in l) / denm


def to_pols_normalized(ll):
    pols = [R(to_pol(l)) for l in ll]
    l = flatten([a.coefficients() for a in pols])
    idl = K.ideal(l)
    a = idl.gens_reduced()[0]
    return [p / a for p in pols]


def tpl_to_rust_code(tpl, names):
    if len(tpl) == 1:
        return "into(%s) * %s" % (pol_to_rust_code(tpl[0]), names[0])
    else:
        l = [(pol_to_rust_code(a), b) for a, b in zip(tpl, names)]
        return " + ".join(["into(%s) * %s" % (a, b) for a, b in l if a])


def pol_to_rust_code(pl):
    if pl == 0:
        return None

    def gmp_code(i):
        if i == 1:
            return None
        elif i.abs() < 4294967296:
            return "(%s)" % i
        else:
            return '&Mpz::from_str_radix("%s", 10).unwrap()' % i

    def pow_code(f, i):
        if i == 0:
            return None
        if i == 1:
            return f
        else:
            return "%s.pow(%s)" % (f, i)

    def term(a, b, c, v):
        l = [pow_code("g2", a), pow_code("g5", b), pow_code("g6", c), gmp_code(ZZ(v))]
        return " * ".join([x for x in l if x is not None])
    return " + ".join([term(a, b, c, v) for (a, b, c), v in pl.dict().items()])


def syzygy(f, g):
    f = R(f)
    g = R(g)
    d = gcd(f.lt(), g.lt())
    return (f.lt() / d) * g - (g.lt() / d) * f

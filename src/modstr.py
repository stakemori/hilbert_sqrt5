# -*- coding: utf-8; mode: sage -*-
from itertools import takewhile

from sage.all import (ZZ, FreeModule, PolynomialRing, QuadraticField,
                      cached_method, flatten, gcd, load)
from sage.libs.singular.function import singular_function

K = QuadraticField(5)
R = PolynomialRing(K, names='g2, g5, g6')
g2, g5, g6 = R.gens()


smodule = singular_function("module")
sideal = singular_function("ideal")
squotient = singular_function("quotient")
smres = singular_function("mres")
slist = singular_function("list")
sintersect = singular_function("intersect")
ssyz = singular_function("syz")


def load_wts_brs(i):
    brs = load("/home/sho/work/rust/hilbert_sqrt5/data/brackets/str%s_brs.sobj" % i)
    wts = load("/home/sho/work/rust/hilbert_sqrt5/data/brackets/str%s_weights.sobj" % i)
    return FormsData(wts, [to_pol_over_z(p) for p in brs])


class FormsData(object):

    def __init__(self, weights, brackets):
        self._forms = list(enumerate(weights))
        self._brackets_dict = {}
        keys = [(a, b) for a in self._forms for b in self._forms if a[0] < b[0]]
        self._brackets_dict = {(a, b): br for (a, b), br in zip(keys, brackets)}
        for a in self._forms:
            for b in self._forms:
                if a[0] > b[0]:
                    self._brackets_dict[(a, b)] = -self._brackets_dict[(b, a)]

    @property
    def forms(self):
        return self._forms

    @property
    def brackets_dict(self):
        return self._brackets_dict

    def num_coeffs(self, x, basis):
        f, g = basis
        d = self.brackets_dict
        return [-d[(g, x)], d[(f, x)]]

    # Up to 50 this returns some.
    @cached_method
    def relatively_prime_3forms_maybe(self):
        l = ((f, g, h) for f in self.forms for g in self.forms if f[0] < g[0]
             for h in self.forms if g[0] < h[0]
             if R(gcd(self.brackets_dict[(f, g)],
                      self.brackets_dict[(f, h)],)).degree() == 0)
        for a in l:
            return a
        return None

    def relatively_prime_4forms_maybe(self):
        l = ((f, g, h, j) for f in self.forms for g in self.forms if f[0] < g[0]
             for h in self.forms if g[0] < h[0]
             for j in self.forms if h[0] < j[0]
             if R((gcd(self.brackets_dict[(f, g)],
                       self.brackets_dict[(h, j)],))).degree() == 0)
        for t in l:
            return t
        return None


def min_reol_maybe_with3gens(data):
    forms = data.relatively_prime_3forms_maybe()
    d = data.brackets_dict
    if forms is None:
        return None
    F = FreeModule(R, 2)
    f, g, h = forms
    e0, e1 = F.gens()
    a = d[(g, h)]
    b = -d[(f, h)]
    c = d[(f, g)]
    f = e0 * c
    g = e1 * c
    h = -(a * e0 + b * e1)
    n = smodule(f, g, h, b * e0, b * e1)
    sn = ssyz(n)
    m = apply(smodule, [F(x[-2] * e0 + x[-1] * e1) for x in sn])
    return list(takewhile(lambda l: any(x != 0 for x in l), slist(smres(m, 0))))


def gens():
    return R.gens()


def to_pol(l):
    sqrt5 = K.gen()
    return sum(g2**a * g5**b * g6**c * (ZZ(s) + ZZ(t) * sqrt5) / ZZ(2) for (a, b, c), s, t in l)


def to_pol_over_z(tpl):
    l, dnm = tpl
    dnm = ZZ(dnm)
    return sum(g2**a * g5**b * g6**c * ZZ(s) / ZZ(2) for (a, b, c), s in l) / dnm


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

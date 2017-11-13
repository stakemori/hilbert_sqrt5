# -*- coding: utf-8; mode: sage -*-
from itertools import takewhile
from pickle import Pickler
from os.path import join
from sage.all import (ZZ, FreeModule, PolynomialRing, QuadraticField,
                      TermOrder, cached_method, flatten, gcd, load, QQ, cached_function)
from sage.libs.singular.function import singular_function

K = QuadraticField(5)
Monomial_Wts = (6, 5, 2)
R = PolynomialRing(K, names='g6, g5, g2', order=TermOrder('wdegrevlex', Monomial_Wts))
g6, g5, g2 = R.gens()

DATA_DIR = "/home/sho/work/rust/hilbert_sqrt5/data/brackets"

smodule = singular_function("module")
sideal = singular_function("ideal")
squotient = singular_function("quotient")
smres = singular_function("mres")
slist = singular_function("list")
sintersect = singular_function("intersect")
ssyz = singular_function("syz")


def diag_res(f):
    R_el = PolynomialRing(QQ, "E4, E6")
    E4, E6 = R_el.gens()
    Delta = (E4**3 - E6**2) / 1728
    d = {g2: E4, g5: 0, g6: 2 * Delta}
    return f.subs(d)


@cached_function
def load_wts_brs(i):
    brs = load(join(DATA_DIR, "str%s_brs.sobj" % i))
    wts = load(join(DATA_DIR, "str%s_weights.sobj" % i))
    return FormsData(wts, [to_pol_over_z(p) for p in brs])


def degree(p):
    p = R(p)
    return int(p.weighted_degree(Monomial_Wts))


def degree_vec(v, wts):
    return next(int(degree(p) + w) for p, w in zip(v, wts) if p != 0)


def to_unicode(a):
    return unicode(str(a), 'utf-8')


def to_string_monom_formal(pl):
    pl = R(pl)
    pl = pl.change_ring(ZZ)
    return [((k[0], k[1], k[2]), to_unicode(a)) for (k, a) in pl.dict().items()]


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

    # Up to 50 this returns some.
    @cached_method
    def relatively_prime_3forms_maybe(self):
        l = ((f, g, h) for f in self.forms for g in self.forms if f[0] < g[0]
             for h in self.forms if g[0] < h[0]
             if self.brackets_dict[(f, g)] != 0 and
             R(gcd(self.brackets_dict[(f, g)],
                   self.brackets_dict[(f, h)],)).degree() == 0)
        for a in l:
            return a
        return None

    def relatively_prime_4forms_maybe(self):
        l = ((f, g, h, j) for f in self.forms for g in self.forms if f[0] < g[0]
             for h in self.forms if g[0] < h[0]
             for j in self.forms if h[0] < j[0]
             if self.brackets_dict[(f, g)] != 0 and
             R((gcd(self.brackets_dict[(f, g)],
                    self.brackets_dict[(h, j)],))).degree() == 0)
        for t in l:
            return t
        return None

    def weight_of_basis(self):
        f, g, _ = self.relatively_prime_3forms_maybe()
        c = self.brackets_dict[(f, g)]
        c_deg = degree(c)
        return (f[1] - c_deg, g[1] - c_deg)


def min_resol_to_primitive(m_rel):
    def to_string_monoms(l):
        return [[to_string_monom_formal(p) for p in v] for v in l]
    return ([to_string_monoms(l) for l in m_rel[0]],
            m_rel[1],
            (m_rel[2][0], m_rel[2][1], to_string_monom_formal(m_rel[2][2])))


@cached_function
def load_star_norms(i):
    return [to_pol_over_z(d) for d in load(join(DATA_DIR, "str%s_star_norms.sobj" % i))]


@cached_function
def load_cand_dnm(i):
    l = load(join(DATA_DIR, "str%s_cand.sobj" % i))
    return to_pol_over_z_wo_dnm(l[-1][-1])


@cached_function
def load_cand_wts(i):
    l = load(join(DATA_DIR, "str%s_cand.sobj" % i))
    return l[1]


def check_construction(i):
    l = load_star_norms(i)
    assert all(a != 0 for a in l)
    dnm = load_cand_dnm(i)
    assert dnm != 0
    return all((dnm**2).divides(a) for a in l)


def save_min_resol_prim(i):
    data = load_wts_brs(i)
    resl = min_reol_maybe_with3gens(data)
    resl_prim = min_resol_to_primitive(resl)
    fname = join(DATA_DIR, "str%s_cand.sobj" % i)
    with open(fname, "w") as fp:
        Pickler(fp, 2).dump(resl_prim)


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
    n = smodule(f, g, h)
    idl = sideal(b)
    m = squotient(n, idl)
    wts = data.weight_of_basis()
    mls = list(takewhile(lambda l: any(x != 0 for x in l), slist(smres(m, 0))))
    mls = [[list(v) for v in list(l)] for l in mls]
    wts_of_mls = []
    wts_of_mls.append([degree_vec(v, wts) for v in mls[0]])
    for i, l in enumerate(mls[1:]):
        wts = wts_of_mls[i]
        wts_of_mls.append([degree_vec(v, wts) for v in l])
    return (mls, wts_of_mls, (forms[0], forms[1], c))


def gens():
    return R.gens()


def to_pol(l):
    sqrt5 = K.gen()
    return sum(g2**a * g5**b * g6**c * (ZZ(s) + ZZ(t) * sqrt5) / ZZ(2) for (a, b, c), s, t in l)


def to_pol_over_z(tpl):
    l, dnm = tpl
    dnm = ZZ(dnm)
    return sum(g2**a * g5**b * g6**c * ZZ(s) / ZZ(2) for (a, b, c), s in l) / dnm


def to_pol_over_z_wo_dnm(l):
    return sum(g2**a * g5**b * g6**c * ZZ(s) / ZZ(2) for (a, b, c), s in l)


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


def c_km2_1_01(k):
    k1, k2 = [a % 4 for a in k]
    if (k1, k2) in [(0, 0), (2, 2)]:
        return 1
    if k1 % 2 == 1 or k2 % 2 == 1:
        return 0
    if (k1, k2) in [(0, 2), (2, 0)]:
        return -1


def c_km2_1_11(k):
    k1, k2 = [a % 3 for a in k]
    if (k1, k2) in [(0, 0), (2, 2)]:
        return 1
    if k1 == 1 or k2 == 1:
        return 0
    if (k1, k2) in [(0, 2), (2, 0)]:
        return -1


def c_km2_1_rho(k, rho):
    k1, k2 = [a % 5 for a in k]
    k = (k1, k2)
    if k in [(0, 0), (2, 2), (4, 3), (3, 4)]:
        return 1
    if k in [(3, 0), (4, 2)]:
        return rho
    if k in [(2, 4), (0, 3)]:
        return rho.galois_conjugate()
    if k1 == 1 or k2 == 1:
        return 0
    if k in [(2, 3), (0, 4)]:
        return -(rho.galois_conjugate())
    if k in [(3, 2), (4, 0)]:
        return -rho
    if k in [(0, 2), (2, 0), (3, 3), (4, 4)]:
        return -1


def dimension_cuspforms_sqrt5(k1, k2):
    '''
    Return the dimension of hilbert cusp forms of
    weight (k1, k2) where k1 > 2 and k2 > 2.

    cf. K. Takase, on the trace formula of the Hecke operators and the special
    values of the second L-functions attached to the Hilbert modular forms.
    manuscripta math. 55, 137 --  170 (1986).
    '''
    k = (k1, k2)
    F = QuadraticField(5)
    rho = (1 + F.gen()) / 2
    a = c_km2_1_rho(k, rho)
    return ((k1 - 1) * (k2 - 1) / 60 + c_km2_1_01(k) / 4 + c_km2_1_11(k) / 3 +
            (a + F(a).galois_conjugate()) / 5)

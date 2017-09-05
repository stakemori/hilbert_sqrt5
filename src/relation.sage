from pickle import Pickler

K = QuadraticField(5)
sqrt5 = K.gen()


def to_unicode(a):
    return unicode(str(a), 'utf-8')


def relation(vv):
    vv = [[K(ZZ(a) + ZZ(b) * sqrt5) / ZZ(2) for a, b in w] for w in vv]
    m = matrix(vv)
    try:
        rels = m.left_kernel().basis()
        rels_normlalized = []
        for rel in rels:
            rel = [K(a) for a in rel]
            I = K.ideal(rel)
            a = I.gens_reduced()[0]
            rel = [b / a for b in rel]
            rel_t = [(to_unicode(2 * b[0]), to_unicode(2 * b[1])) for b in rel]
            rels_normlalized.append(rel_t)
        with open(data_name, "w") as f:
            Pickler(f, 2).dump(rels_normlalized)
    finally:
        return m.rank()


print relation(load(data_name))

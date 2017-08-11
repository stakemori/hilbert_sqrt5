K = QuadraticField(5)
sqrt5 = K.gen()

data_name = "./data/rust_python_data.sobj"
data_res_name = "./data/rust_python_data_res.sobj"


def relation(tpl):
    v, vv = tpl
    v = [(ZZ(a) + ZZ(b) * sqrt5) / ZZ(2) for a, b in v]
    vv = [[K(ZZ(a)) for a in w] for w in vv]
    vv.insert(0, v)
    m = matrix(vv)
    rel = m.left_kernel().basis()[0]
    rel = [K(a) for a in rel]
    I = K.ideal(rel)
    a = I.gens_reduced()[0]
    rel = [b / a for b in rel]
    rel_t = [(unicode(str(2 * b[0]), 'utf-8'), unicode(str(2 * b[1]), 'utf-8')) for b in rel]
    save(rel_t, data_res_name)
    return m.nrows() - m.rank()


print relation(load(data_name))

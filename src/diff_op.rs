use elements::HmfGen;
use misc::Sqrt5Elt;
use gmp::mpz::Mpz;
use libc::{c_ulong, c_long};
use std::cmp::min;
use eisenstein::eisenstein_series;
use theta_chars::g5_normalized;
use misc::PowGen;
use fcvec;

impl PartialEq for Sqrt5Elt<Mpz> {
    fn eq(&self, other: &Self) -> bool {
        self.rt == other.rt && self.ir == other.ir
    }
}

impl Sqrt5Elt<Mpz> {
    pub fn new() -> Self {
        Sqrt5Elt {
            rt: Mpz::new(),
            ir: Mpz::new(),
        }
    }

    pub fn set_mul_mut(&mut self, other: &Self, tmp: &mut Mpz) {
        let &mut Sqrt5Elt {
            rt: ref mut a,
            ir: ref mut b,
        } = self;
        let &Sqrt5Elt {
            rt: ref c,
            ir: ref d,
        } = other;
        tmp.set(b);
        *b *= c;
        b.addmul_mut(a, d);
        *b >>= 1;

        *tmp *= 5 as c_ulong;
        *a *= c;
        a.addmul_mut(tmp, d);
        *a >>= 1;
    }

    fn set(&mut self, other: &Self) {
        let &mut Sqrt5Elt {
            rt: ref mut a,
            ir: ref mut b,
        } = self;
        let &Sqrt5Elt {
            rt: ref c,
            ir: ref d,
        } = other;
        a.set(c);
        b.set(d);
    }

    pub fn pow_mut(&mut self, expt: &Vec<char>, tmp_elt: &mut Self, tmp: &mut Mpz) {
        tmp_elt.set(self);
        self.set_one();
        for &c in expt.iter() {
            if c == '0' {
                tmp_elt.square(tmp);
            } else if c == '1' {
                self.set_mul_mut(tmp_elt, tmp);
                tmp_elt.square(tmp);
            }
        }
    }

    fn set_one(&mut self) {
        self.rt.set_ui(2);
        self.ir.set_ui(0);
    }

    /// self = self^2
    fn square(&mut self, tmp: &mut Mpz) {
        let &mut Sqrt5Elt {
            rt: ref mut a,
            ir: ref mut b,
        } = self;
        tmp.set(b);
        *b *= a as &Mpz;
        a.set_pow_ui(2);
        tmp.set_pow_ui(2);
        *tmp *= 5 as c_ulong;
        *a += tmp as &Mpz;
        *a >>= 1;
    }

    pub fn minus_norm_mut(self: &Self, res: &mut Mpz) {
        let &Sqrt5Elt {
            rt: ref a,
            ir: ref b,
        } = self;
        res.mul_mut(b, b);
        *res *= -5 as c_long;
        res.addmul_mut(a, a);
        *res >>= 2;
        res.negate();
    }
}


fn expt_to_chars(expt: usize) -> Vec<char> {
    let s = format!("{:b}", expt);
    let bts = s.into_bytes();
    let strs = bts.iter().rev().map(|&i| i as char).collect::<Vec<char>>();
    strs
}


pub fn monom_g2_g6_g10(prec: usize, expt1: usize, expt2: usize, expt3: usize) -> HmfGen {
    let mut tmp = HmfGen::new(prec);
    let mut res = HmfGen::new(prec);
    let g2 = eisenstein_series(2, prec);
    let g6 = eisenstein_series(6, prec);
    let mut g10 = g5_normalized(prec);
    g10.square();
    res.pow_mut(&g2, expt1);
    tmp.pow_mut(&g6, expt2);
    if expt2 != 0 {
        res *= &tmp;
    }
    tmp.pow_mut(&g10, expt3);
    if expt3 != 0 {
        res *= &tmp;
    }
    res
}

pub fn g15_normalized(prec: usize) -> HmfGen {
    let g2 = eisenstein_series(2, prec);
    let g6 = eisenstein_series(6, prec);
    let g5 = g5_normalized(prec);
    let mut tmp1 = HmfGen::new(prec);
    let mut tmp2 = HmfGen::new(prec);
    g15_term(&mut tmp1, 2, &g2, &g5, &g6);
    g15_term(&mut tmp2, 5, &g5, &g2, &g6);
    tmp1 -= &tmp2;
    g15_term(&mut tmp2, 6, &g6, &g2, &g5);
    tmp1 += &tmp2;
    let a = Mpz::from_si(-86400);
    assert!(tmp1.is_divisible_by_const(&a));
    tmp1 /= &a;
    tmp1
}

fn g15_term(res: &mut HmfGen, wt: c_ulong, f: &HmfGen, g: &HmfGen, h: &HmfGen) {
    diff_mul_mut_ir(res, (1, 0), (0, 1), &g, &h);
    *res *= wt;
    *res *= f;
}

/// Rt part of ∂^(a1 + a2)/(∂^a1∂^a2)f1 * ∂^(b1 + b2)/(∂^b1∂^b2)f2, where
/// (a1, a2) = expt1, (b1, b2) = expt2.
#[allow(dead_code)]
fn diff_mul_mut_rt(
    res: &mut HmfGen,
    expt1: (usize, usize),
    expt2: (usize, usize),
    f1: &HmfGen,
    f2: &HmfGen,
) {
    let prec = f1.prec;
    let ref mut tmp_f1 = HmfGen::new(prec);
    let ref mut tmp_f2 = HmfGen::new(prec);
    let ref mut tmp_f3 = HmfGen::new(prec);
    diff_mut(tmp_f1, res, expt1, &f1);
    diff_mut(tmp_f2, tmp_f3, expt2, &f2);
    *res *= 5 as c_ulong;
    *res *= tmp_f3 as &HmfGen;
    tmp_f3.mul_mut(tmp_f1, tmp_f2);
    *res += tmp_f3;
    *res >>= 1;
}

/// Similar to diff_mul_mut_rt for ir part.
fn diff_mul_mut_ir(
    res: &mut HmfGen,
    expt1: (usize, usize),
    expt2: (usize, usize),
    f1: &HmfGen,
    f2: &HmfGen,
) {
    let prec = f1.prec;
    let ref mut tmp_f1 = HmfGen::new(prec);
    let ref mut tmp_f2 = HmfGen::new(prec);
    let ref mut tmp_f3 = HmfGen::new(prec);
    diff_mut(tmp_f1, tmp_f3, expt1, &f1);
    diff_mut(tmp_f2, res, expt2, &f2);
    *res *= tmp_f1 as &HmfGen;
    tmp_f1.mul_mut(tmp_f2, tmp_f3);
    *res += tmp_f1;
    *res >>= 1;
}

fn diff_mut(res_rt: &mut HmfGen, res_ir: &mut HmfGen, expt: (usize, usize), f: &HmfGen) {
    let norm_expt = min(expt.0, expt.1);
    if norm_expt > 0 {
        diff_mut_minus_norm(res_rt, norm_expt, f);
        diff_mut_minus_norm(res_ir, norm_expt, f);
    }
    if expt.0 >= expt.1 {
        diff_mut_rt(res_rt, res_ir, expt.0 - expt.1, f);
    } else {
        diff_mut_ir(res_rt, res_ir, expt.1 - expt.0, f);
    }
}

fn diff_mut_rt(res_rt: &mut HmfGen, res_ir: &mut HmfGen, expt: usize, f: &HmfGen) {
    let mut tmp_elt = Sqrt5Elt::<Mpz>::new();
    let mut a_pow = Sqrt5Elt::<Mpz>::new();
    let mut tmp = Mpz::new();
    let expt = expt_to_chars(expt);
    v_u_bd_iter!((f.u_bds, v, u, bd) {
        a_pow.ir.set_ui(v as u64);
        a_pow.rt.set_si(u);
        a_pow.pow_mut(&expt, &mut tmp_elt, &mut tmp);
        res_rt.fcvec.fc_ref_mut(v, u, bd).mul_mut(&a_pow.rt, f.fcvec.fc_ref(v, u, bd));
        res_ir.fcvec.fc_ref_mut(v, u, bd).mul_mut(&a_pow.ir, f.fcvec.fc_ref(v, u, bd));
    }
    );
}

fn diff_mut_ir(res_rt: &mut HmfGen, res_ir: &mut HmfGen, expt: usize, f: &HmfGen) {
    let mut tmp_elt = Sqrt5Elt::<Mpz>::new();
    let mut a_pow = Sqrt5Elt::<Mpz>::new();
    let eps = if is_even!(expt) {
        1 as c_long
    } else {
        -1 as c_long
    };
    let mut tmp = Mpz::new();
    let expt = expt_to_chars(expt);
    v_u_bd_iter!((f.u_bds, v, u, bd) {
        a_pow.ir.set_ui(v as u64);
        a_pow.rt.set_si(u);
        a_pow.pow_mut(&expt, &mut tmp_elt, &mut tmp);
        a_pow.rt *= eps;
        a_pow.ir *= -eps;
        res_rt.fcvec.fc_ref_mut(v, u, bd).mul_mut(&a_pow.rt, f.fcvec.fc_ref(v, u, bd));
        res_ir.fcvec.fc_ref_mut(v, u, bd).mul_mut(&a_pow.ir, f.fcvec.fc_ref(v, u, bd));
    }
    );
}

fn diff_mut_minus_norm(res: &mut HmfGen, expt: usize, f: &HmfGen) {
    let mut a = Sqrt5Elt::<Mpz>::new();
    let mut tmp = Mpz::new();
    let expt = expt as u64;
    v_u_bd_iter!((f.u_bds, v, u, bd) {
        a.ir.set_ui(v as u64);
        a.rt.set_si(u);
        a.minus_norm_mut(&mut tmp);
        tmp.set_pow_ui(expt);
        res.fcvec.fc_ref_mut(v, u, bd).mul_mut(&tmp, f.fcvec.fc_ref(v, u, bd));
    });
}

/// set set = f/g15
pub fn divide_by_g15(res: &mut HmfGen, f: &HmfGen, g15: &HmfGen) {
    let prec = g15.prec;
    res.prec = prec - 2;
    let mut f_cloned = f.clone();
    let mut tmp = HmfGen::new(prec);
    res.decrease_prec(prec - 2);
    res.set_zero();
    let ref u_bds = f_cloned.u_bds;
    for v in 3..(prec + 1) {
        for i in 3..(v + 1) {
            fcvec::mul_mut(
                &mut tmp.fcvec.vec[v],
                &g15.fcvec.vec[i],
                &f_cloned.fcvec.vec[v - i + 2],
                i,
                v - i,
                u_bds.vec[i],
                u_bds.vec[v - i + 2],
                u_bds.vec[v],
                u_bds,
                0,
                0,
                0,
            );
            fcvec::sub_assign(&mut f_cloned.fcvec.vec[v], &tmp.fcvec.vec[v], v, u_bds);
        }
    }
    for (v, &bd1) in u_bds.vec.iter().enumerate().skip(2) {
        let v_i = v as i64;
        let bd = u_bds.vec[v-2] as i64;
        let bd1 = bd1 as i64;
        for u in u_iter!(v_i, bd) {
            res.fcvec.fc_ref_mut(v-2, u, bd).set(
                f_cloned.fcvec.fc_ref_mut(v, u, bd1)
            );
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[ignore]
    #[test]
    fn g15() {
        let f = g15_normalized(10);
        println!("{}", f);
    }

    #[test]
    fn test_mul_pow() {
        let mut a = Sqrt5Elt {
            rt: Mpz::from_ui(3),
            ir: Mpz::from_ui(5),
        };

        let b = Sqrt5Elt {
            rt: Mpz::from_ui(1),
            ir: Mpz::from_ui(1),
        };

        let mut tmp_elt = Sqrt5Elt::<Mpz>::new();
        let mut tmp = Mpz::new();

        a.set_mul_mut(&b, &mut tmp);
        assert_eq!(
            a,
            Sqrt5Elt {
                rt: Mpz::from_ui(14),
                ir: Mpz::from_ui(4),
            }
        );

        let mut a = Sqrt5Elt {
            rt: Mpz::from_ui(3),
            ir: Mpz::from_ui(5),
        };

        a.pow_mut(&expt_to_chars(10), &mut tmp_elt, &mut tmp);
        assert_eq!(
            a,
            Sqrt5Elt {
                rt: Mpz::from_ui(322355827),
                ir: Mpz::from_ui(142989825),
            }
        );

        a.minus_norm_mut(&mut tmp);
        assert_eq!(tmp.to_str_radix(10), "-420707233300201");
    }

    #[test]
    fn test_divide_byg15() {
        let prec = 15;
        let g15 = g15_normalized(prec);
        let mut e2 = eisenstein_series(2, prec);
        let f = &g15 * &e2;
        let mut res = HmfGen::new(prec - 2);
        divide_by_g15(&mut res, &f, &g15);
        e2.decrease_prec(prec - 2);
        assert_eq!(e2, res);
    }
}

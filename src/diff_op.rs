use elements::{HmfGen, weight_div};
use misc::Sqrt5Elt;
use gmp::mpz::Mpz;
use libc::{c_ulong, c_long};
use std::cmp::min;
use eisenstein::eisenstein_series;
use theta_chars::g5_normalized;
use misc::PowGen;
use fcvec;
use bignum::Sqrt5Mpz;
use bignum::BigNumber;
use std::ops::{SubAssign, MulAssign, AddAssign, ShrAssign};

impl PartialEq for Sqrt5Elt<Mpz> {
    fn eq(&self, other: &Self) -> bool {
        self.rt == other.rt && self.ir == other.ir
    }
}

/// TODO: Remove duplicate of code.
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


pub fn monom_g2_g6_g10(prec: usize, expt1: usize, expt2: usize, expt3: usize) -> HmfGen<Mpz> {
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

pub fn g15_normalized(prec: usize) -> HmfGen<Mpz> {
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
    tmp1.weight = Some((15, 15));
    tmp1
}

fn g15_term(res: &mut HmfGen<Mpz>, wt: c_ulong, f: &HmfGen<Mpz>, g: &HmfGen<Mpz>, h: &HmfGen<Mpz>) {
    diff_mul_mut_ir(res, (1, 0), (0, 1), &g, &h);
    *res *= wt;
    *res *= f;
}

/// Rt part of ∂^(a1 + a2)/(∂^a1∂^a2)f1 * ∂^(b1 + b2)/(∂^b1∂^b2)f2, where
/// (a1, a2) = expt1, (b1, b2) = expt2.
fn diff_mul_mut_rt(
    res: &mut HmfGen<Mpz>,
    expt1: (usize, usize),
    expt2: (usize, usize),
    f1: &HmfGen<Mpz>,
    f2: &HmfGen<Mpz>,
) {
    let prec = f1.prec;
    let ref mut tmp_f1 = HmfGen::new(prec);
    let ref mut tmp_f2 = HmfGen::new(prec);
    let ref mut tmp_f3 = HmfGen::new(prec);
    diff_mut(tmp_f1, res, expt1, &f1);
    diff_mut(tmp_f2, tmp_f3, expt2, &f2);
    *res *= 5 as c_ulong;
    *res *= tmp_f3 as &HmfGen<Mpz>;
    tmp_f3.mul_mut(tmp_f1, tmp_f2);
    *res += tmp_f3;
    *res >>= 1;
}

/// Similar to diff_mul_mut_rt for ir part.
fn diff_mul_mut_ir(
    res: &mut HmfGen<Mpz>,
    expt1: (usize, usize),
    expt2: (usize, usize),
    f1: &HmfGen<Mpz>,
    f2: &HmfGen<Mpz>,
) {
    let prec = f1.prec;
    let ref mut tmp_f1 = HmfGen::new(prec);
    let ref mut tmp_f2 = HmfGen::new(prec);
    let ref mut tmp_f3 = HmfGen::new(prec);
    diff_mut(tmp_f1, tmp_f3, expt1, &f1);
    diff_mut(tmp_f2, res, expt2, &f2);
    *res *= tmp_f1 as &HmfGen<Mpz>;
    tmp_f1.mul_mut(tmp_f2, tmp_f3);
    *res += tmp_f1;
    *res >>= 1;
}

fn diff_mut(
    res_rt: &mut HmfGen<Mpz>,
    res_ir: &mut HmfGen<Mpz>,
    expt: (usize, usize),
    f: &HmfGen<Mpz>,
) {
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

fn diff_mut_rt(res_rt: &mut HmfGen<Mpz>, res_ir: &mut HmfGen<Mpz>, expt: usize, f: &HmfGen<Mpz>) {
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

fn diff_mut_ir(res_rt: &mut HmfGen<Mpz>, res_ir: &mut HmfGen<Mpz>, expt: usize, f: &HmfGen<Mpz>) {
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

fn diff_mut_minus_norm(res: &mut HmfGen<Mpz>, expt: usize, f: &HmfGen<Mpz>) {
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
fn divide_by_g15<T>(res: &mut HmfGen<T>, f: &HmfGen<T>, g15: &HmfGen<T>)
where
    T: BigNumber + Clone,
    for<'a> T: SubAssign<&'a T>,
{
    let prec = f.prec;
    assert_eq!(prec, g15.prec);
    res.prec = prec - 2;
    let mut f_cloned = f.clone();
    let mut tmp = HmfGen::new(prec);
    res.decrease_prec(prec - 2);
    res.weight = weight_div(f.weight, g15.weight);
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
        let bd = u_bds.vec[v - 2] as i64;
        let bd1 = bd1 as i64;
        for u in u_iter!(v_i, bd) {
            res.fcvec.fc_ref_mut(v - 2, u, bd).set_g(
                f_cloned
                    .fcvec
                    .fc_ref_mut(v, u, bd1),
            );
        }
    }
}

#[derive(Debug)]
pub struct NotHhmError {}

/// Return <f>.
pub fn bracket_proj<T>(f: &HmfGen<T>, g15: &HmfGen<T>) -> Result<HmfGen<T>, NotHhmError>
where
    T: BigNumber + Clone + MulAssign<c_ulong> + ShrAssign<usize>,
    for<'a> T: SubAssign<&'a T>,
    for<'a> T: AddAssign<&'a T>,
{
    let mut res = HmfGen::new(f.prec - 2);
    let mut g = f.clone();
    let mut tmp = T::new_g();
    match f.weight {
        None => Err(NotHhmError {}),
        Some((k1, k2)) => {
            if is_even!((k1 + k2) >> 1) {
                for (v, &bd) in g.u_bds.vec.iter().enumerate() {
                    let bd = bd as i64;
                    let v_i = v as i64;
                    g.fcvec.fc_ref_mut(v, 0, bd).set_ui_g(0);
                    for u in u_iter_pos!(v_i, bd) {
                        tmp.set_g(g.fcvec.fc_ref(v, -u, bd));
                        *g.fcvec.fc_ref_mut(v, u, bd) -= &tmp;
                        tmp.set_g(g.fcvec.fc_ref(v, u, bd));
                        tmp.negate_g();
                        g.fcvec.fc_ref_mut(v, -u, bd).set_g(&tmp);
                    }
                }
            } else {
                for (v, &bd) in g.u_bds.vec.iter().enumerate() {
                    let bd = bd as i64;
                    let v_i = v as i64;
                    *g.fcvec.fc_ref_mut(v, 0, bd) *= 2 as c_ulong;
                    for u in u_iter_pos!(v_i, bd) {
                        tmp.set_g(g.fcvec.fc_ref(v, -u, bd));
                        *g.fcvec.fc_ref_mut(v, u, bd) += &tmp;
                        tmp.set_g(g.fcvec.fc_ref(v, u, bd));
                        g.fcvec.fc_ref_mut(v, -u, bd).set_g(&tmp);
                    }
                }

            }
            divide_by_g15(&mut res, &g, &g15);
            res >>= 1;
            Ok(res)
        }
    }
}

macro_rules! define_rankin_cohen {
    ($fun: ident, $mul_fun: ident) => {
        pub fn $fun(m: usize,
                    f: &HmfGen<Mpz>,
                    g: &HmfGen<Mpz>) -> Result<HmfGen<Mpz>, NotHhmError> {
            assert_eq!(f.prec, g.prec);
            let mut res = HmfGen::new(f.prec);
            let mut tmp = HmfGen::new(f.prec);
            let mut tmp_z = Mpz::new();
            let mut tmp_z1 = Mpz::new();
            let mut tmp_z2 = Mpz::new();
            if !f.weight.is_none() && !g.weight.is_none() {
                let (k1, k2) = f.weight.unwrap();
                let (l1, l2) = g.weight.unwrap();
                for i in 0..(m + 1) {
                    tmp_z1.set_ui((m + k2 - 1) as c_ulong);
                    tmp_z.bin_ui_mut(&tmp_z1, (m - i) as c_ulong);
                    tmp_z1.set_ui((m + l2 - 1) as c_ulong);
                    tmp_z2.bin_ui_mut(&tmp_z1, i as c_ulong);
                    tmp_z *= &tmp_z2;
                    $mul_fun(&mut tmp, (0, i), (0, m - i), &f, &g);
                    tmp *= &tmp_z;
                    res += &tmp;
                }
                res.weight = Some((k1 + l1, k2 + l2 + 2 * m));
                Ok(res)
            } else {
                Err(NotHhmError{})
            }
        }
    }
}

define_rankin_cohen!(rankin_cohen_rt, diff_mul_mut_rt);
define_rankin_cohen!(rankin_cohen_ir, diff_mul_mut_ir);
// TODO: make this function generic
pub fn rankin_cohen_sqrt5(
    m: usize,
    f: &HmfGen<Mpz>,
    g: &HmfGen<Mpz>,
) -> Result<HmfGen<Sqrt5Mpz>, NotHhmError> {
    let mut res = HmfGen::<Sqrt5Mpz>::new(f.prec);
    let res_rt = rankin_cohen_rt(m, &f, &g)?;
    let res_ir = rankin_cohen_ir(m, &f, &g)?;
    res.weight = res_rt.weight;
    v_u_bd_iter!((f.u_bds, v, u, bd) {
        res.fcvec.fc_ref_mut(v, u, bd).rt.set(res_rt.fcvec.fc_ref(v, u, bd));
        res.fcvec.fc_ref_mut(v, u, bd).ir.set(res_ir.fcvec.fc_ref(v, u, bd));
    });
    Ok(res)
}

pub fn star_op<T>(res: &mut HmfGen<T>, f: &HmfGen<T>)
where
    T: BigNumber,
{
    let (k1, k2) = f.weight.unwrap();
    v_u_bd_iter!((f.u_bds, v, u, bd) {
        res.fcvec.fc_ref_mut(v, u, bd).set_g(f.fcvec.fc_ref(v, -u, bd));
    });
    if !is_even!((k1 + k2) >> 1) {
        res.negate();
    }
    res.weight = Some((k2, k1));
}
/// Return <f g.star()>
pub fn bracket_inner_prod<T>(
    f: &HmfGen<T>,
    g: &HmfGen<T>,
    g15: &HmfGen<T>,
) -> Result<HmfGen<T>, NotHhmError>
where
    T: BigNumber + Clone + MulAssign<c_ulong> + ShrAssign<usize>,
    for<'a> T: AddAssign<&'a T>,
    for<'a> T: SubAssign<&'a T>,
{
    let mut res: HmfGen<T> = HmfGen::new(f.prec);
    star_op(&mut res, &g);
    res *= f;
    bracket_proj(&res, g15)
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
        let mut g5 = g5_normalized(prec);
        {
            let mut e4 = eisenstein_series(4, prec);
            let f = &g15 * &e4;
            let mut res = HmfGen::new(prec - 2);
            divide_by_g15(&mut res, &f, &g15);
            e4.decrease_prec(prec - 2);
            assert_eq!(res.weight, e4.weight);
            assert_eq!(e4, res);
        }

        let mut res = HmfGen::new(prec - 2);
        let g = &g15 * &g5;
        divide_by_g15(&mut res, &g, &g15);
        g5.decrease_prec(prec - 2);
        assert_eq!(g5, res);
    }
}

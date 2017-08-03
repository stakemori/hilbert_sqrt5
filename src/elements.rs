use std;
use std::fmt;
use gmp::mpz::Mpz;
use misc::{pow_mut, PowGen};
use libc::{c_ulong, c_long};
use std::ops::{AddAssign, MulAssign, DivAssign, SubAssign, ShlAssign, ShrAssign, Mul, Sub, Neg,
               Add};
use std::cmp::min;
use misc::Sqrt5Elt;

type Weight = Option<(usize, usize)>;

/// struct for hilbert modualr form over Q(sqrt(5))
/// this corresponds finite sum of the q-expansion of the form
/// Σ a(u, v) exp(2piTr 1/sqrt(5) (u + v * sqrt(5))/2)
/// where v <= prec.
/// a(u, v) = fc[v][a], where a = u + u_bds[v]
#[derive(Debug, Clone)]
pub struct HmfGen {
    pub prec: usize,
    pub fcvec: FcVec,
    pub weight: Weight,
    // vth element of u_bds.vec is (sqrt(5) * v).floor()
    pub u_bds: UBounds,
}


macro_rules! is_even {
    ($expr: expr) => {$expr & 1 == 0}
}

macro_rules! u_iter {
    ($v: expr, $bd: ident) => {
        {
            (-$bd..($bd+1)).filter(|&x| is_even!(x+$v))
        }
    }
}

macro_rules! u_iter_pos {
    ($v: expr, $bd: ident) => {
        {
            (1..($bd+1)).filter(|&x| is_even!(x+$v))
        }
    }
}

macro_rules! v_u_bd_iter {
    (($u_bds: expr, $v: ident, $u: ident, $bd: ident) $body:expr) =>
    {
        for ($v, &$bd) in $u_bds.vec.iter().enumerate() {
            let $bd = $bd as i64;
            let v_i64 = $v as i64;
            for $u in u_iter!(v_i64, $bd) {
                $body
            }
        };
    }
}

macro_rules! v_u_bd_iter_non_const {
    (($u_bds: expr, $v: ident, $u: ident, $bd: ident) $body:expr) =>
    {
        for ($v, &$bd) in $u_bds.vec.iter().enumerate().skip(1) {
            let $bd = $bd as i64;
            let v_i64 = $v as i64;
            for $u in u_iter!(v_i64, $bd) {
                $body
            }
        };
    }
}

impl fmt::Display for HmfGen {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut vec = Vec::new();
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            let a = self.fcvec.fc_ref(v, u, bd);
            if !a.is_zero() {
                vec.push(format!("({}, {}): {}", u, v, a));
            }
        }
        );
        write!(f, "{}", vec.join("\n"))
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct FcVec {
    pub vec: Vec<Vec<Mpz>>,
}

impl FcVec {
    pub fn fc_ref(&self, v: usize, u: i64, bd: i64) -> &Mpz {
        debug_assert!(u + bd >= 0);
        self.vec[v].get((u + bd) as usize).unwrap()
    }

    pub fn fc_ref_mut(&mut self, v: usize, u: i64, bd: i64) -> &mut Mpz {
        debug_assert!(u + bd >= 0);
        self.vec[v].get_mut((u + bd) as usize).unwrap()
    }

    fn new(u_bds: &UBounds) -> FcVec {
        let vec = u_bds
            .vec
            .iter()
            .map(|&bd| (0..(2 * bd + 1)).map(|_| Mpz::from_ui(0)).collect())
            .collect();
        FcVec { vec: vec }
    }
}

#[derive(Debug, Clone)]
pub struct UBounds {
    pub vec: Vec<usize>,
}

impl UBounds {
    pub fn new(prec: usize) -> UBounds {
        assert!(5 * prec * prec < std::usize::MAX);
        assert!(prec < std::i64::MAX as usize);
        let mut u_bds = Vec::new();
        let sqrt5 = 5_f64.sqrt();
        for v in 0..(prec + 1) {
            u_bds.push((sqrt5 * v as f64).floor() as usize);
        }
        UBounds { vec: u_bds }
    }

    pub fn take(&self, n: usize) -> UBounds {
        let v = self.vec.iter().map(|&x| x).take(n).collect();
        UBounds { vec: v }
    }
}

impl PowGen for HmfGen {
    fn set_one(&mut self) {
        v_u_bd_iter_non_const!((self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_ui(0);
        });
        self.fcvec.fc_ref_mut(0, 0, 0).set_ui(1);
        self.weight = Some((0, 0));
    }

    fn square(&mut self) {
        let f = self.clone();
        self.weight = weight_pow(self.weight, 2);
        let mut tmp = Mpz::from_ui(0);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            _mul_mut_tmp(&mut tmp, u, v, &f.fcvec, &f.fcvec, &self.u_bds);
            self.fcvec.fc_ref_mut(v, u, bd).set(&tmp);
            })
    }
}

fn weight_mul(a: Weight, b: Weight) -> Weight {
    a.and_then(|x| b.and_then(|y| Some((x.0 + y.0, x.1 + y.1))))
}

fn weight_add(a: Weight, b: Weight) -> Weight {
    a.and_then(|x| b.and_then(|y| if x == y { Some(x) } else { None }))
}

fn weight_pow(a: Weight, n: usize) -> Weight {
    a.and_then(|(k1, k2)| Some((k1 * n, k2 * n)))
}

pub fn weight_div(a: Weight, b: Weight) -> Weight {
    a.and_then(|x| b.and_then(|y| Some((x.0 - y.0, x.1 - y.1))))
}

impl HmfGen {
    /// Return 0 q-expantion
    pub fn new(prec: usize) -> HmfGen {
        let u_bds = UBounds::new(prec);
        let fcvec = FcVec::new(&u_bds);
        HmfGen {
            weight: None,
            prec: prec,
            fcvec: fcvec,
            u_bds: u_bds,
        }
    }

    /// Decrease prec to prec.
    pub fn decrease_prec(&mut self, prec: usize) {
        self.u_bds = self.u_bds.take(prec + 1);
        self.prec = prec;
    }

    pub fn one(prec: usize) -> HmfGen {
        let mut f = Self::new(prec);
        f.fcvec.fc_ref_mut(0, 0, 0).set_ui(1);
        f
    }

    /// set self = f1 + f2
    pub fn add_mut(&mut self, f1: &HmfGen, f2: &HmfGen) {
        let prec = min(f1.prec, f2.prec);
        self.decrease_prec(prec);
        self.weight = weight_add(f1.weight, f2.weight);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::add_mut(
                self.fcvec.fc_ref_mut(v, u, bd),
                f1.fcvec.fc_ref(v, u, bd),
                f2.fcvec.fc_ref(v, u, bd),
            );
        })
    }

    pub fn is_divisible_by_const(&self, a: &Mpz) -> bool {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            if !self.fcvec.fc_ref(v, u, bd).is_multiple_of(&a) {
                return false;
            }
        });
        true
    }

    /// set self = f1 - f2
    pub fn sub_mut(&mut self, f1: &HmfGen, f2: &HmfGen) {
        let prec = min(f1.prec, f2.prec);
        self.decrease_prec(prec);
        self.weight = weight_add(f1.weight, f2.weight);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::sub_mut(
                self.fcvec.fc_ref_mut(v, u, bd),
                f1.fcvec.fc_ref(v, u, bd),
                f2.fcvec.fc_ref(v, u, bd),
            );
        })
    }

    /// set self = f1 * f2
    pub fn mul_mut(&mut self, f1: &HmfGen, f2: &HmfGen) {
        let mut tmp = Mpz::from_ui(0);
        let prec = min(f1.prec, f2.prec);
        self.decrease_prec(prec);
        self.weight = weight_mul(f1.weight, f2.weight);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            _mul_mut_tmp(&mut tmp, u, v, &f1.fcvec, &f2.fcvec, &self.u_bds);
            self.fcvec.fc_ref_mut(v, u, bd).set(&tmp);
        })
    }

    /// self = f * a
    pub fn mul_mut_by_const(&mut self, f: &HmfGen, a: &Mpz) {
        self.weight = f.weight;
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::mul_mut(self.fcvec.fc_ref_mut(v, u, bd), f.fcvec.fc_ref(v, u, bd), a)
            })
    }

    pub fn pow_mut(&mut self, f: &HmfGen, a: usize) {
        if a == 0 {
            self.set_one();
        } else if a == 1 {
            self.set(&f);
        } else {
            pow_mut(self, f, a)
        };
        self.weight = weight_pow(f.weight, a);
    }

    pub fn is_zero(&self) -> bool {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            if !self.fcvec.fc_ref(v, u, bd).is_zero() {
                return false;
            }
        });
        true
    }

    pub fn fourier_coefficient(&self, v: usize, u: i64) -> Mpz {
        let bd = self.u_bds.vec[v] as i64;
        let mut a = Mpz::new();
        a.set(self.fcvec.fc_ref(v, u, bd));
        a
    }

    pub fn fourier_coefficients(&self, vec: &Vec<(usize, i64)>) -> Vec<Mpz> {
        vec.iter()
            .map(|&(v, u)| self.fourier_coefficient(v, u))
            .collect()
    }

    pub fn set(&mut self, other: &Self) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set(other.fcvec.fc_ref(v, u, bd));
        })
    }

    pub fn set_zero(&mut self) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_ui(0);
        })
    }

    pub fn negate(&mut self) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).negate();
        })
    }

    pub fn into_sqrt5_coeff(&self) -> Sqrt5Elt<Self> {
        let mut f = self.clone();
        f <<= 1;
        let g = Self::new(self.prec);
        Sqrt5Elt::<HmfGen> { rt: f, ir: g }
    }
}

impl<'a> DivAssign<&'a Mpz> for HmfGen {
    fn div_assign(&mut self, num: &Mpz) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            self.fcvec.fc_ref_mut(v, u, bd).set_divexact(&num);
        })
    }
}

impl<'a> AddAssign<&'a HmfGen> for HmfGen {
    fn add_assign(&mut self, other: &HmfGen) {
        self.weight = weight_add(self.weight, other.weight);
        let prec = min(self.prec, other.prec);
        self.decrease_prec(prec);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::add_assign(
                self.fcvec.fc_ref_mut(v, u, bd),
                other.fcvec.fc_ref(v, u, bd),
            );

        })
    }
}

impl<'a> SubAssign<&'a HmfGen> for HmfGen {
    fn sub_assign(&mut self, other: &HmfGen) {
        self.weight = weight_add(self.weight, other.weight);
        let prec = min(self.prec, other.prec);
        self.decrease_prec(prec);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) -= other.fcvec.fc_ref(v, u, bd);
        })
    }
}

impl<'a, 'b> Add<&'a HmfGen> for &'b HmfGen {
    type Output = HmfGen;
    fn add(self, other: &HmfGen) -> HmfGen {
        let prec = min(self.prec, other.prec);
        let mut res = HmfGen::new(prec);
        res.add_mut(self, other);
        res
    }
}

impl<'a, 'b> Mul<&'a HmfGen> for &'b HmfGen {
    type Output = HmfGen;
    fn mul(self, other: &HmfGen) -> HmfGen {
        let prec = min(self.prec, other.prec);
        let mut res = HmfGen::new(prec);
        res.mul_mut(self, other);
        res
    }
}

impl<'a, 'b> Mul<&'a Mpz> for &'b HmfGen {
    type Output = HmfGen;
    fn mul(self, other: &Mpz) -> HmfGen {
        let mut res = HmfGen::new(self.prec);
        res.mul_mut_by_const(self, other);
        res
    }
}

impl<'a, 'b> Sub<&'a HmfGen> for &'b HmfGen {
    type Output = HmfGen;
    fn sub(self, other: &HmfGen) -> HmfGen {
        let prec = min(self.prec, other.prec);
        let mut res = HmfGen::new(prec);
        res.sub_mut(self, other);
        res
    }
}

impl<'a> Neg for &'a HmfGen {
    type Output = HmfGen;
    fn neg(self) -> HmfGen {
        let mut res = HmfGen::new(self.prec);
        let mone = Mpz::from_si(-1);
        res.mul_mut_by_const(&self, &mone);
        res
    }
}

impl PartialEq for HmfGen {
    fn eq(&self, other: &HmfGen) -> bool {
        let prec = min(self.prec, other.prec);
        v_u_bd_iter!((self.u_bds.take(prec + 1), v, u, bd) {
            if self.fcvec.fc_ref(v, u, bd) != other.fcvec.fc_ref(v, u, bd) {
                return false;
            }
        }
        );
        true
    }
}


/// set (v, u) th F.C. of fc_vec1 * fc_vec2 to a.
/// This function take care the case when fc_vec2 is sparse.
fn _mul_mut_tmp(a: &mut Mpz, u: i64, v: usize, fc_vec1: &FcVec, fc_vec2: &FcVec, u_bds: &UBounds) {
    a.set_ui(0);
    let mut tmp = Mpz::new();
    for v2 in 0..(v + 1) {
        let bd2 = u_bds.vec[v2] as i64;
        let v2_i64 = v2 as i64;
        for u2 in u_iter!(v2_i64, bd2) {
            if !fc_vec2.fc_ref(v2, u2, bd2).is_zero() {
                let u1 = u - u2;
                let v1 = v - v2;
                let u1abs = u1.abs() as usize;
                if u1abs * u1abs <= 5 * v1 * v1 {
                    let bd1 = u_bds.vec[v1] as i64;
                    tmp.mul_mut(fc_vec2.fc_ref(v2, u2, bd2), fc_vec1.fc_ref(v1, u1, bd1));
                    Mpz::add_assign(a, &tmp);
                }
            }
        }
    }
}

impl<'a> MulAssign<&'a HmfGen> for HmfGen {
    fn mul_assign(&mut self, other: &HmfGen) {
        let prec = min(self.prec, other.prec);
        self.weight = weight_mul(self.weight, other.weight);
        self.decrease_prec(prec);
        // We need cloned self.
        let f = self.clone();
        let mut tmp = Mpz::from_ui(0);
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            _mul_mut_tmp(&mut tmp, u, v, &f.fcvec, &other.fcvec, &self.u_bds);
            self.fcvec.fc_ref_mut(v, u, bd).set(&tmp);
            })
    }
}

impl<'a> MulAssign<&'a Mpz> for HmfGen {
    fn mul_assign(&mut self, other: &Mpz) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::mul_assign(self.fcvec.fc_ref_mut(v, u, bd), other);
        }
        )
    }
}

impl MulAssign<c_ulong> for HmfGen {
    fn mul_assign(&mut self, other: c_ulong) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) *= other;
        }
        );
    }
}

impl MulAssign<c_long> for HmfGen {
    fn mul_assign(&mut self, other: c_long) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            *self.fcvec.fc_ref_mut(v, u, bd) *= other;
        }
        );
    }
}

impl ShlAssign<usize> for HmfGen {
    fn shl_assign(&mut self, other: usize) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::shl_assign(self.fcvec.fc_ref_mut(v, u, bd), other);
        }
        );
    }
}


impl ShrAssign<usize> for HmfGen {
    fn shr_assign(&mut self, other: usize) {
        v_u_bd_iter!((self.u_bds, v, u, bd) {
            Mpz::shr_assign(self.fcvec.fc_ref_mut(v, u, bd), other);
        }
        );
    }
}

impl Sqrt5Elt<HmfGen> {
    /// Return 0 q-expansion
    pub fn new(prec: usize) -> Self {
        let f = HmfGen::new(prec);
        let g = HmfGen::new(prec);
        Self { rt: f, ir: g }
    }

    /// Decrease prec to prec.
    pub fn decrease_prec(&mut self, prec: usize) {
        self.rt.decrease_prec(prec);
        self.ir.decrease_prec(prec);
    }

    pub fn u_bds(&self) -> &UBounds {
        &self.rt.u_bds
    }

    pub fn weight(&self) -> Weight {
        self.rt.weight
    }

    pub fn prec(&self) -> usize {
        self.rt.prec
    }
}

impl fmt::Display for Sqrt5Elt<HmfGen> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut vec = Vec::new();
        v_u_bd_iter!((self.u_bds(), v, u, bd) {
            let rt = self.rt.fcvec.fc_ref(v, u, bd);
            let ir = self.ir.fcvec.fc_ref(v, u, bd);
            if !rt.is_zero() | !ir.is_zero() {
                vec.push(format!("({}, {}): ({}, {})", u, v, rt, ir));
            }
        }
        );
        write!(f, "{}", vec.join("\n"))
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weight_add_mul() {
        let a: Weight = Some((1, 2));
        let b: Weight = Some((3, 2));
        assert_eq!(weight_mul(a, b), Some((4, 4)));
        assert_eq!(weight_mul(a, None), None);
        assert_eq!(weight_add(a, b), None);
        assert_eq!(weight_add(a, a), a);
    }
}

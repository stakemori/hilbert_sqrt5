use libc::{c_ulong, c_long};
use gmp::mpz::Mpz;
use std::fmt;
use std::ops::{AddAssign, SubAssign, ShlAssign, ShrAssign, MulAssign};

pub trait BigNumber {
    fn is_zero_g(&self) -> bool;
    fn from_ui_g(x: c_ulong) -> Self;
    fn from_si_g(x: c_long) -> Self;
    fn set_ui_g(&mut self, x: c_ulong);
    fn set_si_g(&mut self, x: c_long);
    fn set_g(&mut self, other: &Self);
    fn add_mut_g(&mut self, x: &Self, y: &Self);
    fn mul_mut_g(&mut self, x: &Self, y: &Self);
    fn sub_mut_g(&mut self, x: &Self, y: &Self);
    fn is_multiple_of_g(&self, x: &Self, tmpelt: &mut Self, tmp: &mut Mpz) -> bool;
    fn new_g() -> Self;
    fn negate_g(&mut self);
    fn set_divexact_g(&mut self, x: &Self, tmp: &mut Mpz);
    fn mul_assign_g(&mut self, other: &Self, tmp: &mut Mpz);
    fn addmul_mut_g(&mut self, x: &Self, y: &Self, tmp: &mut Mpz);
}

impl BigNumber for Mpz {
    fn is_zero_g(&self) -> bool {
        self.is_zero()
    }

    fn from_ui_g(x: c_ulong) -> Self {
        Mpz::from_ui(x)
    }

    fn from_si_g(x: c_long) -> Self {
        Mpz::from_si(x)
    }

    fn set_ui_g(&mut self, x: c_ulong) {
        self.set_ui(x);
    }

    fn set_si_g(&mut self, x: c_long) {
        self.set_si(x);
    }

    fn set_g(&mut self, other: &Mpz) {
        self.set(other)
    }

    fn add_mut_g(&mut self, x: &Mpz, y: &Mpz) {
        self.add_mut(x, y);
    }

    fn mul_mut_g(&mut self, x: &Mpz, y: &Mpz) {
        self.mul_mut(x, y);
    }

    fn sub_mut_g(&mut self, x: &Mpz, y: &Mpz) {
        self.sub_mut(x, y);
    }

    fn is_multiple_of_g(&self, x: &Mpz, _tmpelt: &mut Self, _tmp: &mut Mpz) -> bool {
        self.is_multiple_of(x)
    }

    fn new_g() -> Mpz {
        Mpz::new()
    }

    fn negate_g(&mut self) {
        self.negate()
    }

    fn set_divexact_g(&mut self, x: &Mpz, _tmp: &mut Mpz) {
        self.set_divexact(x);
    }

    fn mul_assign_g(&mut self, other: &Mpz, _tmp: &mut Mpz) {
        *self *= other;
    }

    fn addmul_mut_g(&mut self, x: &Mpz, y: &Mpz, _tmp: &mut Mpz) {
        self.addmul_mut(x, y);
    }
}

/// (rt + ir sqrt(5))/2
#[derive(Debug, Clone)]
pub struct Sqrt5Mpz {
    pub rt: Mpz,
    pub ir: Mpz,
}

impl Sqrt5Mpz {
    pub fn conj_mut(&mut self) {
        self.ir.negate();
    }

    pub fn norm(&self, res: &mut Mpz) {
        let &Self {
            rt: ref a,
            ir: ref b,
        } = self;
        res.mul_mut(b, b);
        *res *= -5 as c_long;
        res.addmul_mut(a, a);
        *res >>= 2;
    }

    pub fn from_sisi(rt: c_long, ir: c_long) -> Self {
        Self {rt: Mpz::from_si(rt), ir: Mpz::from_si(ir)}
    }
}

impl BigNumber for Sqrt5Mpz {
    fn is_zero_g(&self) -> bool {
        self.rt.is_zero() && self.ir.is_zero()
    }

    fn from_ui_g(x: c_ulong) -> Self {
        let mut a = Mpz::from_ui(x);
        a <<= 1;
        Self {
            rt: a,
            ir: Mpz::zero(),
        }
    }

    fn from_si_g(x: c_long) -> Self {
        let mut a = Mpz::from_si(x);
        a <<= 1;
        Self {
            rt: a,
            ir: Mpz::zero(),
        }
    }

    fn set_ui_g(&mut self, x: c_ulong) {
        self.rt.set_ui(x);
        self.rt <<= 1;
        self.ir.set_ui(0);
    }

    fn set_si_g(&mut self, x: c_long) {
        self.rt.set_si(x);
        self.rt <<= 1;
        self.ir.set_ui(0);
    }

    fn set_g(&mut self, other: &Self) {
        self.ir.set(&other.ir);
        self.rt.set(&other.rt);
    }

    fn add_mut_g(&mut self, x: &Self, y: &Self) {
        self.ir.add_mut(&x.ir, &y.ir);
        self.rt.add_mut(&x.rt, &y.rt);
    }

    fn sub_mut_g(&mut self, x: &Self, y: &Self) {
        self.ir.sub_mut(&x.ir, &y.ir);
        self.rt.sub_mut(&x.rt, &y.rt);
    }

    fn mul_mut_g(&mut self, x: &Self, y: &Self) {
        self.rt.mul_mut(&x.ir, &y.ir);
        self.rt *= 5 as c_ulong;
        self.rt.addmul_mut(&x.rt, &y.rt);
        self.ir.mul_mut(&x.rt, &y.ir);
        self.ir.addmul_mut(&x.ir, &y.rt);
        self.rt >>= 1;
        self.ir >>= 1;
    }

    fn new_g() -> Self {
        Self {
            rt: Mpz::new(),
            ir: Mpz::new(),
        }
    }

    fn negate_g(&mut self) {
        self.rt.negate();
        self.ir.negate();
    }

    fn set_divexact_g(&mut self, x: &Self, tmp: &mut Mpz) {
        // (x.conj() * y/y.norm()).conj()
        self.conj_mut();
        self.mul_assign_g(&x, tmp);
        x.norm(tmp);
        self.ir.set_divexact(tmp);
        self.rt.set_divexact(tmp);
        self.conj_mut();
    }

    fn is_multiple_of_g(&self, x: &Self, tmpelt: &mut Self, tmp: &mut Mpz) -> bool {
        tmpelt.set_g(self);
        tmpelt.conj_mut();
        tmpelt.mul_assign_g(&x, tmp);
        x.norm(tmp);
        if tmpelt.ir.is_multiple_of(&tmp) && tmpelt.rt.is_multiple_of(&tmp) {
            tmpelt.ir /= tmp as &Mpz;
            tmpelt.rt /= tmp as &Mpz;
            tmp.add_mut(&tmpelt.ir, &tmpelt.rt);
            // TODO: Use macro mpz_even_p.
            tmp.is_multiple_of_ui(2)
        } else {
            false
        }
    }

    fn mul_assign_g(&mut self, other: &Self, tmp: &mut Mpz) {
        let &mut Self {
            rt: ref mut a,
            ir: ref mut b,
        } = self;
        let &Self {
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

    fn addmul_mut_g(&mut self, x: &Self, y: &Self, tmp: &mut Mpz) {
        tmp.mul_mut(&x.ir, &y.ir);
        *tmp *= 5 as c_ulong;
        tmp.addmul_mut(&x.rt, &y.rt);
        *tmp >>= 1;
        self.rt += tmp as &Mpz;
        tmp.mul_mut(&x.rt, &y.ir);
        tmp.addmul_mut(&x.ir, &y.rt);
        *tmp >>= 1;
        self.ir += tmp as &Mpz;
    }
}

impl fmt::Display for Sqrt5Mpz {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.ir.is_zero() {
            write!(f, "{}", &self.rt >> 1)
        } else {
            write!(f, "({}, {})", self.rt, self.ir)
        }
    }
}

impl<'a> AddAssign<&'a Sqrt5Mpz> for Sqrt5Mpz {
    fn add_assign(&mut self, other: &Sqrt5Mpz) {
        self.rt += &other.rt;
        self.ir += &other.ir;
    }
}

impl<'a> SubAssign<&'a Sqrt5Mpz> for Sqrt5Mpz {
    fn sub_assign(&mut self, other: &Sqrt5Mpz) {
        self.rt -= &other.rt;
        self.ir -= &other.ir;
    }
}

impl MulAssign<c_long> for Sqrt5Mpz {
    fn mul_assign(&mut self, other: c_long) {
        self.ir *= other;
        self.rt *= other;
    }
}

impl MulAssign<c_ulong> for Sqrt5Mpz {
    fn mul_assign(&mut self, other: c_ulong) {
        self.ir *= other;
        self.rt *= other;
    }
}

impl ShlAssign<usize> for Sqrt5Mpz {
    fn shl_assign(&mut self, other: usize) {
        self.rt <<= other;
        self.ir <<= other;
    }
}

impl ShrAssign<usize> for Sqrt5Mpz {
    fn shr_assign(&mut self, other: usize) {
        self.rt >>= other;
        self.ir >>= other;
    }
}

impl PartialEq for Sqrt5Mpz {
    fn eq(&self, other: &Self) -> bool {
        self.rt == other.rt && self.ir == other.ir
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_divexact_g() {
        let a = Sqrt5Mpz{rt: Mpz::from_ui(1), ir: Mpz::from_ui(3)};
        let mut b = Sqrt5Mpz{rt: Mpz::from_ui(22), ir: Mpz::from_si(0)};
        let mut tmp = Mpz::new();
        b.set_divexact_g(&a, &mut tmp);
        assert_eq!(b, Sqrt5Mpz{rt: Mpz::from_si(-1), ir: Mpz::from_si(3)});
    }

    #[test]
    fn test_set_fun() {
        let a = Sqrt5Mpz::from_si_g(3);
        let mut b = Sqrt5Mpz::new_g();
        b.set_si_g(3);
        assert_eq!(a, b);
    }

    #[test]
    fn test_addmul_mut() {
        let mut tmp = Mpz::new();
        let a = Sqrt5Mpz::from_sisi(3, 5);
        let b = Sqrt5Mpz::from_sisi(7, 1);
        let mut res = Sqrt5Mpz::from_sisi(4, 6);
        res.addmul_mut_g(&a, &b, &mut tmp);
        assert_eq!(res, Sqrt5Mpz::from_sisi(27, 25));
    }

    #[test]
    fn test_mul_mut() {
        let mut res = Sqrt5Mpz::from_sisi(2, 4);
        let mut tmp = Mpz::new();
        let mut a = Sqrt5Mpz{rt: Mpz::from_si(5), ir: Mpz::from_si(1)};
        let b = Sqrt5Mpz{rt: Mpz::from_si(3), ir: Mpz::from_si(7)};
        res.mul_mut_g(&a, &b);
        a.mul_assign_g(&b, &mut tmp);
        assert_eq!(res, a);
    }

    #[test]
    fn test_is_multiple_of() {
        let ref mut tmpelt = Sqrt5Mpz::new_g();
        let ref mut tmp = Mpz::new();
        let a = Sqrt5Mpz::from_ui_g(10);
        let b = Sqrt5Mpz::from_ui_g(2);
        assert!(a.is_multiple_of_g(&b, tmpelt, tmp));
        let b = Sqrt5Mpz::from_ui_g(3);
        assert!(!a.is_multiple_of_g(&b, tmpelt, tmp));
        let a = Sqrt5Mpz{rt: Mpz::from_si(2), ir: Mpz::from_si(4)};
        let b = Sqrt5Mpz::from_ui_g(2);
        assert!(!a.is_multiple_of_g(&b, tmpelt, tmp));
        let a = Sqrt5Mpz::from_ui_g(33);
        let b = Sqrt5Mpz{rt: Mpz::from_si(-1), ir: Mpz::from_si(3)};
        assert!(a.is_multiple_of_g(&b, tmpelt, tmp));
    }
}

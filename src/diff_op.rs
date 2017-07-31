extern crate libc;

use misc::Sqrt5Elt;
use gmp::mpz::Mpz;
use self::libc::{c_ulong, c_long};

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

    pub fn pow_mut(&mut self, expt: usize, tmp_elt: &mut Self, tmp: &mut Mpz) {
        tmp_elt.set(self);
        self.set_one();
        let s = format!("{:b}", expt);
        let bts = s.into_bytes();
        let strs = bts.iter().rev().map(|&i| i as char).collect::<Vec<char>>();
        for &c in strs.iter() {
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

#[cfg(test)]
mod tests {
    use super::*;

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

        a.pow_mut(10, &mut tmp_elt, &mut tmp);
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
}

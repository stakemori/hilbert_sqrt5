extern crate libc;

use misc::Sqrt5Elt;
use gmp::mpz::Mpz;
use self::libc::{c_ulong, c_long};

impl Sqrt5Elt<Mpz> {
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
        a.addmul_mut(b, d);
        *a >>= 1;
    }

    pub fn norm_mut(self: &Self, res: &mut Mpz) {
        let &Sqrt5Elt {
            rt: ref a,
            ir: ref b,
        } = self;
        res.mul_mut(b, b);
        *res *= -5 as c_long;
        res.addmul_mut(a, a);
        *res >>= 2;
    }
}

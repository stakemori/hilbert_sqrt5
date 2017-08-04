use libc::{c_ulong, c_long};
use gmp::mpz::Mpz;

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
    fn is_multiple_of_g(&self, x: &Self) -> bool;
    fn new_g() -> Self;
    fn negate_g(&mut self);
    fn set_divexact_g(&mut self, x: &Self);
}

impl BigNumber for Mpz {
    fn is_zero_g(&self) -> bool {
        self.is_zero()
    }

    fn from_ui_g(x: c_ulong)  -> Self {
        Mpz::from_ui(x)
    }

    fn from_si_g(x: c_long)  -> Self {
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

    fn is_multiple_of_g(&self, x: &Mpz) -> bool {
        self.is_multiple_of(x)
    }

    fn new_g() -> Mpz {
        Mpz::new()
    }

    fn negate_g(&mut self) {
        self.negate()
    }

    fn set_divexact_g(&mut self, x: &Mpz) {
        self.set_divexact(x);
    }
}


#![crate_name = "hilbert_sqrt5"]
#![warn(deprecated)]
extern crate gmp;
extern crate libc;
#[macro_use]
pub mod elements;
pub mod theta_chars;
pub mod eisenstein;
pub mod diff_op;
pub mod misc;
mod fcvec;
mod bignum;

#[cfg(test)]
mod test;

#![crate_name = "hilbert_sqrt5"]
#![warn(deprecated)]
extern crate gmp;

#[macro_use]
pub mod elements;
pub mod theta_chars;
pub mod eisenstein;
pub mod diff_op;
mod misc;

#[cfg(test)]
mod test;

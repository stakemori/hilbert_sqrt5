#![feature(test)]

extern crate test;
extern crate hilbert_sqrt5;

use hilbert_sqrt5::misc::prime_sieve;

use test::Bencher;

#[bench]
fn prime_sieve_bench(b: &mut Bencher) {
    b.iter(|| { prime_sieve(10000); })
}

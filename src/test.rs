use std::time::Instant;
use theta_chars::theta;
use eisenstein::eisenstein_series;
use misc::prime_sieve;
use gmp::mpz::Mpz;


// Taken from http://qiita.com/pseudo_foxkeh/items/5d5226e3ffa27631e80d
macro_rules! measure_time {
  ( $x:expr) => {
    {
      let start = Instant::now();
      let result = $x;
      let end = start.elapsed();
      println!("{}.{:03} seconds passed", end.as_secs(), end.subsec_nanos() / 1_000_000);
      result
    }
  };
}
mod theta {
    use super::*;

    #[test]
    fn theta_fun() {
        measure_time!(theta(10));
    }
}

mod eisen {
    use super::*;

    #[test]
    fn test_eisen() {
        measure_time!(eisenstein_series(4, &Mpz::one(), &Mpz::from_ui(240), 30));
    }
}

mod misc {
    use super::*;

    #[test]
    fn prime_sieve_test() {
        let n = 1009;
        let v = prime_sieve(n);
        assert_eq!(v.len(), 169);
    }
}

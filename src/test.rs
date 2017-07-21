use std::time::Instant;
use theta_chars::theta;

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

#[ignore]
#[test]
fn theta_fun() {
    measure_time!(theta(10));
}

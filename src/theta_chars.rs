/// Return vec of (x, y) s.t. (x^2 + 2uv + 5v^2)/4 <= r and x equiv a, y equiv b
/// mod 2, (x-a)/2 equiv (x-b)/2 mod 2.
pub fn points_in_ellipse(a: i64, b: i64, r: i64) -> Vec<(i64, i64)> {
    let mut vec: Vec<(i64, i64)> = Vec::new();
    let n1 = (r as f64).sqrt().floor();
    let dab = a - b;
    debug_assert!(n1 < ::std::i64::MAX as f64);
    let n1 = n1 as i64;
    for y in -n1..(n1 + 1) {
        let n2 = (2_f64 * ((r - y * y) as f64).sqrt()).floor();
        debug_assert!(n2 < ::std::i64::MAX as f64);
        let n2 = n2 as i64;
        for x in (-n2 - y)..(n2 - y + 1) {
            debug_assert!(x * x + 2 * x * y + 5 * y * y <= 4 * r);
            if (x - a) & 1 == 0 && (y - b) & 1 == 0 && (x - y - dab) & 0b11 == 0 {
                vec.push((x, y))
            }
        }
    }
    vec
}

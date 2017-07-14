/// Return vec of (u, v) s.t. (u^2 + 2uv + 5v^2)/4 <= r and u equiv a, v equiv b
/// mod 2.
pub fn points_in_ellipse(a: i64, b: i64, r: i64) -> Vec<(i64, i64)> {
    let mut vec: Vec<(i64, i64)> = Vec::new();
    let n1 = (r as f64).sqrt().floor() as i64;
    for y in -n1..(n1 + 1) {
        let n2 = (2_f64 * ((r - y * y) as f64).sqrt()).floor() as i64;
        for x in (-n2 - y)..(n2 - y + 1) {
            debug_assert!(x * x + 2 * x * y + 5 * y * y <= 4 * r);
            if (x - a) & 1 == 0 && (y - b) & 1 == 0 {
                vec.push((x, y))
            }
        }
    }
    vec
}

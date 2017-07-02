/// Return vec of (u, v) s.t. (u^2 + 5v^2)/2 <= r^2 and u equiv a, v equiv b mod
/// 2.
pub fn points_in_ellipse(a: i64, b: i64, r: i64) -> Vec<(i64, i64)> {
    let mut vec: Vec<(i64, i64)> = Vec::new();
    let r1 = (2_f64.sqrt() * r as f64 / 5_f64.sqrt()).floor() as i64;
    let t_r_sq = 2 * r * r;
    for v in -r1..(r1 + 1) {
        let f_v_sq = 5 * v * v;
        let r2 = ((t_r_sq - f_v_sq) as f64).sqrt().floor() as i64;
        for u in -r2..(r2 + 1) {
            if (u - a) & 1 == 0 && (v - b) & 1 == 0 && u * u + f_v_sq <= t_r_sq {
                vec.push((u, v))
            }
        }
    }
    vec
}

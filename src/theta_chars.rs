/// Implemention of Theta in Satz 2, Gundlach, Die Bestimmung der Functionen zur
/// Hilbertschen Modulargruppe des Zahlkörpers Q(sqrt(5)).

use elements::HmfGen;
use gmp::mpz::Mpz;
use std::ops::{AddAssign, SubAssign};

/// Return vec of (x, y) s.t. (x^2 + 2xy + 5y^2)/4 <= r and x equiv a, y equiv b
/// mod 2, (x-a)/2 equiv (x-b)/2 mod 2.
fn points_in_ellipse(a: i64, b: i64, r: i64) -> Vec<(i64, i64)> {
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

/// corresponds to (rt + ir sqrt(5))/2.
#[derive(Debug)]
struct Sqrt5Elt<T> {
    rt: T,
    ir: T,
}

fn theta_char(alpha: &Sqrt5Elt<i64>, beta: &Sqrt5Elt<i64>, prec: usize) -> HmfGen {
    let mut res = HmfGen::new(prec);
    for &(x, y) in points_in_ellipse(alpha.rt, alpha.ir, prec as i64).iter() {
        let nu_1 = (x - alpha.rt) >> 1;
        let nu_2 = (y - alpha.ir) >> 1;
        let s = (x * x + 5 * y * y) >> 1;
        let t = x * y;
        let u = (s + 5 * t) >> 1;
        let v = (s + t) >> 1;
        let trc = (nu_1 * beta.ir + nu_2 * beta.rt) >> 1;
        let v = v as usize;
        let bd = res.u_bds.vec[v] as i64;
        if trc & 1 == 0 {
            Mpz::add_assign(res.fcvec.fc_ref_mut(v, u, bd), 1);
        } else {
            Mpz::sub_assign(res.fcvec.fc_ref_mut(v, u, bd), 1);
        }
    }
    res
}

fn theta_0_eps(prec: usize) -> HmfGen {
    let one = &Sqrt5Elt { rt: 2, ir: 0 };
    let zero = &Sqrt5Elt { rt: 0, ir: 0 };
    let eps = &Sqrt5Elt { rt: 1, ir: 1 };
    let eps_star = &Sqrt5Elt { rt: 1, ir: -1 };
    let ary = [
        (one, zero),
        (eps, zero),
        (eps_star, zero),
        (zero, zero),
        (eps_star, eps),
        (eps, eps_star),
        (zero, one),
        (one, one),
        (zero, eps_star),
        (zero, eps),
    ];
    let prec8 = prec * 8;
    let mut res = theta_char(ary[0].0, ary[0].1, prec8);
    for &(alpha, beta) in (&ary[1..10]).iter() {
        res *= &theta_char(alpha, beta, prec8);
    }
    res
}

/// Return Theta in Satz 2 in Gundlach's paper.
pub fn theta(prec: usize) -> HmfGen {
    let f = theta_0_eps(prec);
    let mut res = HmfGen::new(prec);
    v_u_bd_iter!((res.u_bds, v, u, bd) {
        let v1 = v << 3;
        let u1 = u << 3;
        let bd1 = f.u_bds.vec[v1] as i64;
        let a = f.fcvec.fc_ref(v1, u1, bd1);
        if !a.is_zero() {
            res.fcvec.fc_ref_mut(v, u, bd).set(a);
        }
    });
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[ignore]
    #[test]
    fn test_points_in_ellipse() {
        let r = 10000;
        assert_eq!(points_in_ellipse(0, 0, r).len(), 7845);
        assert_eq!(points_in_ellipse(0, 1, r).len(), 7854);
        assert_eq!(points_in_ellipse(1, 0, r).len(), 7842);
        assert_eq!(points_in_ellipse(1, 1, r).len(), 7860);
    }

    #[test]
    fn test_theta_0_eps() {
        let f = theta_0_eps(5);
        v_u_bd_iter!((f.u_bds, v, u, bd) {
            let a = f.fcvec.fc_ref(v, u, bd);
            if !a.is_zero() {
                assert!(u & 0b111 == 0 && v & 0b111 == 0);
            }
        });
    }
}

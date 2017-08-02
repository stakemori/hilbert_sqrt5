/// Implemention of Theta in Satz 2, Gundlach, Die Bestimmung der Functionen zur
/// Hilbertschen Modulargruppe des Zahlkörpers Q(sqrt(5)).
use elements::HmfGen;
use eisenstein::eisenstein_series;
use misc::{Sqrt5Elt, PowGen};
use gmp::mpz::Mpz;
use std::ops::{AddAssign, SubAssign};
use fcvec;

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
            if is_even!(x - a) && is_even!(y - b) && (x - y - dab) & 0b11 == 0 {
                vec.push((x, y))
            }
        }
    }
    vec
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
        if is_even!(trc) {
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

fn theta_squared(prec: usize) -> HmfGen {
    let e2 = eisenstein_series(2, prec);
    let e6 = eisenstein_series(6, prec);
    let mut res = eisenstein_series(10, prec);
    let mut tmp = HmfGen::new(prec);
    tmp.pow_mut(&e2, 2);
    let mut f10_2 = &tmp * &e6;
    tmp.square();
    tmp *= &e2;
    let mut f10_1 = tmp.clone();
    f10_2 *= &Mpz::from_si(-11465);
    f10_1 *= &Mpz::from_ui(355404);
    res += &f10_2;
    res += &f10_1;
    res /= &Mpz::from_ui(5443200000);
    res
}


/// Return normalized cusp form of weight 5 that is propotional to the return
/// value of theta(prec).
pub fn g5_normalized(prec: usize) -> HmfGen {
    let prec = prec + 1;
    let g10 = theta_squared(prec);
    let mut f10 = HmfGen::new(prec);
    divide_by_squared(&mut f10, &g10);

    let mut tmp = HmfGen::new(prec);
    let mut res = HmfGen::new(prec - 1);
    {
        let bd1 = res.u_bds.vec[1] as i64;
        res.fcvec.fc_ref_mut(1, 1, bd1 as i64).set_ui(1);
        res.fcvec.fc_ref_mut(1, -1, bd1 as i64).set_si(-1);
    }

    let u_bds = g10.u_bds;
    for v in 3..(prec + 1) {
        let v_d = if is_even!(v) { v >> 1 } else { (v >> 1) + 1 };
        if is_even!(v) {
            fcvec::mul_mut(
                &mut tmp.fcvec.vec[v],
                &f10.fcvec.vec[v_d + 1],
                &f10.fcvec.vec[v_d + 1],
                v_d,
                v_d,
                u_bds.vec[v_d + 1],
                u_bds.vec[v_d + 1],
                u_bds.vec[v],
                &u_bds,
                1,
                1,
                0,
            );
            fcvec::sub_assign(&mut f10.fcvec.vec[v], &tmp.fcvec.vec[v], v, &u_bds);
        }
        fcvec::shr_assign(&mut f10.fcvec.vec[v], v, &u_bds, 1);
        for i in 2..v_d {
            fcvec::mul_mut(
                &mut tmp.fcvec.vec[v],
                &f10.fcvec.vec[i + 1],
                &f10.fcvec.vec[v - i + 1],
                i,
                v - i,
                u_bds.vec[i + 1],
                u_bds.vec[v - i + 1],
                u_bds.vec[v],
                &u_bds,
                1,
                1,
                0,
            );
            fcvec::sub_assign(&mut f10.fcvec.vec[v], &tmp.fcvec.vec[v], v, &u_bds);
        }
    }

    let mut tmp_z = Mpz::new();
    for v in 2..prec {
        let bd = u_bds.vec[v] as i64;
        let bd_1 = u_bds.vec[v + 1] as i64;
        let v_i = v as i64;
        if is_even!(v_i + bd) {
            res.fcvec.fc_ref_mut(v, bd, bd).set(f10.fcvec.fc_ref(
                v + 1,
                bd - 1,
                bd_1,
            ));
        }
        for u in (1..bd).filter(|i| is_even!(i + v_i)) {
            tmp_z.set(f10.fcvec.fc_ref(v + 1, u - 1, bd_1));
            tmp_z -= f10.fcvec.fc_ref(v + 1, u + 1, bd_1);
            res.fcvec.fc_ref_mut(v, u, bd).set(&tmp_z);
        }
        for u in (1..(bd + 1)).filter(|i| is_even!(i + v_i)) {
            tmp_z.set(res.fcvec.fc_ref(v, u, bd));
            tmp_z.negate();
            res.fcvec.fc_ref_mut(v, -u, bd).set(&tmp_z);
        }
    }
    res.weight = Some((5, 5));
    res
}

/// set f = g / (q_1 - q_1^(-1))^2, where q_1 = e(u).
fn divide_by_squared(f: &mut HmfGen, g: &HmfGen) {
    let prec = g.prec;
    for v in 2..(prec + 1) {
        let bd = g.u_bds.vec[v] - 2;
        let v_mod2 = if is_even!(v) { 0 as usize } else { 1 as usize };
        let bd_i = bd as i64;
        let mut tmp = Mpz::new();
        if is_even!(bd + v_mod2) {
            f.fcvec.fc_ref_mut(v, bd_i, bd_i + 2).set(g.fcvec.fc_ref(
                v,
                bd_i + 2,
                bd_i + 2,
            ));
        }

        let bd_i = bd as i64;
        for u in (0..bd).filter(|&x| is_even!(x + v_mod2 + bd - 1)) {
            let u_i = bd_i - 1 - (u as i64);
            tmp.set(g.fcvec.fc_ref(v, u_i + 2, bd_i + 2));
            if u_i + 4 <= bd_i {
                tmp -= f.fcvec.fc_ref(v, u_i + 4, bd_i + 2);
            }
            if u_i + 2 <= bd_i {
                tmp += f.fcvec.fc_ref(v, u_i + 2, bd_i + 2);
                tmp += f.fcvec.fc_ref(v, u_i + 2, bd_i + 2);
            }
            f.fcvec.fc_ref_mut(v, u_i, bd_i + 2).set(&tmp);
        }

        for u in (1..(bd + 1)).filter(|&x| is_even!(x + v_mod2)) {
            let u_i = u as i64;
            tmp.set(f.fcvec.fc_ref(v, u_i, bd_i + 2));
            f.fcvec.fc_ref_mut(v, -u_i, bd_i + 2).set(&tmp);
        }
    }
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

    #[ignore]
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

    fn print_vth_cf(f: &HmfGen, v: usize) {
        let bd = f.u_bds.vec[v] as i64;
        let v_i = v as i64;
        let mut res = Vec::new();
        for u in (-bd..(bd + 1)).filter(|u| is_even!(v_i + u)) {
            let a = f.fcvec.fc_ref(v, u, bd);
            if !a.is_zero() {
                res.push(format!("({}) * q1**({})", a, u));
            }
        }
        println!("{}", res.join(" + "));
    }


    #[ignore]
    #[test]
    fn test_divide() {
        let mut f = HmfGen::new(10);
        let g = theta_squared(10);
        divide_by_squared(&mut f, &g);
        for v in 2..11 {
            println!("{}", v);
            print_vth_cf(&g, v);
            print_vth_cf(&f, v);
        }
    }

    #[test]
    fn test_mul_mut() {
        let mut f = HmfGen::new(10);
        let g = eisenstein_series(2, 10);
        let h = eisenstein_series(4, 10);
        let v_g = 2;
        let v_h = 4;
        fcvec::mul_mut(
            &mut f.fcvec.vec[v_g + v_h],
            &g.fcvec.vec[v_g],
            &h.fcvec.vec[v_h],
            v_g,
            v_h,
            f.u_bds.vec[v_g],
            f.u_bds.vec[v_h],
            f.u_bds.vec[v_g + v_h],
            &f.u_bds,
            0,
            0,
            0,
        );
        print_vth_cf(&g, v_g);
        print_vth_cf(&h, v_h);
        print_vth_cf(&f, v_g + v_h);
    }
}

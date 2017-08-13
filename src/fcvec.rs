use gmp::mpz::Mpz;
use elements::UBounds;
use std::ops::{SubAssign, ShrAssign, ShlAssign};
use bignum::BigNumber;


#[allow(dead_code)]
pub fn shl_assign<T>(f_vec: &mut Vec<T>, v: usize, u_bds: &UBounds, a: usize)
where
    T: ShlAssign<usize>,
{
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd + i] <<= a;
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd - i] <<= a;
    }
}

pub fn shr_assign<T>(f_vec: &mut Vec<T>, v: usize, u_bds: &UBounds, a: usize)
where
    T: ShrAssign<usize>,
{
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd + i] >>= a;
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd - i] >>= a;
    }
}

pub fn sub_assign<T>(f_vec: &mut Vec<T>, g_vec: &Vec<T>, v: usize, u_bds: &UBounds)
where
    for<'a> T: SubAssign<&'a T>,
{
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        T::sub_assign(&mut f_vec[bd + i], &g_vec[bd + i]);
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        let idx = bd - i;
        T::sub_assign(&mut f_vec[idx], &g_vec[idx])
    }
}

/// set v_g + v_h coefficient (Laurant polynomial of e(u)) to the product of
/// g[v_g], h[v_h].
/// g_vec[gap_g + i] is defined only when i.abs() <= bd_g.
pub fn mul_mut<T>(
    f_vec: &mut Vec<T>,
    g_vec: &Vec<T>,
    h_vec: &Vec<T>,
    v_g: usize,
    v_h: usize,
    gap_g: usize,
    gap_h: usize,
    gap_gh: usize,
    u_bds: &UBounds,
    parity_g: usize,
    parity_h: usize,
    parity_gh: usize,
) where
    T: BigNumber,
{
    let bd_g = u_bds.vec[v_g];
    let bd_h = u_bds.vec[v_h];
    let bd_gh = u_bds.vec[v_g + v_h];

    for i in (0..(bd_gh + 1)).filter(|&x| is_even!(v_g + v_h + x + parity_gh)) {
        f_vec[(gap_gh + i) as usize].set_ui_g(0);
        f_vec[(gap_gh - i) as usize].set_ui_g(0);
    }
    let mut tmp = Mpz::new();
    // naive implementation of polynomial multiplication
    // i -> i - bd_g
    for i in (0..(2 * bd_g + 1)).filter(|&x| is_even!(v_g + x + bd_g + parity_g)) {
        for j in (0..(2 * bd_h + 1)).filter(|&x| is_even!(v_h + x + bd_h + parity_h)) {
            f_vec[gap_gh + i + j - bd_g - bd_h].addmul_mut_g(
                &g_vec[gap_g + i - bd_g],
                &h_vec[gap_h + j - bd_h],
                &mut tmp,
            );
        }
    }
}

/// Assuming f_vec is divisible by g_vec in integral coefficients, set h_vec to
/// f_vec/g_vec.
pub fn div_mut<T>(
    f_vec: &Vec<T>,
    g_vec: &Vec<T>,
    h_vec: &mut Vec<T>,
    v_g: usize,
    v_h: usize,
    gap_g: usize,
    gap_h: usize,
    gap_gh: usize,
    u_bds: &UBounds,
) where
    T: BigNumber,
    for<'a> T: SubAssign<&'a T>,
{
    let gap_gh = gap_gh as i64;
    let gap_g = gap_g as i64;
    let gap_h = gap_h as i64;
    let bd_g = u_bds.vec[v_g] as i64;
    let bd_h = u_bds.vec[v_h] as i64;
    for i in -bd_h..(bd_h + 1) {
        h_vec[(gap_h + i) as usize].set_ui_g(0);
    }
    let mut tmp = Mpz::new();
    // initial index of g
    let n = ((-bd_g)..(bd_g + 1))
        .rev()
        .filter(|&i| !g_vec[(gap_g + i) as usize].is_zero_g())
        .next()
        .unwrap();

    let mut tmp_elt = T::new_g();
    h_vec[(gap_h + bd_h) as usize].set_g(&f_vec[(gap_gh + n + bd_h) as usize]);
    h_vec[(gap_h + bd_h) as usize].set_divexact_g(&g_vec[(gap_g + n) as usize], &mut tmp);

    for m in ((-bd_h)..bd_h).rev() {
        tmp_elt.set_g(&f_vec[(gap_gh + n + m) as usize]);
        for i in ((n + m - bd_h)..n).filter(|&x| -bd_g <= x) {
            tmp_elt.submul_mut_g(
                &g_vec[(gap_g + i) as usize],
                &h_vec[(gap_h + n + m - i) as usize],
                &mut tmp,
            );

        }
        h_vec[(gap_h + m) as usize].set_g(&tmp_elt);
        h_vec[(gap_h + m) as usize].set_divexact_g(&g_vec[(gap_g + n) as usize], &mut tmp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use theta_chars::g5_normalized;
    use elements::HmfGen;

    // fn print_vth_cf(f: &HmfGen<Mpz>, v: usize) {
    //     let bd = f.u_bds.vec[v] as i64;
    //     let v_i = v as i64;
    //     let mut res = Vec::new();
    //     for u in (-bd..(bd + 1)).filter(|u| is_even!(v_i + u)) {
    //         let a = f.fcvec.fc_ref(v, u, bd);
    //         if !a.is_zero() {
    //             res.push(format!("({}) * q1**({})", a, u));
    //         }
    //     }
    //     println!("{}", res.join(" + "));
    // }

    #[test]
    fn test_divide() {
        let prec = 10;
        let a = g5_normalized(prec);
        let mut b = HmfGen::new(prec);
        b.pow_mut(&a, 3);
        let c = &a * &b;
        let mut h = HmfGen::<Mpz>::new(prec);
        let u_bds = b.u_bds;
        div_mut(
            &c.fcvec.vec[4],
            &a.fcvec.vec[1],
            &mut h.fcvec.vec[3],
            1,
            3,
            u_bds.vec[1],
            u_bds.vec[3],
            u_bds.vec[4],
            &u_bds,
        );
        let v = vec![0, 0, 0, -1, 0, 3, 0, -3, 0, 1, 0, 0, 0];
        let v: Vec<_> = v.iter().map(|&x| Mpz::from_si(x)).collect();
        assert_eq!(h.fcvec.vec[3], v);
    }
}

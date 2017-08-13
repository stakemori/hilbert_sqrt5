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

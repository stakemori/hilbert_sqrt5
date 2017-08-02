use gmp::mpz::Mpz;
use elements::UBounds;
use std::ops::SubAssign;

#[allow(dead_code)]
pub fn shl_assign(f_vec: &mut Vec<Mpz>, v: usize, u_bds: &UBounds, a: usize) {
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd + i] <<= a;
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd - i] <<= a;
    }
}

pub fn shr_assign(f_vec: &mut Vec<Mpz>, v: usize, u_bds: &UBounds, a: usize) {
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd + i] >>= a;
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        f_vec[bd - i] >>= a;
    }
}

pub fn sub_assign(f_vec: &mut Vec<Mpz>, g_vec: &Vec<Mpz>, v: usize, u_bds: &UBounds) {
    let bd = u_bds.vec[v];
    for i in (0..(bd + 1)).filter(|x| is_even!(x + v)) {
        Mpz::sub_assign(&mut f_vec[bd + i], &g_vec[bd + i]);
    }
    for i in (1..(bd + 1)).filter(|x| is_even!(x + v)) {
        let idx = bd - i;
        Mpz::sub_assign(&mut f_vec[idx], &g_vec[idx])
    }
}

/// set v_g + v_h coefficient (Laurant polynomial of e(u)) to the product of
/// g[v_g], h[v_h].
pub fn mul_mut(
    f_vec: &mut Vec<Mpz>,
    g_vec: &Vec<Mpz>,
    h_vec: &Vec<Mpz>,
    v_g: usize,
    v_h: usize,
    gap_g: usize,
    gap_h: usize,
    gap_gh: usize,
    u_bds: &UBounds,
    parity_g: usize,
    parity_h: usize,
    parity_gh: usize,
) {
    let bd_g = u_bds.vec[v_g];
    let bd_h = u_bds.vec[v_h];
    let bd_gh = u_bds.vec[v_g + v_h];

    for i in (0..(bd_gh + 1)).filter(|&x| is_even!(v_g + v_h + x + parity_gh)) {
        f_vec[(i + gap_gh) as usize].set_ui(0);
        f_vec[(gap_gh - i) as usize].set_ui(0);
    }

    // naive implementation of polynomial multiplication
    // i -> i - bd_g
    for i in (0..(2 * bd_g + 1)).filter(|&x| is_even!(v_g + x + bd_g + parity_g)) {
        for j in (0..(2 * bd_h + 1)).filter(|&x| is_even!(v_h + x + bd_h + parity_h)) {
            f_vec[gap_gh + i + j - bd_g - bd_h].addmul_mut(
                &g_vec[i + gap_g - bd_g],
                &h_vec[j + gap_h - bd_h],
            );
        }
    }
}

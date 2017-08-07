use gmp::mpz::Mpz;
use eisenstein::eisenstein_series;
use theta_chars::g5_normalized;
use elements::HmfGen;
use serde_pickle;
use std::fs::File;
use std::io::Write;
use std::io::Read;

// R = C[g2, g5, g6]

/// Return a vector of (a, b, c) s.t. 2*a + 5*b + 6*c = k.
fn tpls_of_wt(k: usize) -> Vec<(usize, usize, usize)> {
    let mut res = Vec::new();
    let c_max = k / 6;
    for c in 0..(c_max + 1) {
        let b_max = (k - 6 * c) / 5;
        for b in 0..(b_max + 1) {
            let rem = k - (6 * c + 5 * b);
            if is_even!(rem) {
                res.push((rem >> 1, b, c));
            }
        }
    }
    res
}

fn monom_g2_g5_g6(prec: usize, expts: (usize, usize, usize)) -> HmfGen<Mpz> {
    let mut tmp = HmfGen::new(prec);
    let mut res = HmfGen::new(prec);
    let g2 = eisenstein_series(2, prec);
    let g5 = g5_normalized(prec);
    let g6 = eisenstein_series(6, prec);
    let (e2, e5, e6) = expts;
    res.pow_mut(&g2, e2);
    tmp.pow_mut(&g6, e6);
    if e6 > 0 {
        res *= &tmp;
    }
    tmp.pow_mut(&g5, e5);
    if e5 > 0 {
        res *= &tmp;
    }
    res
}

fn save_as_pickle_z(vec: &Vec<Mpz>, path_name: &str) {
    let vec: Vec<String> = vec.iter().map(|x| x.to_str_radix(10)).collect();
    let v = serde_pickle::to_vec(&vec, false).unwrap();
    let mut buffer = File::create(path_name).unwrap();
    buffer.write(&v).unwrap();
}

fn load_pickle_z(path_name: &str) -> Vec<Mpz> {
    let file = File::open(path_name).unwrap();
    let buf: Vec<u8> = file.bytes().map(|x| x.unwrap()).collect();
    let v: Vec<String> = serde_pickle::from_slice(&buf).unwrap();
    v.iter().map(|x| Mpz::from_str_radix(x, 10).unwrap()).collect()
}

pub fn monoms_of_g2_g5_g6(k: usize, prec: usize) -> Vec<HmfGen<Mpz>> {
    tpls_of_wt(k)
        .iter()
        .map(|&x| monom_g2_g5_g6(prec, x))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tpls_of_wt() {
        assert_eq!(tpls_of_wt(30).len(), 13);
        assert_eq!(tpls_of_wt(100).len(), 99);
    }

    #[test]
    fn test_pickle() {
        let v: Vec<Mpz> = (0..10000).map(|x| Mpz::from_ui(x)).collect();
        save_as_pickle_z(&v, "/home/sho/foo.sobj");
        let w = load_pickle_z("/home/sho/foo.sobj");
        assert_eq!(v, w);
    }

}

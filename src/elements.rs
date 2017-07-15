extern crate gmp;

use self::gmp::mpz::Mpz;
use std;

/// struct for hilbert modualr form over Q(sqrt(5))
/// this corresponds finite sum of the q-expansion of the form
/// Σ a(u, v) exp(2piTr 1/sqrt(5) (u + v * sqrt(5))/2)
/// where v <= prec and ε = (1 + sqrt(5))/2.
/// a(u, v) = fc[v][a], where a = u + u_bds[v]
pub struct HmfGen {
    pub prec: usize,
    pub fcvec: FcVec,
    // vth element of u_bds.vec is (sqrt(5) * v).floor()
    u_bds: UBounds,
}

pub struct FcVec {
    pub vec: Vec<Vec<Mpz>>,
}

impl FcVec {
    fn fc_inner(&self, v: usize, u: i64, bd: i64) -> &Mpz {
        debug_assert!(u + bd >= 0);
        &self.vec[v][(u + bd) as usize]
    }

    fn fc_inner_mut(&mut self, v: usize, u: i64, bd: i64) -> &mut Mpz {
        debug_assert!(u + bd >= 0);
        &mut self.vec[v][(u + bd) as usize]
    }

    fn new(u_bds: &UBounds) -> FcVec {
        let vec = u_bds
            .vec
            .iter()
            .map(|&bd| (0..(2 * bd + 1)).map(|_| Mpz::new()).collect())
            .collect();
        FcVec { vec: vec }
    }
}

struct UBounds {
    vec: Vec<usize>,
}

impl UBounds {
    fn new(prec: usize) -> UBounds {
        assert!(5 * prec * prec < std::usize::MAX);
        let mut u_bds = Vec::new();
        let sqrt5 = 5_f64.sqrt();
        for v in 0..prec {
            u_bds.push((sqrt5 * v as f64).floor() as usize);
        }
        UBounds { vec: u_bds }
    }
}

impl HmfGen {
    pub fn new(prec: usize) -> HmfGen {
        let u_bds = UBounds::new(prec);
        let fcvec = FcVec::new(&u_bds);
        HmfGen {
            prec: prec,
            fcvec: fcvec,
            u_bds: u_bds,
        }
    }

    /// set self = f1 + f2
    pub fn add_mut(&mut self, f1: &HmfGen, f2: &HmfGen) {
        for (v, &bd) in self.u_bds.vec.iter().enumerate() {
            let bd = bd as i64;
            for u in -bd..(bd + 1) {
                Mpz::add_mut(
                    self.fcvec.fc_inner_mut(v, u, bd),
                    f1.fcvec.fc_inner(v, u, bd),
                    f2.fcvec.fc_inner(v, u, bd),
                );
            }
        }
    }

    /// set self = f1 * f2
    pub fn mul_mut(&mut self, f1: &HmfGen, f2: &HmfGen) {
        let mut tmp = Mpz::from_ui(0);
        for (v, &bd) in self.u_bds.vec.iter().enumerate() {
            let bd = bd as i64;
            for u in -bd..(bd + 1) {
                tmp.set_ui(0);

                for v1 in 0..v {
                    let bd1 = self.u_bds.vec[v1] as i64;
                    for u1 in -bd1..(bd1 + 1) {
                        let u2 = u - u1;
                        let v2 = v - v1;
                        let bd2 = self.u_bds.vec[v2] as i64;
                        let u2abs = u2.abs() as usize;
                        if u2abs * u2abs <= 5 * v2 * v2 {
                            tmp.add_mut(
                                f1.fcvec.fc_inner(v1, u1, bd1),
                                f2.fcvec.fc_inner(v2, u2, bd2),
                            );
                        }
                    }
                }

                self.fcvec.fc_inner_mut(v, u, bd).set(&tmp);
            }
        }
    }
}

/// Implemention of Eisenstein series Gundlach, Die Bestimmung der Functionen zur
/// Hilbertschen Modulargruppe des Zahlkörpers Q(sqrt(5)).

use elements::{HmfGen, UBounds};
use gmp::mpz::Mpz;
use misc::{Sqrt5Elt, prime_sieve};
use std::ops::MulAssign;


pub fn eisenstein_series(k: u64, l_val_num: &Mpz, l_val_denom: &Mpz, prec: usize) -> HmfGen {
    assert!((prec as f64) < (::std::i64::MAX as f64) * 5_f64.sqrt());
    let u_bds = UBounds::new(prec);
    let mut res = HmfGen::new(prec);

    // Set constant term
    res.fcvec.fc_ref_mut(0, 0, 0).set(&l_val_num);

    let expt = k - 1;

    for (v, &bd) in u_bds.vec.iter().enumerate().skip(1) {
        let bd = bd as i64;
        let v_i = v as i64;
        for u in u_iter!(v_i, bd) {
            res.fcvec.fc_ref_mut(v, u, bd).set_ui(1);
        }
    }

    split_primes_factor(prec, expt, &mut res, &u_bds);
    innert_primes_factor(prec, expt, &mut res, &u_bds);
    sqrt5_factor(prec, expt, &mut res, &u_bds);

    for (v, &bd) in u_bds.vec.iter().enumerate().skip(1) {
        let bd = bd as i64;
        let v_i = v as i64;
        for u in u_iter!(v_i, bd) {
            Mpz::mul_assign(res.fcvec.fc_ref_mut(v, u, bd), l_val_denom);
        }
    }
    res
}

fn mul_factor(a: &mut Mpz, term_m1: &Mpz, term_last_m1: &Mpz) {
    Mpz::mul_assign(a, term_m1);
    debug_assert!(a.is_multiple_of(&term_last_m1));
    a.set_divexact(&term_last_m1);
}

fn split_primes_factor(prec: usize, expt: u64, res: &mut HmfGen, u_bds: &UBounds) {
    let split_prec = {
        if is_even!(prec) {
            let v1 = prec >> 1;
            5 * v1 * v1
        } else {
            (5 * prec * prec - 1) >> 2
        }
    };

    let mut term = Mpz::new();
    let mut term_m1 = Mpz::new();
    let mut term_last = Mpz::new();
    let mut term_last_m1 = Mpz::new();
    let mut mult = Mpz::new();

    for &p in prime_sieve(split_prec).iter().filter(|&x| {
        (x % 5 == 1) | (x % 5 == 4)
    })
    {
        let mut p_pow = p as i64;

        mult.set_ui(p as u64);
        mult.set_pow_ui(expt);
        term_last.set(&mult);
        term.mul_mut(&term_last, &term_last);


        let elt_p = split_prime_gen(p as i64);
        let mut elt = elt_p.clone();
        let mut p_pow_z = Mpz::new();
        let mut elt_ir_inv = Mpz::new();
        let mut elt_ir_z = Mpz::new();

        while p_pow as usize <= split_prec {
            p_pow_z.set_si(p_pow);
            term_m1.sub_ui_mut(&term, 1);
            term_last_m1.sub_ui_mut(&term_last, 1);

            for (v, &bd) in u_bds.vec.iter().enumerate().skip(1) {
                let bd = bd as i64;
                elt_ir_z.set_si(elt.ir);
                elt_ir_inv.invert_mod_mut(&elt_ir_z, &p_pow_z);
                let mut u_init_val = (elt.rt * elt_ir_inv.into_ui().unwrap() as i64 * v as i64)
                    .abs() % p_pow;
                if u_init_val > (p_pow + 1) >> 1 {
                    u_init_val = p_pow - u_init_val;
                }
                if u_init_val == 0 && is_even!(v) {
                    mul_factor(res.fcvec.fc_ref_mut(v, 0, bd), &term_m1, &term_last_m1);
                }
                let v_i64 = v as i64;
                for u in (0..)
                    .map(|i| u_init_val + p_pow * i)
                    .take_while(|&i| i <= bd)
                    .filter(|i| is_even!(i - v_i64))
                {
                    mul_factor(res.fcvec.fc_ref_mut(v, u, bd), &term_m1, &term_last_m1);
                    mul_factor(res.fcvec.fc_ref_mut(v, -u, bd), &term_m1, &term_last_m1);
                }
                for u in (0..)
                    .map(|i| p_pow - u_init_val + p_pow * i)
                    .take_while(|&i| i <= bd)
                    .filter(|i| is_even!(i - v_i64))
                {
                    mul_factor(res.fcvec.fc_ref_mut(v, u, bd), &term_m1, &term_last_m1);
                    mul_factor(res.fcvec.fc_ref_mut(v, -u, bd), &term_m1, &term_last_m1);
                }
            }
            p_pow *= p as i64;
            term_last.set(&term);
            term *= &mult;
            elt = &elt * &elt_p;
        }
    }
}


fn innert_primes_factor(prec: usize, expt: u64, res: &mut HmfGen, u_bds: &UBounds) {

    let mut term = Mpz::new();
    let mut term_m1 = Mpz::new();
    let mut term_last = Mpz::new();
    let mut term_last_m1 = Mpz::new();
    let mut mult = Mpz::new();

    for &p in prime_sieve(prec).iter().filter(
        |&x| (x % 5 == 2) | (x % 5 == 3),
    )
    {
        let mut p_pow = p;
        mult.set_ui(p as u64);
        mult.set_pow_ui(2 * expt);
        term_last.set(&mult);
        term.mul_mut(&term_last, &term_last);

        while p_pow <= prec {
            term_m1.sub_ui_mut(&term, 1);
            term_last_m1.sub_ui_mut(&term_last, 1);
            // p may be 2.
            for (v, v_divp) in (1..).map(|i| (p_pow * i, i)).take_while(
                |&(i, _)| i <= prec,
            )
            {
                let bd = u_bds.vec[v] as i64;
                let v_divp_i64 = v_divp as i64;

                if is_even!(v_divp) {
                    mul_factor(res.fcvec.fc_ref_mut(v, 0, bd), &term_m1, &term_last_m1);
                }

                for (u, _) in (1..)
                    .map(|i| (p_pow as i64 * i, i))
                    .take_while(|&(i, _)| i <= bd)
                    .filter(|&(_, i)| is_even!(i - v_divp_i64))
                {
                    let u = u as i64;
                    let bd = bd as i64;
                    mul_factor(res.fcvec.fc_ref_mut(v, u, bd), &term_m1, &term_last_m1);
                    mul_factor(res.fcvec.fc_ref_mut(v, -u, bd), &term_m1, &term_last_m1);
                }
            }
            p_pow *= p;
            term_last.set(&term);
            term *= &mult;
        }
    }
}

fn sqrt5_factor(prec: usize, expt: u64, res: &mut HmfGen, u_bds: &UBounds) {
    let mut term = Mpz::new();
    let mut term_m1 = Mpz::new();
    let mut term_last = Mpz::new();
    let mut term_last_m1 = Mpz::new();
    let mut mult = Mpz::new();

    let mut p_pow_u = 5 as i64;
    let mut p_pow_v = 1 as usize;

    mult.set_ui(5);
    mult.set_pow_ui(expt);
    term_last.set(&mult);
    term.mul_mut(&term_last, &term_last);

    let mut i: usize = 1;

    while p_pow_v <= prec {
        term_m1.sub_ui_mut(&term, 1);
        term_last_m1.sub_ui_mut(&term_last, 1);
        for v in (1..).map(|i| p_pow_v * i).take_while(|&i| i <= prec) {
            let bd = u_bds.vec[v] as i64;
            let v_i64 = v as i64;
            if is_even!(v) {
                mul_factor(res.fcvec.fc_ref_mut(v, 0, bd), &term_m1, &term_last_m1);
            }
            for u in (1..).map(|i| p_pow_u * i).take_while(|&i| i <= bd).filter(
                |&i| {
                    is_even!(i - v_i64)
                },
            )
            {
                let u = u as i64;
                let bd = bd as i64;
                mul_factor(res.fcvec.fc_ref_mut(v, u, bd), &term_m1, &term_last_m1);
                mul_factor(res.fcvec.fc_ref_mut(v, -u, bd), &term_m1, &term_last_m1);
            }
        }
        term_last.set(&term);
        term *= &mult;
        if is_even!(i) {
            p_pow_u *= 5;
        } else {
            p_pow_v *= 5;
        }
        i += 1;
    }
}

/// Assuming p splits in Q(sqrt5), return an element alpha s.t.
/// (p, alpha) is a prime abote p.
fn split_prime_gen(p: i64) -> Sqrt5Elt<i64> {
    let mut five_sqrt = 0;
    for i in 2..p {
        if (i * i - 5) % p == 0 {
            if is_even!(i) {
                five_sqrt = p - i;
            } else {
                five_sqrt = i;
            }
        }
    }
    Sqrt5Elt {
        rt: five_sqrt,
        ir: 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ops::AddAssign;

    // #[ignore]
    #[test]
    fn test_prime_element() {
        for &p in prime_sieve(1000).iter().filter(
            |&x| (x % 5 == 1) | (x % 5 == 4),
        )
        {
            let p = p as i64;
            let elt = split_prime_gen(p);
            let sqrd = &elt * &elt;
            assert!(elt.norm() % p == 0);
            assert!(sqrd.norm() % (p * p) == 0);
        }
    }

    #[test]
    fn test_eisenstein() {
        let f = eisenstein_series(4, &Mpz::one(), &Mpz::from_ui(240), 30);
        // ellipti c Eisenstein sereis of weight 8
        let ell_eisen = vec![
            "1",
            "480",
            "61920",
            "1050240",
            "7926240",
            "37500480",
            "135480960",
            "395301120",
            "1014559200",
            "2296875360",
            "4837561920",
            "9353842560",
            "17342613120",
            "30119288640",
            "50993844480",
            "82051050240",
            "129863578080",
            "196962563520",
            "296296921440",
            "429058435200",
            "619245426240",
            "864918850560",
            "1206645690240",
            "1634316215040",
            "2219855529600",
            "2929725000480",
            "3885388234560",
            "5023266412800",
            "6527607394560",
            "8279940628800",
            "10584585480960",
        ];
        let mut a = Mpz::new();
        assert_eq!(f.fcvec.fc_ref(0, 0, 0).to_str_radix(10), "1");
        for (v, &bd) in f.u_bds.vec.iter().enumerate().skip(1) {
            a.set_ui(0);
            let bd = bd as i64;
            let v_i = v as i64;
            for u in u_iter!(v_i, bd) {
                Mpz::add_assign(&mut a, f.fcvec.fc_ref(v, u, bd));
            }
            assert_eq!(a.to_str_radix(10), ell_eisen[v].to_string());
        }
    }
}

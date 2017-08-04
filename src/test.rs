use std::time::Instant;
use theta_chars::{theta, g5_normalized};
use elements::HmfGen;
use eisenstein::eisenstein_series;
use misc::prime_sieve;
use gmp::mpz::Mpz;
use diff_op::{g15_normalized, monom_g2_g6_g10, rankin_cohen_sqrt5, bracket_inner_prod};


// Taken from http://qiita.com/pseudo_foxkeh/items/5d5226e3ffa27631e80d
macro_rules! measure_time {
  ( $x:expr) => {
    {
      let start = Instant::now();
      let result = $x;
      let end = start.elapsed();
      println!("{}.{:03} seconds passed", end.as_secs(), end.subsec_nanos() / 1_000_000);
      result
    }
  };
}

mod g15_part {
    use super::*;
    use diff_op::bracket_proj;
    use misc::PowGen;

    #[test]
    fn test_bracket_proj() {
        let prec = 10;
        let g15 = g15_normalized(prec);
        let mut g5 = g5_normalized(prec);
        let mut e20 = eisenstein_series(10, prec);
        e20.square();
        let mut e2 = eisenstein_series(2, prec);
        let f = &e20 + &(&g5 * &g15);
        let g = &e2 * &g15;
        assert_eq!(f.weight, Some((20, 20)));
        let f5 = bracket_proj(&f, &g15).unwrap();
        let f2 = bracket_proj(&g, &g15).unwrap();
        g5.decrease_prec(f5.prec);
        e2.decrease_prec(f2.prec);
        assert_eq!(g5, f5);
        assert_eq!(e2, f2);
    }
}

mod rankin_cohen {
    use super::*;
    use bignum::{Sqrt5Mpz, RealQuadElement};
    use libc::{c_long};

    #[test]
    fn test_rankin_cohen1() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g6 = eisenstein_series(6, prec);
        let g5 = g5_normalized(prec);
        {
            let f = rankin_cohen_sqrt5(1, &g2, &g2).unwrap();
            assert!(f.is_zero());
        }

        let g15: HmfGen<Sqrt5Mpz> = {
            let g15 = g15_normalized(prec);
            From::from(&g15)
        };

        let g7_9 = rankin_cohen_sqrt5(1, &g2, &g5).unwrap();
        let g8_10 = rankin_cohen_sqrt5(1, &g2, &g6).unwrap();
        let g11_13 = rankin_cohen_sqrt5(1, &g5, &g6).unwrap();
        let f2 = bracket_inner_prod(&g7_9, &g8_10, &g15).unwrap();
        assert!(f2.rt_part().is_zero());
        let mut f2 = f2.ir_part();
        assert!(f2.is_divisible_by_const(&Mpz::from_si(-172800)));
        f2 /= &Mpz::from_si(-172800);
        assert_eq!(f2, g2);

        let f6 = bracket_inner_prod(&g8_10, &g11_13, &g15).unwrap();
        assert!(f6.rt_part().is_zero());
        let mut g6 = g6;
        g6 *= -518400 as c_long;
        let f6 = f6.ir_part();
        assert_eq!(f6, g6);
    }

    #[test]
    fn test_rankin_cohen2() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let f = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        let g = rankin_cohen_sqrt5(2, &g2, &g5).unwrap();
        let g15 = g15_normalized(prec);
        let g15: HmfGen<Sqrt5Mpz> = From::from(&g15);
        let h = bracket_inner_prod(&f, &g, &g15).unwrap();
        assert_eq!(h.prec, prec - 2);
        let mut const_form: HmfGen<Sqrt5Mpz> = HmfGen::one(h.prec);
        const_form *= &h.fourier_coefficient(0, 0);
        assert_eq!(h, const_form);
    }
}

mod g15_squared {
    use super::*;
    use misc::PowGen;

    #[allow(dead_code)]
    fn fc_vec(f: &HmfGen<Mpz>, num: usize) -> Vec<Mpz> {
        let mut tpls = Vec::new();
        v_u_bd_iter!((f.u_bds, v, u, bd) {
            if u >= 0 {
                tpls.push((v, u));
            }
        });
        tpls = tpls.into_iter().take(num).collect();
        f.fourier_coefficients(&tpls)
    }

    #[test]
    fn test_g15() {
        let prec = 10;
        let mut g30 = g15_normalized(prec);
        g30.square();
        assert_eq!(g30.weight, Some((30, 30)));
        let f0 = monom_g2_g6_g10(prec, 0, 0, 3);
        let f1 = monom_g2_g6_g10(prec, 2, 1, 2);
        let f2 = monom_g2_g6_g10(prec, 5, 0, 2);
        let f3 = monom_g2_g6_g10(prec, 1, 3, 1);
        let f4 = monom_g2_g6_g10(prec, 4, 2, 1);
        let f5 = monom_g2_g6_g10(prec, 7, 1, 1);
        let f6 = monom_g2_g6_g10(prec, 10, 0, 1);
        let f7 = monom_g2_g6_g10(prec, 0, 5, 0);
        let f8 = monom_g2_g6_g10(prec, 3, 4, 0);
        let f9 = monom_g2_g6_g10(prec, 6, 3, 0);
        let f10 = monom_g2_g6_g10(prec, 9, 2, 0);
        let f11 = monom_g2_g6_g10(prec, 12, 1, 0);
        let f12 = monom_g2_g6_g10(prec, 15, 0, 0);
        let forms = vec![f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, g30];
        for f in forms.iter() {
            assert_eq!(f.weight, Some((30, 30)));
        }
        let v = vec![
            Mpz::from_str_radix("272097792000000000000000", 10).unwrap(),
            Mpz::from_str_radix("251942400000000000", 10).unwrap(),
            Mpz::from_str_radix("-11438184960000000000", 10).unwrap(),
            Mpz::from_str_radix("-972000000", 10).unwrap(),
            Mpz::from_str_radix("172044000000", 10).unwrap(),
            Mpz::from_str_radix("-9963972000000", 10).unwrap(),
            Mpz::from_str_radix("187622244000000", 10).unwrap(),
            Mpz::from_str_radix("1", 10).unwrap(),
            Mpz::from_str_radix("-310", 10).unwrap(),
            Mpz::from_str_radix("38190", 10).unwrap(),
            Mpz::from_str_radix("-2334280", 10).unwrap(),
            Mpz::from_str_radix("70679305", 10).unwrap(),
            Mpz::from_str_radix("-846347082", 10).unwrap(),
            Mpz::from_str_radix("-87071293440000000000", 10).unwrap(),
        ];
        let mut tmp = HmfGen::new(prec);
        let mut tmp1 = HmfGen::new(prec);
        tmp.mul_mut_by_const(&forms[0], &v[0]);
        for (f, a) in forms.iter().zip(v.iter()).skip(1) {
            tmp1.mul_mut_by_const(&f, &a);
            tmp += &tmp1;
        }
        assert!(tmp.is_zero());
    }
}

mod theta_eisen_relatioin {
    use super::*;

    #[test]
    fn test_relation() {
        let prec = 10;
        let mut g5_normalized = theta(prec);
        assert!(g5_normalized.is_divisible_by_const(&Mpz::from_ui(64)));
        g5_normalized /= &Mpz::from_ui(64);
        println!("{}", g5_normalized);
        let e2 = eisenstein_series(2, prec);
        let e6 = eisenstein_series(6, prec);
        let e10 = eisenstein_series(10, prec);
        let mut tmp = HmfGen::new(prec);
        tmp.pow_mut(&e2, 5);
        let f10_1 = tmp.clone();
        tmp.pow_mut(&e2, 2);
        let f10_2 = &tmp * &e6;
        let g5_sq = &g5_normalized * &g5_normalized;
        let mut g10 = &(&e10 + &(&f10_1 * &Mpz::from_ui(355404))) +
            &(&f10_2 * &Mpz::from_si(-11465));
        assert!(g10.is_divisible_by_const(&Mpz::from_ui(5315625)));
        g10 /= &(&Mpz::from_ui(5315625) * &Mpz::from_ui(1024));
        assert_eq!(&g10, &g5_sq);
    }
}

mod elements {
    use super::*;

    #[test]
    fn test_pow() {
        let e2 = eisenstein_series(2, 20);
        let e4 = eisenstein_series(4, 20);
        let mut tmp = HmfGen::new(20);
        let f = &e2 * &e2;
        assert_eq!(f, e4);
        let g = &f * &e2;
        tmp.pow_mut(&e2, 3);
        assert_eq!(tmp, g);
        let g = &e2 * &(&e2 * &(&e2 * &(&e2 * &(&e2 * &e2))));
        tmp.pow_mut(&e2, 6);
        assert_eq!(tmp, g);
        let g = &e2 * &(&e2 * &(&e2 * &(&e2 * &(&e2 * &(&e2 * &e2)))));
        tmp.pow_mut(&e2, 7);
        assert_eq!(tmp, g);
    }

    #[test]
    fn test_op() {
        let e2 = eisenstein_series(2, 20);
        let mut tmp = HmfGen::new(20);
        tmp.mul_mut_by_const(&e2, &Mpz::from_ui(2));
        assert_eq!(tmp, &e2 + &e2);
        assert!(!tmp.is_zero());
        tmp.mul_mut_by_const(&e2, &Mpz::from_ui(0));
        assert!(tmp.is_zero());
        assert!((&e2 - &e2).is_zero());
    }

    #[test]
    fn test_op_time() {
        let f = eisenstein_series(2, 50);
        let g = f.clone();
        measure_time!(&f * &g);
    }
}

mod theta_fast {
    use super::*;

    #[test]
    fn test_g5() {
        let prec = 15;
        let f = g5_normalized(prec);
        let mut g = theta(prec);
        g /= &Mpz::from_ui(64);
        assert_eq!(f.fcvec.vec, g.fcvec.vec);
    }
}

mod theta_char {
    use super::*;

    #[test]
    fn theta_fun() {
        measure_time!(theta(10));
    }

    #[test]
    fn theta_fun1() {
        measure_time!(g5_normalized(50));
    }
}

mod eisen {
    use super::*;
    use std::ops::AddAssign;

    #[test]
    fn test_e2_sqaured() {
        let e2 = eisenstein_series(2, 20);
        let e4 = eisenstein_series(4, 20);
        let mut f = HmfGen::new(20);
        measure_time!(f.mul_mut(&e2, &e2));
        v_u_bd_iter!((f.u_bds, v, u, bd) {
            assert_eq!(f.fcvec.fc_ref(v, u, bd), e4.fcvec.fc_ref(v, u, bd));
        }
        )
    }

    #[test]
    fn test_eisenstein_diagonal() {
        let f = eisenstein_series(4, 100);
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
            "13206054773760",
            "16622537994720",
            "20466207521280",
            "25408170694080",
            "30883295301120",
            "37928302819680",
            "45567301024320",
            "55348538140800",
            "65901003544320",
            "79263452059200",
            "93482051463360",
            "111574531722240",
            "130472933331840",
            "154460002193280",
            "179445684375360",
            "210826791740160",
            "243179097822720",
            "284141508839040",
            "325547470268640",
            "377934525061920",
            "430954088981760",
            "497359813312320",
            "563861347122240",
            "648001367251200",
            "730778303842560",
            "835534141804800",
            "938779856217600",
            "1068112341115200",
            "1194552712713600",
            "1354908992613120",
            "1508516561290560",
            "1703581065815040",
            "1891577921475840",
            "2127684863324640",
            "2353099544288640",
            "2640140770245120",
            "2909141570555520",
            "3252442811405760",
            "3575883878507520",
            "3983945093844480",
            "4365657676028160",
            "4854825057794400",
            "5302751289167040",
            "5878181832137280",
            "6410238301050240",
            "7085041940457600",
            "7703300917232640",
            "8501229457217280",
            "9217876313356800",
            "10145721901078080",
            "10985883644794080",
            "12059184638773440",
            "13025304475021440",
            "14282404979297280",
            "15387897237563520",
            "16831008399807360",
            "18116510095814400",
            "19770889634582400",
            "21231040749854400",
            "23148493284421440",
            "24804559443740160",
            "26987463658955520",
            "28894847844986880",
            "31370103619130880",
            "33520619308435200",
            "36370113132447360",
            "38783176549494720",
            "41995623664654560",
            "44759605202881920",
            "48378548932926240",
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

    #[test]
    fn test_eisen() {
        measure_time!(eisenstein_series(4, 30));
    }
}

mod misc {
    use super::*;

    #[test]
    fn prime_sieve_test() {
        let n = 1009;
        let v = prime_sieve(n);
        assert_eq!(v.len(), 169);
    }
}

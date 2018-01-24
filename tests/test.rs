extern crate hilbert_sqrt5;
extern crate gmp;
extern crate bincode;
extern crate libc;
extern crate csv;

use std::time::Instant;
use hilbert_sqrt5::theta_chars::{theta, g5_normalized};
use hilbert_sqrt5::elements::HmfGen;
use hilbert_sqrt5::eisenstein::{eisenstein_series, f6_normalized};
use hilbert_sqrt5::misc::prime_sieve;
use gmp::mpz::Mpz;
use hilbert_sqrt5::diff_op::{g15_normalized, rankin_cohen_sqrt5, bracket_inner_prod,
                             bracket_inner_prod1};
macro_rules! is_even {
    ($expr: expr) => {($expr) & 1 == 0}
}


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

mod diag_res {
    use hilbert_sqrt5::eisenstein::{eisenstein_series, f6_normalized};
    use hilbert_sqrt5::diff_op::g15_normalized;

    #[test]
    fn test_res() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g6 = f6_normalized(prec);
        let g15 = g15_normalized(prec);
        println!("{:?}", g2.diagonal_restriction());
        println!("{:?}", g6.diagonal_restriction());
        println!("{:?}", g15.diagonal_restriction());
    }
}

mod str_exe {
    use hilbert_sqrt5::structure::*;
    use hilbert_sqrt5::eisenstein::{eisenstein_series, f6_normalized};
    use hilbert_sqrt5::theta_chars::g5_normalized;
    use hilbert_sqrt5::elements::HmfGen;
    use std::fs::File;
    use hilbert_sqrt5::bignum::Sqrt5Mpz;
    use hilbert_sqrt5::bignum::BigNumber;
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_tpls_of_wt() {
        assert_eq!(tpls_of_wt(30).len(), 13);
        assert_eq!(tpls_of_wt(100).len(), 99);
    }

    fn brackets(forms: &[HmfGen<Sqrt5Mpz>]) -> Vec<(PWtPolyZ, Mpz)> {
        forms
            .iter()
            .enumerate()
            .flat_map(|(i, f)| {
                forms
                    .iter()
                    .skip(i + 1)
                    .map(|g| bracket_inner_prod_as_pol_over_z_maybe(f, g).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    #[test]
    fn test_save_rels_up_to_50() {
        for i in 1..61 {
            println!("{}", i);
            let prec = (2 * i + 6) / 5 + 2;
            println!("{}", prec);
            let forms_w_monoms = mixed_weight_forms(i, prec, 6);
            let forms: Vec<_> = forms_w_monoms.clone().into_iter().map(|f| f.0).collect();
            let weight: Vec<_> = forms_w_monoms
                .iter()
                .map(|f| f.0.weight.unwrap().0)
                .collect();
            let monoms: Vec<_> = forms_w_monoms.iter().map(|f_t| (f_t.1, f_t.2)).collect();
            println!("{:?}", weight);
            let ref mut monms_file = File::create(format!("./data/brackets/str{}_monoms.sobj", i))
                .unwrap();
            let ref mut f_wt = File::create(format!("./data/brackets/str{}_weights.sobj", i))
                .unwrap();
            save_as_pickle(&monoms, monms_file);
            save_as_pickle(weight, f_wt);
            let brs = brackets(&forms);
            let ref mut f = File::create(format!("./data/brackets/str{}_brs.sobj", i)).unwrap();
            save_polys_over_z_pickle(&brs, f);
        }
    }

    #[test]
    fn test_into_form() {
        let prec = 10;
        let monoms = monoms_of_g2_g5_f6(40);
        let monoms1 = monoms_of_g2_g5_f6(42);
        let mut map_g2 = HashMap::new();
        let mut map_g5 = HashMap::new();
        let mut map_g6 = HashMap::new();
        for m in monoms.iter() {
            let f = m.into_form(prec);
            let g = m.into_form_cached(prec, &mut map_g2, &mut map_g5, &mut map_g6);
            assert_eq!(f, g);
        }
        for m in monoms1.iter() {
            let f = m.into_form(prec);
            let g = m.into_form_cached(prec, &mut map_g2, &mut map_g5, &mut map_g6);
            assert_eq!(f, g);
        }
    }

    /// Test if weight 2 forms are holomorphic by multplying g5.
    #[test]
    fn test_weight2_hol() {
        let prec = 10;
        let g5 = g5_normalized(prec);
        for i in 3..51 {
            println!("{}", i);
            let cand = {
                let cand_f = File::open(format!("./data/brackets/str{}_cand.sobj", i)).unwrap();
                let monom_f = File::open(format!("./data/brackets/str{}_monoms.sobj", i)).unwrap();
                StrCand::load(i, &cand_f, &monom_f).unwrap()
            };
            if cand.gen_wts.iter().any(|&x| x == 2) {
                let gens = cand.gens_normalized(prec);
                for f in &gens {
                    if f.weight.unwrap().0 == 2 {
                        let g7 = f * &Into::<HmfGen<Sqrt5Mpz>>::into(&g5);
                        assert!(g7.diagonal_restriction().iter().all(|x| x.is_zero_g()));
                    }
                }
            }
            println!("{}: {:?}", i, cand.gen_wts);
        }
    }

    fn write_csv_form(
        f: &HmfGen<Sqrt5Mpz>,
        p: &String,
        num: Option<Sqrt5Mpz>,
        den: Option<Sqrt5Mpz>,
    ) {
        let mut wtr = csv::Writer::from_path(p).unwrap();
        if let Some(num_val) = num {
            wtr.write_record(&["numerator", &format!("{}", num_val)])
                .unwrap();
            wtr.write_record(&["denominator", &format!("{}", den.unwrap())])
                .unwrap();
        }
        wtr.write_record(&["(v, u)", "a(v, u)"]).unwrap();
        for (v, &bd) in f.u_bds.vec.iter().enumerate() {
            let bd = bd as i64;
            let v_i = v as i64;
            for u in (-bd..(bd + 1)).filter(|&x| is_even!(x + v_i)) {
                wtr.write_record(
                    &[
                        format!("{:?}", (v, u)),
                        format!("{}", f.fourier_coefficient(v, u)),
                    ],
                ).unwrap();
            }
        }
        wtr.flush().unwrap();
    }

    #[test]
    fn write_gens_parallel_wt() {
        let prec = 15;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let f6 = f6_normalized(prec);
        let g15 = g15_normalized(prec);
        let gens = vec![
            From::from(&g2),
            From::from(&g5),
            From::from(&f6),
            From::from(&g15),
        ];
        let paths = vec![
            "./forms_csv/A0/g2.csv",
            "./forms_csv/A0/g5.csv",
            "./forms_csv/A0/g6.csv",
            "./forms_csv/A0/g15.csv",
        ];
        for (p, f) in paths.into_iter().zip(gens.iter()) {
            write_csv_form(f, &p.to_string(), None, None);
        }
    }

    #[test]
    fn write_gens_cands() {
        for i in 51..61 {
            let cand = {
                let cand_f = File::open(format!("./data/brackets/str{}_cand.sobj", i)).unwrap();
                let monom_f = File::open(format!("./data/brackets/str{}_monoms.sobj", i)).unwrap();
                StrCand::load(i, &cand_f, &monom_f).unwrap()
            };
            let prec = 15;
            let gens = cand.gens_with_const(prec);
            println!("{}", i);
            for (n, &(ref f, ref num, ref den)) in gens.iter().enumerate() {
                let (w1, w2) = f.weight.unwrap();
                let path = format!("./forms_csv/A{}/gen{}_wt_{}_{}.csv", i, n, w1, w2);
                write_csv_form(f, &path, Some(num.clone()), Some(den.clone()));
            }
        }
    }
}

mod str_test {
    use hilbert_sqrt5::structure::*;
    use hilbert_sqrt5::eisenstein::eisenstein_series;
    use hilbert_sqrt5::bignum::RealQuadElement;
    use super::*;

    #[test]
    fn test_relations_over_z() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let f = &g2 * 2;
        let forms = [f, g2];
        let rels = relations_over_z(&forms);
        assert_eq!(rels.len(), 1);
    }

    #[test]
    fn test_relations() {
        for i in 3..21 {
            let cand = load_cand(i);
            println!("{}", i);
            cand.test_relations(15);
        }
    }

    #[test]
    fn test_bracket_as_form() {
        for i in 3..41 {
            println!("{}", i);
            let prec = 5;
            let cand = load_cand(i);
            let gens = cand.gens_normalized(prec);
            let f = &gens[0];
            let g = &gens[1];
            let prec = (f.weight.unwrap().0 + g.weight.unwrap().1 - 15) / 5 + 5;
            let gens = cand.gens_normalized(prec);
            let f = &gens[0];
            let g = &gens[1];
            let br = bracket_inner_prod1(&f, &g).unwrap();
            let f1 = br.ir_part();
            let f2 = br.rt_part();
            if !f1.is_zero() {
                let pl = r_elt_as_pol_over_z(&f1);
                assert!(pl.is_some());
                // println!("{:?}", pl.unwrap());
            }
            if !f2.is_zero() {
                let pl = r_elt_as_pol_over_z(&f2);
                assert!(pl.is_some());
                // println!("{:?}", pl.unwrap());

            }
        }
    }
}

mod div {
    use super::*;
    use hilbert_sqrt5::elements::{div_mut, div_mut_with_denom};

    #[test]
    fn test_div() {
        let prec = 10;
        let f6 = f6_normalized(prec);
        let g5 = g5_normalized(prec);
        let mut f = &f6 * &g5;
        f *= &Mpz::from_ui(2);
        let mut res = HmfGen::new(prec);
        div_mut(&mut res, &f, &g5);
        let mut g6 = f6.clone();
        g6.decrease_prec(res.prec);
        assert_eq!(&g6 * &Mpz::from_ui(2), res);
    }

    #[test]
    fn test_div2() {
        let prec = 10;
        let h = eisenstein_series(2, prec);
        let g = g15_normalized(prec);
        let f = &g * &h;
        let mut res = HmfGen::new(prec);
        div_mut(&mut res, &f, &g);
        let mut h = h.clone();
        h.decrease_prec(res.prec);
        assert_eq!(h, res);
    }

    #[test]
    fn test_div_with_dnm() {
        let prec = 10;
        let g = eisenstein_series(6, prec);
        let h = g5_normalized(prec);
        let f = &g * &h;
        let mut res = HmfGen::new(prec);
        let dnm = div_mut_with_denom(&mut res, &f, &g, true);
        println!("{}", dnm);
        assert_eq!(res, &h * &dnm);
        let a = res.gcd();
        res /= &a;
        assert_eq!(res, h);
    }
}

mod structure {
    use super::*;
    use hilbert_sqrt5::structure::{relation, monoms_of_g2_g5_f6};
    use hilbert_sqrt5::bignum::{Sqrt5Mpz, RealQuadElement};

    #[test]
    fn test_relation() {
        let prec = 5;
        let f = eisenstein_series(6, prec);
        let f: HmfGen<Sqrt5Mpz> = From::from(&f);
        let forms_monom = monoms_of_g2_g5_f6(6);
        let v = {
            let forms: Vec<_> = forms_monom
                .iter()
                .map(|x| From::from(&x.into_form(prec)))
                .collect();
            relation(10, &f, &forms)
        };
        let mut f6 = HmfGen::new(prec);
        let mut tmp = HmfGen::new(prec);
        for (f, a) in forms_monom.iter().zip(v.iter().skip(1)) {
            let mut a: Mpz = a.rt_part();
            a >>= 1;
            if !a.is_zero() {
                tmp.mul_mut_by_const(&f.into_form(prec), &a);
                f6 += &tmp;
            }
        }
        f6.negate();
        let f6: HmfGen<Sqrt5Mpz> = From::from(&f6);
        println!("{:?}", v);
        assert_eq!(f, f6);
    }
}

mod serialize {
    use super::*;

    use bincode::{serialize, deserialize, Infinite};

    #[test]
    fn test_serialize() {
        let f = eisenstein_series(2, 30);
        let serialized = serialize(&f, Infinite).unwrap();
        let deserialized: HmfGen<Mpz> = deserialize(&serialized).unwrap();
        assert_eq!(f, deserialized);
    }
}

mod g15_part {
    use super::*;
    use hilbert_sqrt5::diff_op::bracket_proj;
    use hilbert_sqrt5::misc::PowGen;

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
    use hilbert_sqrt5::bignum::{Sqrt5Mpz, RealQuadElement};
    use libc::c_long;
    use hilbert_sqrt5::misc::PowGen;
    use hilbert_sqrt5::structure::{relation, monoms_of_g2_g5_f6};
    use hilbert_sqrt5::elements::div_mut;

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

    #[test]
    fn test_rankin_cohen4() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g7_15 = rankin_cohen_sqrt5(4, &g2, &g5).unwrap();
        let g2_sqrt5: HmfGen<Sqrt5Mpz> = From::from(&g2);
        let mut g5_13 = HmfGen::new(prec);
        div_mut(&mut g5_13, &g7_15, &g2_sqrt5);
        assert_eq!(g7_15, &g2_sqrt5 * &g5_13);
        assert!(g5_13.is_divisible_by_const(&Sqrt5Mpz::from_sisi(35, -15)));
        g5_13 /= &Sqrt5Mpz::from_sisi(35, -15);
        println!("{}", g5_13);
        assert_eq!(g7_15.weight, Some((7, 15)));
        let mut g8_16 = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        g8_16.square();
        assert_eq!(g8_16.weight, Some((8, 16)));
        let h = bracket_inner_prod1(&g7_15, &g8_16).unwrap();
        let forms_monom = monoms_of_g2_g5_f6(8);
        let v = {
            let forms: Vec<_> = forms_monom
                .iter()
                .map(|x| From::from(&x.into_form(prec)))
                .collect();
            relation(60, &h, &forms)
        };
        let indics: Vec<_> = forms_monom.iter().map(|x| x.idx).collect();
        println!("{:?}", v);
        println!("{:?}", indics);
    }
}

mod g15_squared {
    use super::*;
    use hilbert_sqrt5::misc::PowGen;
    use hilbert_sqrt5::structure::{relation_monom, MonomFormal};
    use hilbert_sqrt5::bignum::{Sqrt5Mpz, BigNumber};

    // #[allow(dead_code)]
    // fn fc_vec(f: &HmfGen<Mpz>, num: usize) -> Vec<Mpz> {
    //     let mut tpls = Vec::new();
    //     v_u_bd_iter!((f.u_bds, v, u, bd) {
    //         if u >= 0 {
    //             tpls.push((v, u));
    //         }
    //     });
    //     tpls = tpls.into_iter().take(num).collect();
    //     f.fourier_coefficients(&tpls)
    // }

    #[test]
    fn test_g15_relation() {
        let prec = 10;
        let mut g30 = g15_normalized(prec);
        g30.square();
        assert_eq!(g30.weight, Some((30, 30)));
        let g30: HmfGen<Sqrt5Mpz> = From::from(&g30);
        let (_, v) = relation_monom(100, &g30);
        let mut f30 = MonomFormal::eval(&v, prec);
        f30.negate();
        assert!(f30.is_divisible_by_const(&Sqrt5Mpz::from_ui_g(16)));
        f30 >>= 4;
        assert_eq!(g30, f30);
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

    #[test]
    fn test_e2_sqaured() {
        let e2 = eisenstein_series(2, 20);
        let e4 = eisenstein_series(4, 20);
        let mut f = HmfGen::new(20);
        measure_time!(f.mul_mut(&e2, &e2));
        assert_eq!(f, e4);
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
        let v = f.diagonal_restriction();
        let w: Vec<String> = ell_eisen
            .iter()
            .map(|&x| x.to_string())
            .take(v.len())
            .collect();
        let v: Vec<String> = v.iter().map(|x| x.to_str_radix(10)).collect();
        assert_eq!(v, w);
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


mod paper {
    use super::*;
    use hilbert_sqrt5::bignum::Sqrt5Mpz;
    use hilbert_sqrt5::structure::*;
    use hilbert_sqrt5::bignum::BigNumber;
    use hilbert_sqrt5::bignum::RealQuadElement;
    #[test]
    fn test_generators_parallel() {
        let prec = 2;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let g15 = g15_normalized(prec);
        print_fc(&g2);
        print_fc(&g5);
        print_fc(&g6);
        print_fc(&g15);

        println!("diag: {:?}", g2.diagonal_restriction());
        println!("diag: {:?}", g5.diagonal_restriction());
        println!("diag: {:?}", g6.diagonal_restriction());
        println!("diag: {:?}", g15.diagonal_restriction());
    }

    fn print_fc<T>(f: &HmfGen<T>)
    where
        T: BigNumber + ::std::fmt::Debug,
    {
        println!(
            "{:?}",
            f.fourier_coefficients(&vec![(0, 0), (1, 1), (1, -1), (2, 4), (2, 2), (2, 0)])
        );
    }

    fn g_4_8(prec: usize) -> HmfGen<Sqrt5Mpz> {
        let g2 = eisenstein_series(2, prec);
        let mut g_4_8 = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        g_4_8 /= &Into::<Sqrt5Mpz>::into(&Mpz::from_ui(360));
        assert_eq!(g_4_8.fourier_coefficient(1, 1), From::from((6, -2)));
        g_4_8
    }

    #[test]
    fn test_gens_a2() {
        let prec = 3;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let mut g_4_8 = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        let mut g_7_11 = rankin_cohen_sqrt5(2, &g2, &g5).unwrap();
        println!("{}", g_4_8);
        println!("\n");
        println!("{}", g_7_11);
        println!("\n");

        g_4_8 /= &Into::<Sqrt5Mpz>::into(&Mpz::from_ui(360));
        g_7_11 <<= 1;
        g_7_11 /= &Into::<Sqrt5Mpz>::into(&Mpz::from_ui(3));
        let br = bracket_inner_prod1(&g_4_8, &g_7_11).unwrap();
        println!("br: {}", br);

        println!("diag {:?}", g_4_8.diagonal_restriction());
        println!("diag {:?}", g_7_11.diagonal_restriction());

        assert_eq!(g_4_8.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(g_4_8.fourier_coefficient(1, 1), From::from((6, -2)));
        assert_eq!(g_4_8.fourier_coefficient(1, -1), From::from((6, 2)));
        assert_eq!(g_4_8.fourier_coefficient(2, 4), From::from((36, -16)));
        assert_eq!(g_4_8.fourier_coefficient(2, 2), From::from((-240, 80)));
        assert_eq!(g_4_8.fourier_coefficient(2, 0), From::from((120, 0)));

        assert_eq!(g_7_11.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(g_7_11.fourier_coefficient(1, 1), From::from((6, -2)));
        assert_eq!(g_7_11.fourier_coefficient(1, -1), From::from((-6, -2)));
        assert_eq!(g_7_11.fourier_coefficient(2, 4), From::from((-36, 16)));
        assert_eq!(g_7_11.fourier_coefficient(2, 2), From::from((-240, 80)));
        assert_eq!(g_7_11.fourier_coefficient(2, 0), From::from((0, 960 * 2)));
    }

    #[test]
    fn test_gens_a1() {
        let prec = 5;

        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let mut g_7_9 = rankin_cohen_sqrt5(1, &g2, &g5).unwrap();
        let mut g_8_10 = rankin_cohen_sqrt5(1, &g2, &g6).unwrap();
        let mut g_11_14 = rankin_cohen_sqrt5(1, &g5, &g6).unwrap();
        g_7_9 *= &Into::<Sqrt5Mpz>::into((0, -2));
        g_8_10 *= &Into::<Sqrt5Mpz>::into((0, -2));
        g_11_14 *= &Into::<Sqrt5Mpz>::into((0, -4));
        print!("g_7_9:");
        print_fc(&g_7_9);
        print!("g_8_10:");
        print_fc(&g_8_10);
        print!("g_11_14:");
        print_fc(&g_11_14);

        let f7 = g_7_9;
        let f8 = g_8_10;
        let f11 = g_11_14;

        println!("diag {:?}", f7.diagonal_restriction());
        println!("diag {:?}", f8.diagonal_restriction());
        println!("diag {:?}", f11.diagonal_restriction());

        assert_eq!(f7.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(f7.fourier_coefficient(1, 1), From::from((-10, 2)));
        assert_eq!(f7.fourier_coefficient(1, -1), From::from((10, 2)));
        assert_eq!(f7.fourier_coefficient(2, 4), From::from((20, -8)));
        assert_eq!(f7.fourier_coefficient(2, 2), From::from((2000, -400)));
        assert_eq!(f7.fourier_coefficient(2, 0), From::from((0, 840 * 2)));

        assert_eq!(f8.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(f8.fourier_coefficient(1, 1), From::from((-10, 2)));
        assert_eq!(f8.fourier_coefficient(1, -1), From::from((-10, -2)));
        assert_eq!(f8.fourier_coefficient(2, 4), From::from((-20, 8)));
        assert_eq!(f8.fourier_coefficient(2, 2), From::from((2000, -400)));
        assert_eq!(f8.fourier_coefficient(2, 0), From::from((2 * 3300, 0)));

        assert_eq!(f11.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(f11.fourier_coefficient(1, 1), From::from((0, 0)));
        assert_eq!(f11.fourier_coefficient(1, -1), From::from((0, 0)));
        assert_eq!(f11.fourier_coefficient(2, 4), From::from((0, 0)));
        assert_eq!(f11.fourier_coefficient(2, 2), From::from((10, -2)));
        assert_eq!(f11.fourier_coefficient(2, 0), From::from((0, -44)));

        // let g2: HmfGen<Sqrt5Mpz> = From::from(&eisenstein_series(2, prec));
        // let f11_1 = &f7 * &g2.pow(2);
        // println!("{}", f11_1);
        let br7_8 = bracket_inner_prod1(&f7, &f8).unwrap();
        let br8_11 = bracket_inner_prod1(&f8, &f11).unwrap();
        let br11_7 = bracket_inner_prod1(&f11, &f7).unwrap();
        println!("br7_8: {}", br7_8.fourier_coefficient(0, 0));
        println!("br8_11: {}", br8_11.fourier_coefficient(1, 1));
        println!("br11_7: {}", br11_7.fourier_coefficient(1, 1));
    }

    #[test]
    fn test_gens_a3() {
        let prec = 5;

        let cand = load_cand(3);
        let gens = cand.gens_normalized(prec);
        let g2: HmfGen<Sqrt5Mpz> = From::from(&eisenstein_series(2, prec));
        let g5: HmfGen<Sqrt5Mpz> = From::from(&g5_normalized(prec));

        let g_4_8 = g_4_8(prec);

        {
            let g2 = eisenstein_series(2, prec);
            let g6 = f6_normalized(prec);
            let mut g_6_12 = rankin_cohen_sqrt5(1, &g_4_8.ir_part(), &g2).unwrap();
            g_6_12 *= &Into::<Sqrt5Mpz>::into((0, 2));
            g_6_12 += &rankin_cohen_sqrt5(1, &g_4_8.rt_part(), &g2).unwrap();
            g_6_12 /= &Into::<Sqrt5Mpz>::into((4, 0));
            g_6_12 *= &Into::<Sqrt5Mpz>::into((0, -2));
            g_6_12 /= &Into::<Sqrt5Mpz>::into((8, 0));
            print_fc(&g_6_12);

            let mut g_8_14 = rankin_cohen_sqrt5(3, &g2, &g6).unwrap();
            g_8_14 *= &Into::<Sqrt5Mpz>::into((0, -2));
            g_8_14 /= &Into::<Sqrt5Mpz>::into((8, 0));
            print_fc(&g_8_14);
        }

        let mut f3 = gens[0].clone();
        let mut f6 = gens[1].clone();
        f3 *= &Into::<Sqrt5Mpz>::into((-10, 4));
        f6 *= &Into::<Sqrt5Mpz>::into((10, -4));
        let f8: HmfGen<Sqrt5Mpz> = &(&(&f3 * &g5) * &Into::<Sqrt5Mpz>::into((840 * 2, 0))) -
            &f6 * &g2;
        assert_eq!(f6.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(f6.fourier_coefficient(1, 1), From::from((10, -4)));
        assert_eq!(f6.fourier_coefficient(1, -1), From::from((10, 4)));
        assert_eq!(f6.fourier_coefficient(2, 4), From::from((85 * 2, -38 * 2)));
        assert_eq!(
            f6.fourier_coefficient(2, 2),
            From::from((-2200 * 2, 880 * 2))
        );
        assert_eq!(f6.fourier_coefficient(2, 0), From::from((-1050 * 2, 0)));
        println!("{:?}", f6.diagonal_restriction());

        assert_eq!(f8.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(f8.fourier_coefficient(1, 1), From::from((-10, 4)));
        assert_eq!(f8.fourier_coefficient(1, -1), From::from((-10, -4)));
        assert_eq!(f8.fourier_coefficient(2, 4), From::from((-85 * 2, 38 * 2)));
        assert_eq!(
            f8.fourier_coefficient(2, 2),
            From::from((-2600 * 2, 1040 * 2))
        );
        assert_eq!(f8.fourier_coefficient(2, 0), From::from((8250 * 2, 0)));
        println!("{:?}", f8.diagonal_restriction());

        assert_eq!(f3.fourier_coefficient(0, 0), From::from((0, 0)));
        assert_eq!(f3.fourier_coefficient(1, 1), From::from((-10, 4)));
        assert_eq!(f3.fourier_coefficient(1, -1), From::from((10, 4)));
        assert_eq!(f3.fourier_coefficient(2, 4), From::from((170, -76)));
        assert_eq!(f3.fourier_coefficient(2, 2), From::from((-100, 40)));
        assert_eq!(f3.fourier_coefficient(2, 0), From::from((0, -120)));
        println!("{:?}", f3.diagonal_restriction());
        let a = &g2 * &f6 + &f8 - &g5 * &(&f3 * &Into::<Sqrt5Mpz>::into((840 * 2, 0)));
        assert!(a.is_zero());

        let br3_6 = bracket_inner_prod1(&f3, &f6).unwrap();
        println!("br:{}", br3_6);
    }
}

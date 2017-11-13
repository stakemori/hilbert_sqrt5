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

    #[test]
    fn test_res() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g6 = f6_normalized(prec);
        println!("{:?}", g2.diagonal_restriction());
        println!("{:?}", g6.diagonal_restriction());
    }
}

mod str_exe {
    use hilbert_sqrt5::structure::*;
    use hilbert_sqrt5::eisenstein::{eisenstein_series, f6_normalized};
    use hilbert_sqrt5::theta_chars::g5_normalized;
    use hilbert_sqrt5::elements::HmfGen;
    use std::fs::File;
    use hilbert_sqrt5::diff_op::{rankin_cohen_sqrt5, star_op};
    use hilbert_sqrt5::bignum::Sqrt5Mpz;
    use hilbert_sqrt5::bignum::BigNumber;
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_tpls_of_wt() {
        assert_eq!(tpls_of_wt(30).len(), 13);
        assert_eq!(tpls_of_wt(100).len(), 99);
    }

    #[test]
    fn relation_slow1() {
        println!("{:?}", Structure1::relations());
        let gens = Structure1::gens(10);
        print_3rel(&relation_slow_3gens(&gens, 50));
    }

    #[test]
    fn check_relations1() {
        assert!(Structure1::check_relations(10));
    }

    #[test]
    fn relation_slow4() {
        let gens = Structure4::gens(10);
        print_3rel(&relation_slow_3gens(&gens, 50));
    }

    #[test]
    fn check_relations4() {
        assert!(Structure4::check_relations(10));
    }

    #[test]
    fn check_relations3() {
        assert!(Structure3::check_relations(10));
    }

    #[test]
    fn check_relations5() {
        assert!(Structure5::check_relations(15));
    }

    #[test]
    fn check_relations6() {
        assert!(Structure6::check_relations(10));
    }

    #[test]
    fn check_relations7() {
        assert!(Structure7::check_relations(15));
    }

    #[test]
    fn test_pickle_gen3() {
        let gens1 = Structure3::gens1(10);
        let rel1 = relation_slow_3gens(&gens1, 50);
        let ref mut f1 = File::create("./data/str3gens1.sobj").unwrap();
        save_as_pickle_rel(&rel1, f1);
        let _gens = Structure3::gens(13);
        // let rel = relation_slow_3gens(&gens, 50);
        // let ref mut f = File::create("./data/str3gens.sobj").unwrap();
        // save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_gens3() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f = rankin_cohen_sqrt5(3, &g2, &g6).unwrap();
        let gens = Structure3::gens(prec);
        let f3 = &gens[0];
        let f6 = &gens[1];
        let g5: HmfGen<Sqrt5Mpz> = From::from(&g5);
        let g2: HmfGen<Sqrt5Mpz> = From::from(&g2);
        let h0 = f3 * &g5;
        let h1 = f6 * &g2;

        let forms = vec![f, h0, h1];
        let rels = relations(50, &forms);
        assert!(!rels[0][0].is_zero_g())
    }

    #[test]
    fn test_pickle_gen6_0() {
        let gens = Structure6::gens(10);
        let rel = relation_slow_3gens(&gens, 50);
        for f in gens.iter() {
            assert!(!f.is_zero());
        }
        let ref mut f = File::create("./data/str6gens.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen6_1() {
        let gens1 = Structure6::gens1(10);
        let rel1 = relation_slow_3gens(&gens1, 50);
        let ref mut f1 = File::create("./data/str6gens1.sobj").unwrap();
        save_as_pickle_rel(&rel1, f1);
    }

    #[test]
    fn test_pickle_gen6_2() {
        let gens = Structure6::gens2(10);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str6gens2.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen5_1() {
        let gens = Structure5::gens1(10);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str5gens1.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen5_2() {
        let gens = Structure5::gens2(10);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str5gens2.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen5_0() {
        let gens = Structure5::gens(10);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str5gens.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen7_1() {
        let gens = Structure7::gens1(10);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str7gens1.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen7_2() {
        let mut gens = Structure7::gens2(10);
        gens.remove(1);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str7gens2.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_pickle_gen7_3() {
        let mut gens = Structure7::gens3(10);
        gens.remove(2);
        let rel = relation_slow_3gens(&gens, 50);
        let ref mut f = File::create("./data/str7gens3.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    #[test]
    fn test_gens7() {
        fn test(
            k: usize,
            prec: usize,
            len: usize,
            f: HmfGen<Sqrt5Mpz>,
            gens: &Vec<HmfGen<Sqrt5Mpz>>,
        ) -> bool {
            let mut forms = forms_generated(k, prec, &gens);
            forms.insert(0, f);
            let rels = relations(len, &forms);
            rels.iter().any(|v| !v[0].is_zero_g())
        }
        let prec = 15;
        let gens = Structure7::gens3(prec);
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f20 = rankin_cohen_sqrt5(7, &g2.pow(2), &(&g5.pow(2) * &g6)).unwrap();
        let f21 = rankin_cohen_sqrt5(7, &g5, &(&g2.pow(5) * &g6)).unwrap();
        let gens2 = Structure2::gens(prec);
        let gens5 = Structure5::gens(prec);
        let f17 = &gens2[1] * &gens5[3];
        let f14 = &gens2[1] * &gens5[2];
        assert!(test(14, prec, 100, f14, &gens));
        assert!(test(17, prec, 100, f17, &gens));
        assert!(test(20, prec, 200, f20, &gens));
        assert!(test(21, prec, 200, f21, &gens));
    }

    #[test]
    fn test_gens5_relation12() {
        let prec = 15;
        let gens = Structure5::gens2(prec);
        let forms_with_monom = forms_generated_with_monom(12, prec, &gens);
        let monoms: Vec<_> = forms_with_monom.iter().map(|x| x.1.clone()).collect();
        let forms: Vec<_> = forms_generated(12, prec, &gens);
        println!("{:?}", forms.len());
        let mut rels = relations(100, &forms);
        println!("{:?}", monoms);
        println!("{:?}", rels.len());
        println!("{:?}", rels);
        let rel: Vec<_> = rels.swap_remove(0)
            .into_iter()
            .zip(monoms.into_iter())
            .map(|x| (x.1, x.0))
            .collect();
        let rel = [
            vec![rel[0].clone(), rel[1].clone()],
            vec![rel[2].clone()],
            vec![rel[3].clone()],
            vec![rel[4].clone()],
        ];
        let ref mut f = File::create("./data/str5rel12.sobj").unwrap();
        save_as_pickle_rel(&rel, f);
    }

    fn save_poly_pickle(x: &(PWtPoly, Sqrt5Mpz), f: &mut File) {
        let v: Vec<_> = x.0
            .iter()
            .map(|x| (x.0.idx, Into::<Sqrt5Wrapper>::into(&x.1)))
            .collect();
        let a: Sqrt5Wrapper = From::from(&x.1);
        save_as_pickle((v, a), f);
    }

    fn save_star_norm_as_poly_pickle(f: &HmfGen<Sqrt5Mpz>, len: usize, fl: &mut File) {
        let mut res = HmfGen::new(f.prec);
        star_op(&mut res, f);
        res *= f;
        let wt = res.weight.unwrap();
        assert_eq!(wt.0, wt.1);
        let poly = r_elt_as_pol(&res, len).unwrap();
        save_poly_pickle(&poly, fl);
    }

    #[test]
    fn test_bracket_inner_pol_gens5() {
        let prec = 15;
        let gens = Structure5::gens(prec);
        let f5_6 = &mut File::create("./data/str5br5_6.sobj").unwrap();
        let f5_7 = &mut File::create("./data/str5br5_7.sobj").unwrap();
        let f6_7 = &mut File::create("./data/str5br6_7.sobj").unwrap();
        let f5_10 = &mut File::create("./data/str5br5_10.sobj").unwrap();
        let f6_10 = &mut File::create("./data/str5br6_10.sobj").unwrap();
        let f7_10 = &mut File::create("./data/str5br7_10.sobj").unwrap();
        let pol5_6 = bracket_inner_prod_as_pol(&gens[0], &gens[1], 50).unwrap();
        let pol5_7 = bracket_inner_prod_as_pol(&gens[0], &gens[2], 50).unwrap();
        let pol6_7 = bracket_inner_prod_as_pol(&gens[1], &gens[2], 50).unwrap();
        let pol5_10 = bracket_inner_prod_as_pol(&gens[0], &gens[3], 50).unwrap();
        let pol6_10 = bracket_inner_prod_as_pol(&gens[1], &gens[3], 50).unwrap();
        let pol7_10 = bracket_inner_prod_as_pol(&gens[2], &gens[3], 50).unwrap();
        save_poly_pickle(&pol5_6, f5_6);
        save_poly_pickle(&pol5_7, f5_7);
        save_poly_pickle(&pol6_7, f6_7);
        save_poly_pickle(&pol5_10, f5_10);
        save_poly_pickle(&pol6_10, f6_10);
        save_poly_pickle(&pol7_10, f7_10);
    }

    #[test]
    fn test_bracket_inner_pol_gens7() {
        let prec = 15;
        let gens = Structure7::gens(prec);
        let f = &mut File::create("./data/str7br5_7.sobj").unwrap();
        let f1 = &mut File::create("./data/str7br5_6.sobj").unwrap();
        let pol = bracket_inner_prod_as_pol(&gens[1], &gens[2], 50).unwrap();
        let pol1 = bracket_inner_prod_as_pol(&gens[0], &gens[1], 50).unwrap();
        save_poly_pickle(&pol, f);
        save_poly_pickle(&pol1, f1);
    }

    #[test]
    fn test_pickle_gen8_1() {
        let prec = 15;
        let gens1 = Structure8::gens1(prec);
        let rel1 = relation_slow_3gens(&gens1, 50);
        let ref mut f1 = File::create("./data/str8gens1.sobj").unwrap();
        save_as_pickle_rel(&rel1, f1);
    }

    #[test]
    fn test_bracket_inner_pol_gens_cand8() {
        let prec = 15;
        let gens = Structure8::gens1(prec);
        for (&(i, j), path) in vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .iter()
            .zip(
                &[
                    "./data/str8br6_7.sobj",
                    "./data/str8br6_8_1.sobj",
                    "./data/str8br6_8_2.sobj",
                    "./data/str8br7_8_1.sobj",
                    "./data/str8br7_8_2.sobj",
                    "./data/str8br8_8.sobj",
                ],
            )
        {
            let f = &mut File::create(path).unwrap();
            let pol = bracket_inner_prod_as_pol(&gens[i], &gens[j], 50).unwrap();
            save_poly_pickle(&pol, f);
        }
    }

    #[test]
    fn test_bracket_innter_pol_cand9() {
        let prec = 15;
        let gens = Structure9::gens1(prec);
        for (&(i, j), path) in vec![(0, 1), (0, 2), (1, 2)].iter().zip(
            &[
                "./data/str9br7_8.sobj",
                "./data/str9br7_11.sobj",
                "./data/str9br8_11.sobj",
            ],
        )
        {
            let f = &mut File::create(path).unwrap();
            let pol = bracket_inner_prod_as_pol(&gens[i], &gens[j], 100).unwrap();
            save_poly_pickle(&pol, f);
        }
    }

    #[test]
    fn test_gens9() {
        let prec = 15;
        let gens = Structure9::gens1(prec as u64);
        let g2 = &eisenstein_series(2, prec);
        let g5 = &g5_normalized(prec);
        let g6 = &f6_normalized(prec);
        let f7 = &gens[0];
        let f8 = &gens[1];
        let into = Into::<HmfGen<Sqrt5Mpz>>::into;
        let f24 = into(
            g2.pow(3) * g5 * g6 * &Mpz::from_str_radix("7157983644", 10).unwrap() +
                g5 * g6.pow(2) * &Mpz::from_str_radix("-26483628326400", 10).unwrap() +
                g2 * g5.pow(3) * &Mpz::from_str_radix("-2559782736000", 10).unwrap(),
        ) * f7 +
            into(
                g2.pow(5) * g6 * (703786) + g2.pow(3) * g5.pow(2) * (-1693956600) +
                    g2.pow(8) * (-323) +
                    g2.pow(2) * g6.pow(2) * &Mpz::from_str_radix("18801456480", 10).unwrap(),
            ) * f8;
        let f23 = into(
            g2.pow(5) * g6 * (153) + g2.pow(3) * g5.pow(2) * (-293580) +
                g5.pow(2) * g6 * (-2354184000) +
                g2.pow(2) * g6.pow(2) * (641898),
        ) * f7 +
            into(g2.pow(2) * g5 * g6 * (864000) + g2.pow(5) * g5 * (170)) * f8;
        let f22 = into(
            g2.pow(2) * g5 * g6 * &Mpz::from_str_radix("-13245444000", 10).unwrap() +
                g5.pow(3) * &Mpz::from_str_radix("8475062400000", 10).unwrap() +
                g2.pow(5) * g5 * (2404602),
        ) * f7 +
            into(
                g2.pow(4) * g6 * (2841210) +
                    g2.pow(2) * g5.pow(2) * &Mpz::from_str_radix("5608440000", 10).unwrap() +
                    g2.pow(7) * (-323),
            ) * f8;
        let f21 = into(
            g2.pow(4) * g6 * (-113682) +
                g2 * g6.pow(2) * &Mpz::from_str_radix("-6162220800", 10).unwrap() +
                g2.pow(2) * g5.pow(2) * (542808000) + g2.pow(7) * (323),
        ) * f7 +
            into(
                g2 * g5 * g6 * &Mpz::from_str_radix("7750080000", 10).unwrap() +
                    g2.pow(4) * g5 * (1524900),
            ) * f8;

        for (f, &pth) in vec![f21, f22, f23, f24].iter().zip(
            &[
                "./data/str9gens_norm_21.sobj",
                "./data/str9gens_norm_22.sobj",
                "./data/str9gens_norm_23.sobj",
                "./data/str9gens_norm_24.sobj",
            ],
        )
        {
            let ref mut fl = File::create(pth).unwrap();
            save_star_norm_as_poly_pickle(f, 250, fl);
        }
    }

    #[test]
    fn test_gens10() {
        let prec = 15;
        let gens = Structure10::gens1(prec as u64);
        let g2 = &eisenstein_series(2, prec);
        let g5 = &g5_normalized(prec);
        let g6 = &f6_normalized(prec);
        let f4 = &gens[0];
        let f7 = &gens[1];
        let into = Into::<HmfGen<Sqrt5Mpz>>::into;
        let f21 = into(
            g2.pow(3) * g5 * g6 * (1798193397) +
                g5 * g6.pow(2) * &Mpz::from_str_radix("9017763955200", 10).unwrap() +
                g2 * g5.pow(3) * &Mpz::from_str_radix("-3754266516000", 10).unwrap(),
        ) * f4 +
            into(
                g2.pow(4) * g6 * (117588672) +
                    g2 * g6.pow(2) * &Mpz::from_str_radix("5189042995200", 10).unwrap() +
                    g2.pow(2) * g5.pow(2) * &Mpz::from_str_radix("-1920796416000", 10).unwrap() +
                    g2.pow(7) * (-59081),
            ) * f7;
        let f22 = into(
            g2.pow(6) * g6 * (77) + g2.pow(3) * g6.pow(2) * (-767004) +
                g2 * g5.pow(2) * g6 * (-42192000) +
                g2.pow(4) * g5.pow(2) * (-396880),
        ) * f4 +
            into(
                g2.pow(2) * g5 * g6 * (-276480000) + g2.pow(5) * g5 * (17600),
            ) * f7;
        let h21 = into(
            g5 * g6.pow(2) * &Mpz::from_str_radix("-1874672633433600000", 10).unwrap() +
                g2.pow(6) * g5 * &Mpz::from_str_radix("24575309759", 10).unwrap() +
                g2 * g5.pow(3) * &Mpz::from_str_radix("789630086044800000", 10).unwrap(),
        ) * f4 +
            into(
                g2.pow(4) * g6 * &Mpz::from_str_radix("-19200158784000", 10).unwrap() +
                    g2 * g6.pow(2) * &Mpz::from_str_radix("-1078732704153600000", 10).unwrap() +
                    g2.pow(2) * g5.pow(2) *
                        &Mpz::from_str_radix("459386484480000000", 10).unwrap() +
                    g2.pow(7) * &Mpz::from_str_radix("12904381820", 10).unwrap(),
            ) * f7;
        let f18 = into(
            g2.pow(4) * g6 * (-690852) +
                g2 * g6.pow(2) * &Mpz::from_str_radix("-37868083200", 10).unwrap() +
                g2.pow(2) * g5.pow(2) * &Mpz::from_str_radix("15933456000", 10).unwrap() +
                g2.pow(7) * (451),
        ) * f4;

        for (f, &pth) in vec![f21, f22, h21, f18].iter().zip(
            &[
                "./data/str10gens_norm_21.sobj",
                "./data/str10gens_norm_22.sobj",
                "./data/str10gens_norm_21_1.sobj",
                "./data/str10gens_norm_18.sobj",
            ],
        )
        {
            let ref mut fl = File::create(pth).unwrap();
            save_star_norm_as_poly_pickle(f, 250, fl);
        }
    }

    #[test]
    fn test_gens5_relation11() {
        let prec = 15;
        let gens = Structure5::gens2(prec);
        let forms_with_monom = forms_generated_with_monom(11, prec, &gens);
        let monoms: Vec<_> = forms_with_monom.iter().map(|x| x.1.idx).collect();
        let forms: Vec<_> = forms_generated(11, prec, &gens);
        println!("{:?}", forms.len());
        let rels = relations(100, &forms);
        println!("{:?}", monoms);
        println!("{:?}", rels.len());
        println!("{:?}", rels);
    }

    // #[test]
    // fn test_gens5_2() {
    //     fn rank(k: usize, prec: usize, len: usize, gens: &Vec<HmfGen<Sqrt5Mpz>>) -> usize {
    //         let forms: Vec<_> = forms_generated(k, prec, gens);
    //         let rels = relations(len, &forms);
    //         forms.len() - rels.len()
    //     }
    //     let prec = 50;
    //     let gens = Structure5::gens2(prec);
    //     // [1, 1, 2, 1, 2, 3, 3, 4, 4, 4, 6, 6, 7, 7, 8, 9, 10, 11, 11, 12, 14, 14, 16, 16, 17, 19]
    //     let v: Vec<_> = (5..31).map(|i| rank(i, prec, 800, &gens)).collect();
    //     println!("{:?}", v);
    // }

    #[test]
    fn test_gens5_3() {
        let prec = 15;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f12 = rankin_cohen_sqrt5(5, &(&g2 * &g5), &g5).unwrap();

        let gens = Structure5::gens(prec);
        let mut forms = forms_generated(12, prec, &gens);
        forms.insert(0, f12);
        let rels = relations(100, &forms);
        assert!(!rels[0][0].is_zero_g());

        let f15 = rankin_cohen_sqrt5(5, &(&g2 * &g5), &(&g6 * &g2)).unwrap();
        let mut forms = forms_generated(15, prec, &gens);
        forms.insert(0, f15);
        let rels = relations(100, &forms);
        assert!(!rels[0][0].is_zero_g());
    }

    #[test]
    fn test_gens6() {
        let prec = 15;
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f11 = rankin_cohen_sqrt5(6, &g5, &g6).unwrap();
        let mut forms = forms_generated(11, prec, &Structure6::gens(prec));
        for f in forms.iter() {
            assert_eq!(f.weight, Some((11, 23)));
        }
        // let rel = relation(50, &f11, &forms);
        // println!("{:?}", rel);
        forms.insert(0, f11);
        let nums = vec![13440, 366, 0, 17820, 80];
        let mut s = HmfGen::new(prec);
        let mut tmp = HmfGen::new(prec);
        for (&ref f, &ref a) in forms.iter().zip(nums.iter()) {
            tmp.mul_mut_by_const(&f, &Sqrt5Mpz::from_si_g(*a));
            s += &tmp;
        }
        assert!(s.is_zero());
    }

    #[allow(dead_code)]
    fn forms_rel_to_rel(k: usize, gens: &[HmfGen<Sqrt5Mpz>], rel: &Vec<Sqrt5Mpz>) -> Relation {
        let monomss = forms_generated_monom(k, &gens);
        let mut rel = rel.clone();
        let mut res = Vec::new();
        for ms in &monomss {
            let mut v = Vec::new();
            for m in ms.iter() {
                let a = rel.remove(0);
                if !a.is_zero_g() {
                    v.push((m.clone(), a));
                }
            }
            res.push(v);
        }
        res
    }

    #[test]
    fn test_gens7_rel() {
        let prec = 15;
        let gens = Structure7::gens3(prec);
        for &k in &[10, 13] {
            let forms = forms_generated(k, prec, &gens);
            let rels = relations(100, &forms);
            let rel = forms_rel_to_rel(k, &gens, &rels[0]);
            if !(rels.is_empty()) {
                println!("{}", k);
                println!("{:?}", rel);
                let ref mut f = File::create(format!("./data/str7rel{}.sobj", k)).unwrap();
                save_as_pickle_rel(&rel, f);
            }
        }
    }

    #[test]
    fn test_g8_gens_construction() {
        let prec = 15;
        let gens = Structure8::gens1(prec);
        let g2 = &eisenstein_series(2, prec);
        let g5 = &g5_normalized(prec);
        let g6 = &f6_normalized(prec);
        let f6 = &gens[0];
        let f7 = &gens[1];
        let into = Into::<HmfGen<Sqrt5Mpz>>::into;
        let f17 = into(g2.pow(3) * g5 * 47 + g5 * g6 * 86400) * f6 +
            into(g2.pow(2) * g6 * (-1890) + g5.pow(2) * 3240000) * f7;
        let f18 = into(g2.pow(3) * g6 * 47 + g6.pow(2) * 86400) * f6 +
            into(g2.pow(3) * g5 * (-1050) + g5 * g6 * 3240000) * f7;

        for (f, &pth) in vec![f17, f18].iter().zip(
            &[
                "./data/str8gens_norm_17.sobj",
                "./data/str8gens_norm_18.sobj",
            ],
        )
        {
            let ref mut fl = File::create(pth).unwrap();
            save_star_norm_as_poly_pickle(f, 250, fl);
        }
    }

    #[test]
    fn test_bracket_innter_pol_cand10() {
        let prec = 15;
        let gens = Structure10::gens1(prec);
        for (&(i, j), path) in vec![(0, 1), (0, 2), (1, 2)].iter().zip(
            &[
                "./data/str10br4_7.sobj",
                "./data/str10br4_11.sobj",
                "./data/str10br7_11.sobj",
            ],
        )
        {
            let f = &mut File::create(path).unwrap();
            let pol = bracket_inner_prod_as_pol(&gens[i], &gens[j], 100).unwrap();
            save_poly_pickle(&pol, f);
        }
    }
    // #[test]
    // fn test_gens7_rel1() {
    //     let prec = 15;
    //     let gens = Structure7::gens(prec);
    //     for k in 5..30 {
    //         let forms = forms_generated(k, prec, &gens);
    //         let rels = relations(200, &forms);
    //         let dim = forms.len() - rels.len();
    //         // [1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 8, 10, 10, 11, 12, 13,
    //         // 14, 15, 16, 17, 18, 20]
    //         println!("{}: {}", k, dim);
    //     }
    // }

    #[test]
    fn test_save_rels() {
        let prec = 15;
        let ref mut f = File::create("./data/rels.sobj").unwrap();
        let v: Vec<_> = (11..20)
            .map(|i| (i, three_forms_rel(i, prec, 50)))
            .take_while(|x| x.1.is_some())
            .map(|x| (x.0, x.1.unwrap()))
            .collect();
        save_as_pickle_3relations(&v, f);
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
    fn test_save_rels6() {
        let prec = 15;
        let gens1 = Structure6::gens1(prec);
        let v = brackets(&gens1);
        let ref mut f = File::create("./data/str6_brs.sobj").unwrap();
        save_polys_over_z_pickle(&v, f);
    }

    #[test]
    fn test_save_rels5() {
        let prec = 15;
        let gens1 = mixed_weight_forms(5, prec, 5);
        let gens = [gens1[0].0.clone(), gens1[1].0.clone(), gens1[4].0.clone()];
        for f in &gens {
            print!("{}, ", f.weight.unwrap().0);
        }
        println!("\n");
        let v = brackets(&gens);
        let ref mut f = File::create("./data/str5_brs.sobj").unwrap();
        save_polys_over_z_pickle(&v, f);
    }

    #[test]
    fn test_save_rels_up_to_50() {
        for i in 3..21 {
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


    fn save_star_norms(diffs: &[u64]) {

        for &i in diffs {
            let cand = {
                let cand_f = File::open(format!("./data/brackets/str{}_cand.sobj", i)).unwrap();
                let monom_f = File::open(format!("./data/brackets/str{}_monoms.sobj", i)).unwrap();
                StrCand::load(i, &cand_f, &monom_f).unwrap()
            };
            let wt_mx = cand.gens_nums_wts()
                .iter()
                .map(|x| (x.0 + x.1) / 5)
                .max()
                .unwrap();
            let prec = wt_mx as usize;
            let gens = cand.gens_nums_as_forms(prec);
            println!("i: {}, prec: {}", i, prec);
            let stars_f = format!("./data/brackets/str{}_star_norms.sobj", i);
            measure_time!(cand.save_star_norms(&gens, &stars_f));
        }
    }

    #[test]
    fn test_save_star_norms0() {
        save_star_norms(&(36..51).filter(|i| i % 2 == 0).collect::<Vec<_>>());
    }
    #[test]
    fn test_save_star_norms1() {
        save_star_norms(&(36..51).filter(|i| i % 2 == 1).collect::<Vec<_>>());
    }

    fn write_csv_form(f: &HmfGen<Sqrt5Mpz>, p: &String) {
        let mut wtr = csv::Writer::from_path(p).unwrap();
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
            write_csv_form(f, &p.to_string());
        }
    }

    #[test]
    fn write_gens_cands() {
        for i in 1..51 {
            let cand = {
                let cand_f = File::open(format!("./data/brackets/str{}_cand.sobj", i)).unwrap();
                let monom_f = File::open(format!("./data/brackets/str{}_monoms.sobj", i)).unwrap();
                StrCand::load(i, &cand_f, &monom_f).unwrap()
            };
            let prec = 15;
            let gens = cand.gens(prec);
            println!("{}", i);
            for (n, f) in gens.iter().enumerate() {
                let (w1, w2) = f.weight.unwrap();
                let path = format!("./forms_csv/A{}/gen{}_wt_{}_{}.csv", i, n, w1, w2);
                write_csv_form(f, &path);
            }
        }
    }
}

mod str_test {
    use hilbert_sqrt5::structure::*;
    use hilbert_sqrt5::eisenstein::eisenstein_series;

    #[test]
    fn test_relations_over_z() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let f = &g2 * 2;
        let forms = [f, g2];
        let rels = relations_over_z(&forms);
        assert_eq!(rels.len(), 1);
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

use std::time::Instant;
use theta_chars::theta;
use elements::HmfGen;
use eisenstein::eisenstein_series;
use misc::prime_sieve;
use gmp::mpz::Mpz;


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

mod theta_eisen_relatioin {
    use super::*;

    #[test]
    fn test_relation() {
        let prec = 10;
        let mut g5 = theta(prec);
        assert!(g5.is_divisible_by_const(&Mpz::from_ui(64)));
        g5 /= &Mpz::from_ui(64);
        println!("{}", g5);
        let e2 = eisenstein_series(2, prec);
        let e6 = eisenstein_series(6, prec);
        let e10 = eisenstein_series(10, prec);
        let mut tmp = HmfGen::new(prec);
        tmp.pow_mut(&e2, 5);
        let f10_1 = tmp.clone();
        tmp.pow_mut(&e2, 2);
        let f10_2 = &tmp * &e6;
        let g5_sq = &g5 * &g5;
        let mut g10 = &(&e10 + &(&f10_1 * &Mpz::from_ui(355404))) + &(&f10_2 * &Mpz::from_si(-11465));
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
}

mod theta {
    use super::*;

    #[test]
    fn theta_fun() {
        measure_time!(theta(10));
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

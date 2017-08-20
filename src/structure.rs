use gmp::mpz::Mpz;
use eisenstein::{eisenstein_series, f6_normalized};
use theta_chars::g5_normalized;
use elements::{HmfGen, div_mut};
use serde_pickle;
use std::fs::File;
use std::io::Write;
use std::io::Read;
use bignum::{RealQuadElement, BigNumber};
use serde;
use std::process::Command;
use std::env;
use bignum::Sqrt5Mpz;
use diff_op::{rankin_cohen_sqrt5, bracket_inner_prod1};
use misc::PowGen;


/// A stupid function that returns a linear relation.
pub fn relation(len: usize, f: &HmfGen<Sqrt5Mpz>, forms: &Vec<HmfGen<Mpz>>) -> Vec<Sqrt5Mpz> {
    let vv: Vec<Vec<Mpz>> = forms.iter().map(|f| f.fc_vector(len)).collect();
    let v = f.fc_vector(len);
    let path_name = "./data/rust_python_data.sobj";
    let res_path_name = "./data/rust_python_data_res.sobj";
    let mut f = File::create(path_name).unwrap();
    save_as_pickle_quadz_vec(&v, &vv, &mut f);
    let corank = sage_command("print load('./src/relation.sage')")
        .lines()
        .next()
        .unwrap()
        .to_string();
    assert_eq!(corank.parse::<i32>().unwrap(), 1);
    let f = File::open(res_path_name).unwrap();
    load_pickle_quadz(&f).unwrap()
}

type PWtPoly = Vec<(MonomFormal, Sqrt5Mpz)>;
type Relation = Vec<PWtPoly>;

pub fn relation_monom(len: usize, f: &HmfGen<Sqrt5Mpz>) -> (Sqrt5Mpz, PWtPoly) {
    let wt = f.weight.unwrap();
    assert_eq!(wt.0, wt.1);
    let forms_monom = monoms_of_g2_g5_f6(wt.0);
    let forms: Vec<_> = forms_monom.iter().map(|x| x.into_form(f.prec)).collect();
    let mut v = relation(len, &f, &forms);
    let a = v.remove(0);
    let vec = forms_monom.into_iter().zip(v.into_iter()).collect();
    (a, vec)
}

fn sage_command(cmd: &'static str) -> String {
    let home = env::home_dir().unwrap();
    let output = Command::new(home.join("bin/sage"))
        .arg("-c")
        .arg(cmd)
        .output()
        .expect("failed to execute process");

    if !output.status.success() {
        panic!("stderr: {}", String::from_utf8_lossy(&output.stderr));
    }
    String::from_utf8_lossy(&output.stdout).into()
}

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

fn monom_g2_g5_f6(prec: usize, expts: (usize, usize, usize)) -> HmfGen<Mpz> {
    let mut tmp = HmfGen::new(prec);
    let mut res = HmfGen::new(prec);
    let g2 = eisenstein_series(2, prec);
    let g5 = g5_normalized(prec);
    let g6 = f6_normalized(prec);
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

fn save_as_pickle_quadz_vec<T>(vec: &Vec<T>, basis_vec: &Vec<Vec<Mpz>>, f: &mut File)
where
    T: RealQuadElement<Mpz>,
{
    let v: Vec<Vec<String>> = basis_vec
        .iter()
        .map(|v| v.iter().map(|x| x.to_str_radix(10)).collect())
        .collect();
    let w: Vec<(String, String)> = vec.iter()
        .map(|a| {
            (a.rt_part().to_str_radix(10), a.ir_part().to_str_radix(10))
        })
        .collect();
    save_as_pickle(&(w, v), f);
}


#[allow(dead_code)]
fn save_as_pickle_quadz<T>(vec: &Vec<T>, f: &mut File)
where
    T: RealQuadElement<Mpz>,
{
    let v: Vec<(String, String)> = vec.iter()
        .map(|x| (x.rt_part(), x.ir_part()))
        .map(|(x, y)| (x.to_str_radix(10), y.to_str_radix(10)))
        .collect();
    save_as_pickle(&v, f);
}

#[allow(dead_code)]
fn load_pickle_quadz<T>(f: &File) -> Result<Vec<T>, serde_pickle::Error>
where
    T: RealQuadElement<Mpz>,
    for<'a> T: From<&'a (Mpz, Mpz)>,
{
    let v: Vec<(String, String)> = try!(load_pickle(f));
    let res = v.iter()
        .map(|&(ref x, ref y)| {
            (
                Mpz::from_str_radix(x, 10).unwrap(),
                Mpz::from_str_radix(y, 10).unwrap(),
            )
        })
        .map(|t| From::from(&t))
        .collect();
    Ok(res)
}

#[allow(dead_code)]
fn save_as_pickle_z(vec: &Vec<Mpz>, f: &mut File) {
    let vec: Vec<String> = vec.iter().map(|x| x.to_str_radix(10)).collect();
    save_as_pickle(&vec, f);
}

#[allow(dead_code)]
fn load_pickle_z(f: &File) -> Result<Vec<Mpz>, serde_pickle::Error> {
    let v: Vec<String> = load_pickle(&f)?;
    let res = v.iter()
        .map(|x| Mpz::from_str_radix(x, 10).unwrap())
        .collect();
    Ok(res)
}

#[allow(dead_code)]
fn save_as_pickle<T>(a: T, f: &mut File)
where
    T: serde::Serialize,
{
    let v = serde_pickle::to_vec(&a, false).unwrap();
    f.write(&v).unwrap();
}

#[allow(dead_code)]
fn load_pickle<'de, T>(f: &File) -> Result<Vec<T>, serde_pickle::Error>
where
    T: serde::Deserialize<'de>,
{
    let buf: Vec<u8> = f.bytes().map(|x| x.unwrap()).collect();
    serde_pickle::from_slice(&buf)
}

/// Corresponds to g2^a * g5^b * f6^c where (a, b, c) = idx.
#[derive(Debug)]
pub struct MonomFormal {
    pub idx: (usize, usize, usize),
}

impl MonomFormal {
    pub fn into_form(&self, prec: usize) -> HmfGen<Mpz> {
        monom_g2_g5_f6(prec, self.idx)
    }
    pub fn eval(v: &PWtPoly, prec: usize) -> HmfGen<Sqrt5Mpz> {
        let mut res = HmfGen::new(prec);
        res.set(&v[0].0.into_form(prec));
        let mut res = From::from(&res);
        res *= &v[0].1;
        for &(ref monom, ref a) in v.iter().skip(1) {
            let mut tmp: HmfGen<Sqrt5Mpz> = From::from(&monom.into_form(prec));
            tmp *= a;
            res += &tmp;
        }
        res
    }
}


pub fn monoms_of_g2_g5_f6(k: usize) -> Vec<MonomFormal> {
    tpls_of_wt(k)
        .iter()
        .map(|&x| MonomFormal { idx: x })
        .collect()
}

fn eval_relation(rel: &Relation, gens: &Vec<HmfGen<Sqrt5Mpz>>) -> HmfGen<Sqrt5Mpz> {
    let prec = gens[0].prec;
    let mut res = HmfGen::<Sqrt5Mpz>::new(prec);
    for (&ref p, f) in rel.iter().zip(gens.iter()) {
        let mut tmp = MonomFormal::eval(p, prec);
        tmp *= f;
        res += &tmp;
    }
    res
}

pub trait Structure {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>>;
    fn relations() -> Option<Vec<Relation>>;

    fn check_relations(prec: usize) -> bool {
        let gens = Self::gens(prec);
        if gens.len() == 2 {
            return true;
        } else {
            let rels = Self::relations().unwrap();
            for rel in rels.iter() {
                if !eval_relation(rel, &gens).is_zero() {
                    return false;
                }
            }
            true
        }
    }
}

pub struct Structure1;

impl Structure for Structure1 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let g7_9 = rankin_cohen_sqrt5(1, &g2, &g5).unwrap();
        let g8_10 = rankin_cohen_sqrt5(1, &g2, &g6).unwrap();
        let g11_13 = rankin_cohen_sqrt5(1, &g5, &g6).unwrap();
        assert_eq!(g7_9.weight, Some((7, 9)));
        assert_eq!(g8_10.weight, Some((8, 10)));
        assert_eq!(g11_13.weight, Some((11, 13)));
        vec![g7_9, g8_10, g11_13]
    }

    fn relations() -> Option<Vec<Relation>> {
        let a0 = (MonomFormal { idx: (0, 0, 1) }, Sqrt5Mpz::from_si_g(-6));
        let a1 = (MonomFormal { idx: (0, 1, 0) }, Sqrt5Mpz::from_si_g(5));
        let a2 = (MonomFormal { idx: (1, 0, 0) }, Sqrt5Mpz::from_si_g(-2));
        Some(vec![vec![vec![a0], vec![a1], vec![a2]]])
    }
}

impl Structure1 {
    #[allow(dead_code)]
    fn relation_slow() {
        let prec = 10;
        let gens = Self::gens(prec);
        let g7_9 = &gens[0];
        let g8_10 = &gens[1];
        let g11_13 = &gens[2];
        let h6 = bracket_inner_prod1(g8_10, g11_13).unwrap();
        let h5 = bracket_inner_prod1(g11_13, g7_9).unwrap();
        let h2 = bracket_inner_prod1(g7_9, g8_10).unwrap();
        let (a6, mut v6) = relation_monom(50, &h6);
        let (a5, mut v5) = relation_monom(50, &h5);
        let (a2, mut v2) = relation_monom(50, &h2);
        let ref mut tmp = Mpz::new();
        v6[1].1.mul_assign_g(&a5, tmp);
        v6[1].1.mul_assign_g(&a2, tmp);
        v5[0].1.mul_assign_g(&a6, tmp);
        v5[0].1.mul_assign_g(&a2, tmp);
        v2[0].1.mul_assign_g(&a6, tmp);
        v2[0].1.mul_assign_g(&a5, tmp);
        println!("{:?}", v6[1]);
        println!("{:?}", v5[0]);
        println!("{:?}", v2[0]);
    }
}

pub struct Structure2;

impl Structure for Structure2 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g4_8 = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        let g7_11 = rankin_cohen_sqrt5(2, &g2, &g5).unwrap();
        assert_eq!(g4_8.weight, Some((4, 8)));
        assert_eq!(g7_11.weight, Some((7, 11)));
        vec![g4_8, g7_11]
    }

    fn relations() -> Option<Vec<Relation>> {
        None
    }
}

pub struct Structure4;

impl Structure for Structure4 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g4_12 = rankin_cohen_sqrt5(4, &g2, &g2).unwrap();
        let g7_15 = rankin_cohen_sqrt5(4, &g2, &g5).unwrap();
        let mut g8_16 = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        let g2_sqrt5: HmfGen<Sqrt5Mpz> = From::from(&g2);
        g8_16.square();
        let mut g5_13 = HmfGen::new(prec);
        div_mut(&mut g5_13, &g7_15, &g2_sqrt5);
        assert_eq!(g4_12.weight, Some((4, 12)));
        assert_eq!(g5_13.weight, Some((5, 13)));
        assert_eq!(g8_16.weight, Some((8, 16)));
        vec![g4_12, g5_13, g8_16]
    }

    fn relations() -> Option<Vec<Relation>> {
        let a0 = (MonomFormal { idx: (0, 0, 1) }, Sqrt5Mpz::from_si_g(1944));
        let a1 = (MonomFormal { idx: (0, 1, 0) }, Sqrt5Mpz::from_si_g(259200));
        let a2 = (MonomFormal { idx: (1, 0, 0) }, Sqrt5Mpz::from_si_g(-7));
        Some(vec![vec![vec![a0], vec![a1], vec![a2]]])
    }
}


#[allow(dead_code)]
fn relation_slow_3gens(gens: &Vec<HmfGen<Sqrt5Mpz>>, len: usize) {
    let g0 = &gens[0];
    let g1 = &gens[1];
    let g2 = &gens[2];
    let h0 = bracket_inner_prod1(g1, g2).unwrap();
    let h1 = bracket_inner_prod1(g2, g0).unwrap();
    let h2 = bracket_inner_prod1(g0, g1).unwrap();
    let (a0, mut v0) = relation_monom(len, &h0);
    let (a1, mut v1) = relation_monom(len, &h1);
    let (a2, mut v2) = relation_monom(len, &h2);
    let ref mut tmp = Mpz::new();
    for &mut (_, ref mut b) in v0.iter_mut() {
        b.mul_assign_g(&a1, tmp);
        b.mul_assign_g(&a2, tmp);
    }
    for &mut (_, ref mut b) in v1.iter_mut() {
        b.mul_assign_g(&a0, tmp);
        b.mul_assign_g(&a2, tmp);
    }
    for &mut (_, ref mut b) in v2.iter_mut() {
        b.mul_assign_g(&a0, tmp);
        b.mul_assign_g(&a1, tmp);
    }
    println!("{:?}", &v0);
    println!("{:?}", &v1);
    println!("{:?}", &v2);
}

impl Structure4 {
    #[allow(dead_code)]
    fn relation_slow() {
        let prec = 10;
        let gens = Self::gens(prec);
        relation_slow_3gens(&gens, 50);
    }
}

pub struct Structure3;

impl Structure for Structure3 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let g7_13 = rankin_cohen_sqrt5(3, &g2, &g5).unwrap();
        let g8_14 = rankin_cohen_sqrt5(3, &g2, &g6).unwrap();
        let g11_17 = rankin_cohen_sqrt5(3, &g5, &g6).unwrap();
        vec![g7_13, g8_14, g11_17]
    }

    fn relations() -> Option<Vec<Relation>> {
        let a0_0 = (MonomFormal { idx: (0, 2, 0) }, Sqrt5Mpz::from_si_g(-5880));
        let a0_1 = (MonomFormal { idx: (2, 0, 1) }, Sqrt5Mpz::from_si_g(7));
        let a1 = (MonomFormal { idx: (2, 1, 0) }, Sqrt5Mpz::from_si_g(7));
        let a2 = (MonomFormal { idx: (3, 0, 0) }, Sqrt5Mpz::from_si_g(-8));
        Some(vec![vec![vec![a0_0, a0_1], vec![a1], vec![a2]]])
    }
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
        let mut f = File::create("/tmp/foo.sobj").unwrap();
        let g = File::open("/tmp/foo.sobj").unwrap();
        save_as_pickle_z(&v, &mut f);
        let w = load_pickle_z(&g).unwrap();
        assert_eq!(v, w);

        let mut f1 = File::create("/tmp/bar.sobj").unwrap();
        let v: Vec<Sqrt5Mpz> = (0..10000)
            .map(|x| From::from(&(Mpz::from_ui(x), Mpz::from_ui(x))))
            .collect();
        save_as_pickle_quadz(&v, &mut f1);
        let g1 = File::open("/tmp/bar.sobj").unwrap();
        let w: Vec<Sqrt5Mpz> = load_pickle_quadz(&g1).unwrap();
        assert_eq!(v, w);

        let mut f = File::create("/tmp/foo.sobj").unwrap();
        let v = vec![
            vec![Mpz::from_si(2), Mpz::from_si(-2)],
            vec![Mpz::from_si(5), Mpz::from_si(-3)],
        ];
        let w = vec![Sqrt5Mpz::from_sisi(2, 4), Sqrt5Mpz::from_sisi(3, 5)];
        save_as_pickle_quadz_vec(&w, &v, &mut f);
    }

    #[test]
    fn relation_slow1() {
        println!("{:?}", Structure1::relations());
        let gens = Structure1::gens(10);
        relation_slow_3gens(&gens, 50);
    }

    #[test]
    fn check_relations1() {
        assert!(Structure1::check_relations(10));
    }

    #[test]
    fn relation_slow4() {
        let gens = Structure4::gens(10);
        relation_slow_3gens(&gens, 50);
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
    fn relation_slow3() {
        let gens = Structure3::gens(10);
        relation_slow_3gens(&gens, 50);
    }
}

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


pub fn rank<T>(len: usize, forms: &Vec<HmfGen<T>>) -> usize
where
    T: RealQuadElement<Mpz> + BigNumber,
{
    let vv: Vec<Vec<T>> = forms.iter().map(|f| f.fc_vector(len)).collect();
    let mut v: Vec<T> = Vec::new();
    for _ in 0..len {
        v.push(T::new_g());
    }
    let path_name = "./data/rust_python_data.sobj";
    let mut f = File::create(path_name).unwrap();
    save_as_pickle_quadz_vec(&v, &vv, &mut f);
    let rank = sage_command("print load('./src/relation.sage')")
        .lines()
        .next()
        .unwrap()
        .to_string();
    let rank = rank.parse::<usize>().unwrap();
    rank
}

/// A stupid function that returns a linear relation.
pub fn relation<T>(len: usize, f: &HmfGen<T>, forms: &Vec<HmfGen<T>>) -> Vec<T>
where
    T: RealQuadElement<Mpz> + BigNumber,
    for<'a> T: From<&'a (Mpz, Mpz)>,
{
    let vv: Vec<Vec<T>> = forms.iter().map(|f| f.fc_vector(len)).collect();
    let v = f.fc_vector(len);
    let path_name = "./data/rust_python_data.sobj";
    let res_path_name = "./data/rust_python_data_res.sobj";
    let mut f = File::create(path_name).unwrap();
    save_as_pickle_quadz_vec(&v, &vv, &mut f);
    let rank = sage_command("print load('./src/relation.sage')")
        .lines()
        .next()
        .unwrap()
        .to_string();
    assert_eq!(rank.parse::<usize>().unwrap(), forms.len());
    let f = File::open(res_path_name).unwrap();
    load_pickle_quadz(&f).unwrap()
}

type PWtPoly = Vec<(MonomFormal, Sqrt5Mpz)>;
type Relation = Vec<PWtPoly>;

pub fn relation_monom(len: usize, f: &HmfGen<Sqrt5Mpz>) -> (Sqrt5Mpz, PWtPoly) {
    let wt = f.weight.unwrap();
    assert_eq!(wt.0, wt.1);
    let forms_monom = monoms_of_g2_g5_f6(wt.0);
    let forms: Vec<_> = forms_monom
        .iter()
        .map(|x| From::from(&x.into_form(f.prec)))
        .collect();
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

fn save_as_pickle_quadz_vec<T>(vec: &Vec<T>, basis_vec: &Vec<Vec<T>>, f: &mut File)
where
    T: RealQuadElement<Mpz>,
{
    let v: Vec<Vec<(String, String)>> = basis_vec
        .iter()
        .map(|v| {
            v.iter()
                .map(|x| {
                    (x.rt_part().to_str_radix(10), x.ir_part().to_str_radix(10))
                })
                .collect()
        })
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

type PolString = Vec<((usize, usize, usize), String, String)>;
#[allow(dead_code)]
fn rel3_to_tuple(rel: &(PWtPoly, PWtPoly, PWtPoly)) -> (PolString, PolString, PolString) {
    let &(ref p0, ref p1, ref p2) = rel;
    let to_vec = |p: &PWtPoly| {
        p.iter()
            .map(|&(ref m, ref a)| {
                (
                    m.idx,
                    a.rt_part().to_str_radix(10),
                    a.ir_part().to_str_radix(10),
                )
            })
            .collect()
    };
    let v0: Vec<_> = to_vec(p0);
    let v1: Vec<_> = to_vec(p1);
    let v2: Vec<_> = to_vec(p2);
    (v0, v1, v2)
}

#[allow(dead_code)]
fn save_as_pickle_3relations(rels: &Vec<(usize, (PWtPoly, PWtPoly, PWtPoly))>, f: &mut File) {
    let v: Vec<_> = rels.iter()
        .map(|&(i, ref rel)| (i, rel3_to_tuple(rel)))
        .collect();
    save_as_pickle(v, f);
}

#[allow(dead_code)]
fn save_as_pickle_rel3(rel: &(PWtPoly, PWtPoly, PWtPoly), f: &mut File) {
    let &(ref p0, ref p1, ref p2) = rel;
    let to_vec = |p: &PWtPoly| {
        p.iter()
            .map(|&(ref m, ref a)| {
                (
                    m.idx,
                    a.rt_part().to_str_radix(10),
                    a.ir_part().to_str_radix(10),
                )
            })
            .collect()
    };
    let v0: Vec<_> = to_vec(p0);
    let v1: Vec<_> = to_vec(p1);
    let v2: Vec<_> = to_vec(p2);
    save_as_pickle((v0, v1, v2), f);
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
fn three_forms(i: usize, prec: usize) -> Option<Vec<HmfGen<Sqrt5Mpz>>> {
    let v = if is_even!(i) {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f0 = rankin_cohen_sqrt5(i, &g2, &g2).unwrap();
        let f1 = rankin_cohen_sqrt5(i, &g2, &g5).unwrap();
        let f2 = rankin_cohen_sqrt5(i, &g5, &g6).unwrap();
        vec![f0, f1, f2]
    } else {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f0 = rankin_cohen_sqrt5(i, &g2, &g5).unwrap();
        let f1 = rankin_cohen_sqrt5(i, &g2, &g6).unwrap();
        let f2 = rankin_cohen_sqrt5(i, &g5, &g6).unwrap();
        vec![f0, f1, f2]
    };
    if v.iter().all(|f| f.is_zero()) {
        None
    } else {
        Some(v)
    }
}

#[allow(dead_code)]
fn three_forms_rel(i: usize, prec: usize, len: usize) -> Option<(PWtPoly, PWtPoly, PWtPoly)> {
    let gens = three_forms(i, prec);
    gens.map(|ref x| relation_slow_3gens(x, len))
}

#[allow(dead_code)]
fn relation_slow_3gens(gens: &Vec<HmfGen<Sqrt5Mpz>>, len: usize) -> (PWtPoly, PWtPoly, PWtPoly) {
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
    (v0, v1, v2)
}

fn print_3rel(rel: (PWtPoly, PWtPoly, PWtPoly)) {
    println!("{:?}", rel.0);
    println!("{:?}", rel.1);
    println!("{:?}", rel.2);
}

impl Structure4 {
    #[allow(dead_code)]
    fn relation_slow() {
        let prec = 10;
        let gens = Self::gens(prec);
        print_3rel(relation_slow_3gens(&gens, 50));
    }
}

pub struct Structure3;

impl Structure3 {
    // This is not gens.
    #[allow(dead_code)]
    fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        // let mut h6 = HmfGen::new(prec);
        // h6.pow_mut(&g2, 3);
        let h0 = rankin_cohen_sqrt5(1, &g2, &g5).unwrap();
        let h1 = rankin_cohen_sqrt5(2, &g2, &g2).unwrap();
        let mut f0 = &h0 * &h1;
        let mut f2 = rankin_cohen_sqrt5(3, &g5, &g6).unwrap();
        f0 *= &Sqrt5Mpz::from_si_g(7);
        f0 += &(&f2 * &Sqrt5Mpz::from_si_g(-1080));
        assert_eq!(f0.weight, Some((11, 17)));
        let f1 = rankin_cohen_sqrt5(3, &g2, &g6).unwrap();
        f2 *= &Sqrt5Mpz::from_si_g(720);
        f2 -= &f0;
        vec![f0, f1, f2]
    }
}

impl Structure for Structure3 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let gens = Self::gens1(prec + 1);
        let mut f0 = HmfGen::new(prec + 1);
        let mut f2 = HmfGen::new(prec + 1);
        {
            let f11_17 = &gens[0];
            let g11_17 = &gens[2];
            let g2 = eisenstein_series(2, prec + 1);
            let g5 = g5_normalized(prec + 1);
            let g6 = f6_normalized(prec + 1);
            div_mut(&mut f0, &f11_17, &From::from(&g5));
            div_mut(&mut f2, &g11_17, &From::from(&(&g2 * &g6)));
        }
        assert_eq!(f0.weight, Some((6, 12)));
        assert_eq!(f2.weight, Some((3, 9)));
        f0.decrease_prec(prec);
        f2.decrease_prec(prec);
        vec![f0, f2]
    }

    fn relations() -> Option<Vec<Relation>> {
        None
    }
}

pub struct Structure6;

impl Structure6 {
    #[allow(dead_code)]
    fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let i = 6;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let k4 = rankin_cohen_sqrt5(i, &g2, &g2).unwrap();
        let mut f0 = HmfGen::new(prec);
        div_mut(&mut f0, &k4, &From::from(&g2));
        let f1 = rankin_cohen_sqrt5(i, &g2, &g5).unwrap();
        let mut gens = Structure3::gens(prec);
        let mut f2 = gens.swap_remove(1);
        f2.square();
        assert_eq!(f2.weight, Some((6, 18)));
        vec![f0, f1, f2]
    }

    #[allow(dead_code)]
    fn gens2(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let i = 6;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f0 = rankin_cohen_sqrt5(i, &g5, &g6).unwrap();
        let f1 = rankin_cohen_sqrt5(i, &g2, &g5).unwrap();
        let mut gens = Structure3::gens(prec);
        let mut f2 = gens.swap_remove(1);
        f2.square();
        assert_eq!(f2.weight, Some((6, 18)));
        vec![f0, f1, f2]
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
    fn relation_slow1() {
        println!("{:?}", Structure1::relations());
        let gens = Structure1::gens(10);
        print_3rel(relation_slow_3gens(&gens, 50));
    }

    #[test]
    fn check_relations1() {
        assert!(Structure1::check_relations(10));
    }

    #[test]
    fn relation_slow4() {
        let gens = Structure4::gens(10);
        print_3rel(relation_slow_3gens(&gens, 50));
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
        print_3rel(relation_slow_3gens(&gens, 50));
    }

    #[test]
    fn test_pickle_gen3() {
        let gens1 = Structure3::gens1(10);
        let rel1 = relation_slow_3gens(&gens1, 50);
        let ref mut f1 = File::create("./data/str3gens1.sobj").unwrap();
        save_as_pickle_rel3(&rel1, f1);
        let _gens = Structure3::gens(13);
        // let rel = relation_slow_3gens(&gens, 50);
        // let ref mut f = File::create("./data/str3gens.sobj").unwrap();
        // save_as_pickle_rel3(&rel, f);
    }

    #[test]
    fn test_gens3() {
        let prec = 10;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f = rankin_cohen_sqrt5(3, &g2, &g6).unwrap();
        let gens = Structure3::gens(prec);
        let f0 = &gens[0];
        let f1 = &gens[1];
        let g5: HmfGen<Sqrt5Mpz> = From::from(&g5);
        let g2: HmfGen<Sqrt5Mpz> = From::from(&g2);
        let h0 = f1 * &g5;
        let h1 = f0 * &g2;
        // {
        //     let forms = vec![h0, h1];
        //     println!("{:?}", relation(50, &f, &forms));
        // }

        let mut res = HmfGen::new(prec);
        let mut tmp = HmfGen::new(prec);
        let nums = vec![
            Sqrt5Mpz::from_si_g(630),
            Sqrt5Mpz::from_si_g(-840),
            Sqrt5Mpz::from_si_g(-1),
        ];
        let forms = vec![f, h0, h1];
        for (f, n) in forms.iter().zip(nums.iter()) {
            tmp.mul_mut_by_const(f, n);
            res += &tmp;
        }
        println!("res: {}", res);
        assert!(res.is_zero());
    }

    #[test]
    fn test_pickle_gen6() {
        let gens1 = Structure6::gens1(10);
        let rel1 = relation_slow_3gens(&gens1, 50);
        let ref mut f1 = File::create("./data/str6gens1.sobj").unwrap();
        save_as_pickle_rel3(&rel1, f1);
    }

    #[test]
    fn test_gens_relation6() {
        let prec = 12;
        let gens1 = Structure6::gens1(prec);
        let f2 = &gens1[0];
        let f7 = &gens1[1];
        let f6 = &gens1[2];
        let f11 = rankin_cohen_sqrt5(6, &g5_normalized(prec), &f6_normalized(prec)).unwrap();
        let g2: HmfGen<Sqrt5Mpz> = From::from(&eisenstein_series(2, prec));
        let g5: HmfGen<Sqrt5Mpz> = From::from(&g5_normalized(prec));
        let mut tmp = HmfGen::new(prec);
        tmp.pow_mut(&g2, 2);
        tmp *= &g5;
        let h0 = f2 * &tmp;
        tmp.pow_mut(&g2, 2);
        let h1 = f7 * &tmp;
        let h2 = f6 * &g5;
        let mut forms = vec![h0, h1, h2];
        for f in forms.iter() {
            assert_eq!(f.weight, Some((11, 23)));
        }
        forms.push(f11);
        println!("{:?}", forms.len() as i64 - rank(50, &forms) as i64);
    }

    #[test]
    fn test_gens_relation6_1() {
        let prec = 12;
        let gens1 = Structure6::gens1(prec);
        let f2 = &gens1[0];
        let f7 = &gens1[1];
        let f6 = &gens1[2];
        let ref f11 = rankin_cohen_sqrt5(6, &g5_normalized(prec), &f6_normalized(prec)).unwrap();
        let g2: HmfGen<Sqrt5Mpz> = From::from(&eisenstein_series(2, prec));
        let g5: HmfGen<Sqrt5Mpz> = From::from(&g5_normalized(prec));
        let g6: HmfGen<Sqrt5Mpz> = From::from(&f6_normalized(prec));
        let mut tmp = HmfGen::new(prec);
        tmp.pow_mut(&g2, 3);
        let h0 = f2 * &(&g5 * &tmp);
        let h1 = f2 * &(&g5 * &g6);
        let h2 = f7 * &tmp;
        let h3 = f7 * &g6;
        let h4 = f6 * &(&g5 * &g2);
        let forms = vec![h0, h1, h2, h3, h4];
        for f in forms.iter() {
            assert_eq!(f.weight, Some((13, 25)));
        }
        let f = f11 * &g2;
        let rel = relation(50, &f, &forms);
        for a in rel.iter() {
            assert!(a.ir_part().is_zero());
        }
        let rel: Vec<_> = rel.iter().map(|a| a.rt_part() >> 1).collect();
        println!("{:?}", rel);
    }

    #[test]
    fn test_brackets6() {
        let prec = 10;
        let gens1 = Structure6::gens1(prec);
        let f2 = &gens1[0];
        let f7 = &gens1[1];
        let f6 = &gens1[2];
        let ref f11 = rankin_cohen_sqrt5(6, &g5_normalized(prec), &f6_normalized(prec)).unwrap();
        let h10 = bracket_inner_prod1(f2, f11).unwrap();
        let h15 = bracket_inner_prod1(f7, f11).unwrap();
        let h14 = bracket_inner_prod1(f6, f11).unwrap();
        let (_, v10) = relation_monom(50, &h10);
        let (_, v15) = relation_monom(50, &h15);
        let (_, v14) = relation_monom(50, &h14);
        let ref mut f = File::create("./data/str6bra.sobj").unwrap();
        save_as_pickle_rel3(&(v10, v15, v14), f);
        // println!("{:?}", v10);
        // println!("{:?}", v15);
        // println!("{:?}", v14);
    }

    #[test]
    fn test_ranks6() {
        for k in 2..100 {
            let forms = forms6(k, 20);
            println!("{}: {}", k, rank(200, &forms));
        }
    }

    #[test]
    fn test_gens6_1() {
        let prec = 15;
        let mut tmp = HmfGen::new(prec);
        let mut tmp1 = HmfGen::new(prec);
        let g2 = eisenstein_series(2, prec);
        let g6 = f6_normalized(prec);
        tmp.pow_mut(&g2, 10);
        tmp1.pow_mut(&g6, 10);
        let f = rankin_cohen_sqrt5(6, &tmp, &tmp1).unwrap();
        let forms = forms6(f.weight.unwrap().0, prec);
        let rel = relation(100, &f, &forms);
        println!("{:?}", rel);
    }

    fn forms6(k: usize, prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut res = Vec::new();
        let gens1 = Structure6::gens1(prec);
        let f2 = &gens1[0];
        let f7 = &gens1[1];
        let f6 = &gens1[2];
        let ref f11 = rankin_cohen_sqrt5(6, &g5_normalized(prec), &f6_normalized(prec)).unwrap();
        fn append(res: &mut Vec<HmfGen<Sqrt5Mpz>>, f: &HmfGen<Sqrt5Mpz>, k: usize, prec: usize) {
            let l = f.weight.unwrap().0;
            if k == l {
                let mut tmp = HmfGen::new(prec);
                tmp.set(f);
                res.push(tmp);
            } else if k > l {
                for a in monoms_of_g2_g5_f6(k - f.weight.unwrap().0).iter().map(
                    |x| {
                        x.into_form(prec)
                    },
                )
                {
                    let mut tmp: HmfGen<Sqrt5Mpz> = From::from(&a);
                    tmp *= f;
                    res.push(tmp);
                }
            }
        }
        append(&mut res, f2, k, prec);
        append(&mut res, f7, k, prec);
        append(&mut res, f6, k, prec);
        append(&mut res, f11, k, prec);
        res
    }

    #[test]
    fn test_save_rels() {
        let prec = 10;
        let ref mut f = File::create("./data/rels.sobj").unwrap();
        let v: Vec<_> = (1..10)
            .map(|i| (i, three_forms_rel(i, prec, 50)))
            .take_while(|x| x.1.is_some())
            .map(|x| (x.0, x.1.unwrap()))
            .collect();
        save_as_pickle_3relations(&v, f);
    }
}

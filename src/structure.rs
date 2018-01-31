use gmp::mpz::Mpz;
use eisenstein::{eisenstein_series, f6_normalized};
use theta_chars::g5_normalized;
use elements::{HmfGen, div_mut_with_denom, initial_term};
use serde_pickle;
use std::fs::File;
use std::io::Write;
use std::io::Read;
use bignum::{RealQuadElement, BigNumber};
use serde;
use std::process::Command;
use std::env;
use bignum::Sqrt5Mpz;
use diff_op::{rankin_cohen_sqrt5, bracket_inner_prod1, star_op};
use rand;
use std::path::Path;
use std::fs::remove_file;
use serde::ser::{Serialize, Serializer};
use serde::{Deserialize, Deserializer};
use libc::c_long;
use flint::fmpz::Fmpz;
use flint::fmpz_mat::FmpzMat;
use std::convert::From;
use std::ops::AddAssign;
use std::fmt;
use std::collections::HashMap;
use flint::traits::*;

pub struct MpzWrapper {
    pub a: Mpz,
}

impl Serialize for MpzWrapper {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        String::serialize(&self.a.to_str_radix(10), serializer)
    }
}

impl<'de> Deserialize<'de> for MpzWrapper {
    fn deserialize<D>(deserializer: D) -> Result<MpzWrapper, D::Error>
    where
        D: Deserializer<'de>,
    {
        let a: String = String::deserialize(deserializer)?;
        let res = Mpz::from_str_radix(&a, 10).unwrap();
        Ok(MpzWrapper { a: res })
    }
}

impl<'b> From<&'b Mpz> for MpzWrapper {
    fn from(a: &Mpz) -> MpzWrapper {
        MpzWrapper { a: a.clone() }
    }
}

impl<'b> From<&'b MpzWrapper> for Mpz {
    fn from(a: &MpzWrapper) -> Mpz {
        a.a.clone()
    }
}

pub struct Sqrt5Wrapper {
    pub a: Sqrt5Mpz,
}

impl<'b> From<&'b Sqrt5Wrapper> for Sqrt5Mpz {
    fn from(a: &Sqrt5Wrapper) -> Sqrt5Mpz {
        a.a.clone()
    }
}

impl<'b> From<&'b Sqrt5Mpz> for Sqrt5Wrapper {
    fn from(a: &Sqrt5Mpz) -> Sqrt5Wrapper {
        Sqrt5Wrapper { a: a.clone() }
    }
}

impl Serialize for Sqrt5Wrapper {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let tpl = (
            self.a.rt_part().to_str_radix(10),
            self.a.ir_part().to_str_radix(10),
        );
        tpl.serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for Sqrt5Wrapper {
    fn deserialize<D>(deserializer: D) -> Result<Sqrt5Wrapper, D::Error>
    where
        D: Deserializer<'de>,
    {
        type StringTuple = (String, String);
        let s: (String, String) = StringTuple::deserialize(deserializer)?;
        let a1 = Mpz::from_str_radix(&s.0, 10).unwrap();
        let a2 = Mpz::from_str_radix(&s.1, 10).unwrap();
        Ok(Sqrt5Wrapper { a: From::from(&(a1, a2)) })
    }
}

/// A stupid function that returns linear relations.
pub fn relations(len: usize, forms: &[HmfGen<Sqrt5Mpz>]) -> Vec<Vec<Sqrt5Mpz>> {
    let vv: Vec<Vec<Sqrt5Mpz>> = forms.iter().map(|f| f.fc_vector(len)).collect();
    let path_name = format!("/tmp/rust_python_data{}.sobj", rand::random::<u64>()).to_string();
    let path = Path::new(&path_name);
    assert!(!path.exists());
    {
        let mut f = File::create(&path_name).unwrap();
        let vv: Vec<Vec<Sqrt5Wrapper>> = vv.iter()
            .map(|v| v.iter().map(From::from).collect())
            .collect();
        save_as_pickle(vv, &mut f);
        let _rank = sage_command(&format!(
            "data_name = '{}'; print load('./src/relation.sage')",
            path_name
        )).lines()
            .next()
            .unwrap()
            .to_string();
    }
    let res = {
        let f = File::open(&path_name).unwrap();
        let vv: Vec<Vec<Sqrt5Wrapper>> = load_pickle(&f).unwrap();
        vv.iter()
            .map(|v| v.iter().map(From::from).collect())
            .collect()
    };
    remove_file(&path).unwrap();
    res
}

pub fn relations_over_z(forms: &[HmfGen<Mpz>]) -> Vec<Vec<Mpz>> {
    let vv: Vec<_> = forms.iter().map(|f| f.fc_vector_u_nonneg()).collect();
    let vv: Vec<Vec<Fmpz>> = vv.iter()
        .map(|v| v.iter().map(|x| From::from(x)).collect())
        .collect();
    let n = forms.len();
    let m = vv[0].len();
    let mut mat = FmpzMat::new(m as i64, n as i64);
    for (i, v) in vv.iter().enumerate() {
        for (j, x) in v.iter().enumerate() {
            mat.set_entry(j as isize, i as isize, x);
        }
    }
    mat.nullspace_basis()
        .iter()
        .map(|v| v.iter().map(Into::<Mpz>::into).collect())
        .collect()
}

// TODO: Remove duplicate of code.
/// Similar to `r_elt_as_pol`.
pub fn r_elt_as_pol_over_z(f: &HmfGen<Mpz>) -> Option<(Vec<(MonomFormal, Mpz)>, Mpz)> {
    let f = f.clone();
    let prec = f.prec;
    let monoms = monoms_of_g2_g5_f6(f.weight.unwrap().0);
    let mut forms: Vec<_> = monoms.iter().map(|x| x.into_form(prec)).collect();
    forms.insert(0, f);
    let mut rels = relations_over_z(&forms);
    if rels.len() == 1 && !rels[0][0].is_zero() {
        let cfs: Vec<_> = rels[0].iter().skip(1).map(|x| -x).collect();
        Some((
            monoms.into_iter().zip(cfs.into_iter()).collect(),
            rels.remove(0).remove(0),
        ))
    } else {
        None
    }
}

pub fn r_elt_as_pol_over_z_cached_gens(
    f: &HmfGen<Mpz>,
    map_g2: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    map_g5: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    map_g6: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
) -> Option<(Vec<(MonomFormal, Mpz)>, Mpz)> {
    let f = f.clone();
    let prec = f.prec;
    let monoms = monoms_of_g2_g5_f6(f.weight.unwrap().0);
    let mut forms: Vec<_> = monoms
        .iter()
        .map(|x| x.into_form_cached(prec, map_g2, map_g5, map_g6))
        .collect();
    forms.insert(0, f);
    let mut rels = relations_over_z(&forms);
    if rels.len() == 1 && !rels[0][0].is_zero() {
        let cfs: Vec<_> = rels[0].iter().skip(1).map(|x| -x).collect();
        Some((
            monoms.into_iter().zip(cfs.into_iter()).collect(),
            rels.remove(0).remove(0),
        ))
    } else {
        None
    }
}

/// f: polynomial of `g2, g5, g6` as q-expansion. Return the corresponding
/// polynomia. The second element is a denominator.
pub fn r_elt_as_pol(f: &HmfGen<Sqrt5Mpz>, len: usize) -> Option<(PWtPoly, Sqrt5Mpz)> {
    let f = f.clone();
    let prec = f.prec;
    let monoms = monoms_of_g2_g5_f6(f.weight.unwrap().0);
    let mut forms: Vec<_> = monoms
        .iter()
        .map(|x| From::from(&x.into_form(prec)))
        .collect();
    forms.insert(0, f);
    let mut rels = relations(len, &forms);
    if rels.len() == 1 && !rels[0][0].is_zero_g() {
        let cfs: Vec<_> = rels[0]
            .iter()
            .skip(1)
            .map(|x| {
                let mut tmp = Sqrt5Mpz::new_g();
                tmp.set_g(x);
                tmp *= -1 as c_long;
                tmp
            })
            .collect();
        Some((
            monoms.into_iter().zip(cfs.into_iter()).collect(),
            rels.remove(0).remove(0),
        ))
    } else {
        None
    }
}

pub fn bracket_inner_prod_as_pol(
    f: &HmfGen<Sqrt5Mpz>,
    g: &HmfGen<Sqrt5Mpz>,
    len: usize,
) -> Option<(PWtPoly, Sqrt5Mpz)> {
    let h = bracket_inner_prod1(f, g).unwrap();
    r_elt_as_pol(&h, len)
}

pub fn bracket_inner_prod_as_pol_over_z_maybe(
    f: &HmfGen<Sqrt5Mpz>,
    g: &HmfGen<Sqrt5Mpz>,
) -> Option<(PWtPolyZ, Mpz)> {
    let h = bracket_inner_prod1(f, g).unwrap();
    if h.is_zero() {
        return Some((vec![], From::from(1)));
    }
    if !h.rt_part().is_zero() {
        None
    } else {
        let h_ir = h.ir_part();
        r_elt_as_pol_over_z(&h_ir)
    }
}

/// A stupid function that returns a linear relation.
pub fn relation(len: usize, f: &HmfGen<Sqrt5Mpz>, forms: &[HmfGen<Sqrt5Mpz>]) -> Vec<Sqrt5Mpz> {
    let mut vv = Vec::new();
    vv.push(f.clone());
    for g in forms.iter() {
        vv.push(g.clone());
    }
    let mut res = relations(len, &vv);
    res.swap_remove(0)
}

pub type PWtPoly = Vec<(MonomFormal, Sqrt5Mpz)>;
pub type PWtPolyZ = Vec<(MonomFormal, Mpz)>;
pub type Relation = Vec<PWtPoly>;

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

fn sage_command(cmd: &String) -> String {
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
pub fn tpls_of_wt(k: usize) -> Vec<(usize, usize, usize)> {
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

fn pow_cached(
    prec: usize,
    expt: usize,
    f: &HmfGen<Mpz>,
    cached_map: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
) -> HmfGen<Mpz> {
    if let Some(v) = cached_map.get(&(prec, expt)) {
        return v.clone();
    }
    let mut tmp = HmfGen::new(prec);
    if let Some(k) = cached_map
        .keys()
        .map(|k| k.clone())
        .filter(|&k| k.0 == prec && k.1 < expt)
        .max()
    {
        tmp.pow_mut(f, expt - k.1);
        let v = {
            cached_map.get(&k).unwrap().clone()
        };
        tmp *= &v;
        cached_map.insert((prec, expt), tmp.clone());
        return tmp;
    }
    tmp.pow_mut(f, expt);
    cached_map.insert((prec, expt), tmp.clone());
    tmp
}

fn monom_g2_g5_g6_cached(
    prec: usize,
    expts: (usize, usize, usize),
    map_g2: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    map_g5: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    map_g6: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
) -> HmfGen<Mpz> {
    let g2 = eisenstein_series(2, prec);
    let g5 = g5_normalized(prec);
    let g6 = f6_normalized(prec);
    let mut g2_pow = pow_cached(prec, expts.0, &g2, map_g2);
    let g5_pow = pow_cached(prec, expts.1, &g5, map_g5);
    let g6_pow = pow_cached(prec, expts.2, &g6, map_g6);
    g2_pow *= &g5_pow;
    g2_pow *= &g6_pow;
    g2_pow
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
fn rel3_to_tuple(rel: &[PWtPoly; 3]) -> (PolString, PolString, PolString) {
    let p0 = &rel[0];
    let p1 = &rel[1];
    let p2 = &rel[2];
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

pub fn save_as_pickle_3relations(rels: &Vec<(usize, [PWtPoly; 3])>, f: &mut File) {
    let v: Vec<_> = rels.iter()
        .map(|&(i, ref rel)| (i, rel3_to_tuple(&rel)))
        .collect();
    save_as_pickle(v, f);
}

pub fn save_as_pickle_rel(rel: &[PWtPoly], f: &mut File) {
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
    let v: Vec<Vec<_>> = rel.iter().map(to_vec).collect();
    save_as_pickle(v, f);
}

pub fn save_as_pickle<T>(a: T, f: &mut File)
where
    T: serde::Serialize,
{
    let v = serde_pickle::to_vec(&a, false).unwrap();
    f.write(&v).unwrap();
}

fn load_pickle<'de, T>(f: &File) -> Result<Vec<T>, serde_pickle::Error>
where
    T: serde::Deserialize<'de>,
{
    let buf: Vec<u8> = f.bytes().map(|x| x.unwrap()).collect();
    serde_pickle::from_slice(&buf)
}

#[derive(Debug)]
pub struct StrCand {
    pub df: u64,
    pub gens_num: Vec<(PWtPolyZ, PWtPolyZ)>,
    pub rels: Vec<Vec<Vec<PWtPolyZ>>>,
    pub gen_wts: Vec<u64>,
    pub rel_wts: Vec<Vec<u64>>,
    pub free_basis_wts: ((usize, u64), (usize, u64)),
    pub gens_dnm: PWtPolyZ,
    pub monoms: ((MonomFormal, MonomFormal), (MonomFormal, MonomFormal)),
}

pub type PWtPolyZRaw = Vec<((usize, usize, usize), String)>;

impl StrCand {
    pub fn load(df: u64, f: &File, monom_fl: &File) -> Result<Self, serde_pickle::Error> {
        let data: (Vec<Vec<Vec<PWtPolyZRaw>>>,
                   Vec<Vec<u64>>,
                   ((usize, u64), (usize, u64), PWtPolyZRaw)) = {
            let buf: Vec<u8> = f.bytes().map(|x| x.unwrap()).collect();
            serde_pickle::from_slice(&buf)?
        };

        let monom_data: Vec<((usize, usize, usize), (usize, usize, usize))> = {
            let buf: Vec<u8> = monom_fl.bytes().map(|x| x.unwrap()).collect();
            serde_pickle::from_slice(&buf)?
        };

        let monoms_vec: Vec<_> = monom_data
            .into_iter()
            .map(|x| (MonomFormal { idx: x.0 }, MonomFormal { idx: x.1 }))
            .collect();
        let free_basis_wts = ((data.2).0, (data.2).1);
        let monoms = (
            monoms_vec[(free_basis_wts.0).0].clone(),
            monoms_vec[(free_basis_wts.1).0].clone(),
        );


        fn to_monom_formal_tpl(
            xs: &Vec<((usize, usize, usize), String)>,
        ) -> Vec<(MonomFormal, Mpz)> {
            xs.iter()
                .map(|x| {
                    (
                        MonomFormal { idx: x.0 },
                        Mpz::from_str_radix(&x.1, 10).unwrap(),
                    )
                })
                .collect()
        }

        let gens_nums: Vec<_> = data.0[0]
            .iter()
            .map(|v| {
                let v: Vec<_> = v.iter().map(to_monom_formal_tpl).collect();
                (v[0].clone(), v[1].clone())
            })
            .collect();
        let rels: Vec<Vec<Vec<PWtPolyZ>>> = data.0
            .clone()
            .into_iter()
            .skip(1)
            .map(|l| {
                l.iter()
                    .map(|v| v.iter().map(to_monom_formal_tpl).collect())
                    .collect()
            })
            .collect();
        let mut wts = data.1.clone();
        let gen_wts = wts.remove(0);
        let rel_wts = wts;
        let gens_dnm = to_monom_formal_tpl(&(data.2).2);
        Ok(StrCand {
            df: df,
            gens_num: gens_nums,
            rels: rels,
            gen_wts: gen_wts,
            rel_wts: rel_wts,
            free_basis_wts: free_basis_wts,
            gens_dnm: gens_dnm,
            monoms: monoms,
        })
    }

    pub fn gens_nums_as_forms(&self, prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        fn to_pwtpoly(x: &PWtPolyZ) -> PWtPoly {
            x.iter().map(|y| (y.0.clone(), From::from(&y.1))).collect()
        }
        let free_basis = self.free_basis_nums(prec);
        self.gens_num
            .iter()
            .map(|nums| {
                let v = [to_pwtpoly(&nums.0), to_pwtpoly(&nums.1)];
                linear_comb(&v, &free_basis)
            })
            .collect()
    }

    pub fn free_basis_nums(&self, prec: usize) -> [HmfGen<Sqrt5Mpz>; 2] {
        let (_, wt_f) = self.free_basis_wts.0;
        let (_, wt_g) = self.free_basis_wts.1;
        let ((ref m0, ref n0), (ref m1, ref n1)) = self.monoms;
        let f = rankin_cohen_sqrt5(self.df as usize, &m0.into_form(prec), &n0.into_form(prec))
            .unwrap();
        let g = rankin_cohen_sqrt5(self.df as usize, &m1.into_form(prec), &n1.into_form(prec))
            .unwrap();
        assert_eq!(f.weight.unwrap().0 as u64, wt_f);
        assert_eq!(g.weight.unwrap().0 as u64, wt_g);
        [f, g]
    }

    pub fn gens_with_const(&self, prec: usize) -> Vec<(HmfGen<Sqrt5Mpz>, Sqrt5Mpz, Sqrt5Mpz)> {
        let v = self.gens_nums_as_forms(prec);
        let mut prec_small = 5;
        let mut prec = prec;
        loop {
            let dnm_form = MonomFormal::eval(&self.gens_dnm, prec_small);
            let a = initial_term(&dnm_form);
            if let Some((v, _, _)) = a {
                prec += v;
                break;
            } else {
                prec_small += 1;
            }
        }
        let mut tmp = Mpz::new();
        let dnm_form = From::from(&MonomFormal::eval(&self.gens_dnm, prec));
        v.iter()
            .map(|f| {
                let mut res = HmfGen::new(prec);
                let mut dnm = div_mut_with_denom(&mut res, f, &dnm_form, true);
                let mut num = From::from(&res.gcd());
                res /= &num;
                let a11 = res.fourier_coefficient(1, 1);
                if !a11.is_zero_g() && res.is_divisible_by_const(&a11) {
                    res /= &a11;
                    num.mul_assign_g(&a11, &mut tmp);
                }
                let a: Sqrt5Mpz = From::from(&(num.content().gcd(&dnm.content())));
                num.set_divexact_g(&a, &mut tmp);
                dnm.set_divexact_g(&a, &mut tmp);
                (res, num, dnm)
            })
            .collect()
    }

    pub fn gens_normalized(&self, prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        self.gens_with_const(prec)
            .iter()
            .map(|x| x.0.clone())
            .collect()
    }

    pub fn gens_nums_wts(&self) -> Vec<(u64, u64)> {
        let (_, wt_f) = self.free_basis_wts.0;
        let (_, wt_g) = self.free_basis_wts.1;
        fn wt(m: &(PWtPolyZ, PWtPolyZ), wt_f: u64, wt_g: u64) -> u64 {
            if m.0.is_empty() {
                m.1.last().unwrap().0.weight() as u64 + wt_g
            } else {
                m.0.last().unwrap().0.weight() as u64 + wt_f
            }
        }
        self.gens_num
            .iter()
            .map(|x| {
                (
                    wt(x, wt_f as u64, wt_g as u64),
                    wt(x, wt_f as u64, wt_g as u64) + 2 * self.df as u64,
                )
            })
            .collect()
    }

    pub fn save_star_norms(&self, gens: &[HmfGen<Sqrt5Mpz>], path: &str) {
        let mut pols = Vec::new();
        for f in gens.iter() {
            let mut nm = HmfGen::new(f.prec);
            star_op(&mut nm, f);
            nm *= f;
            let wt = nm.weight.unwrap();
            assert_eq!(wt.0, wt.1);
            assert!(nm.ir_part().is_zero());
            let poly = r_elt_as_pol_over_z(&nm.rt_part()).unwrap();
            pols.push(poly);
        }
        let mut stars_f = File::create(path).unwrap();
        save_polys_over_z_pickle(&pols, &mut stars_f);
    }

    fn gens_normalized1(&self, prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut gens_wconst = self.gens_with_const(prec);
        let mut tmp = Mpz::new();
        let m = gens_wconst.iter().map(|x| x.2.clone()).fold(
            From::from((2, 0)),
            |acc, x| {
                let mut a = Sqrt5Mpz::new_g();
                a.set_g(&acc);
                a.mul_assign_g(&x, &mut tmp);
                a
            },
        );
        for x in gens_wconst.iter_mut() {
            let mut m1 = m.clone();
            m1.set_divexact_g(&x.2, &mut tmp);
            m1.mul_assign_g(&x.1, &mut tmp);
            x.0 *= &m1;
        }
        gens_wconst.iter().map(|x| x.0.clone()).collect()
    }

    pub fn test_relations(&self, prec: usize) {
        assert!(self.rels.len() <= 1);
        fn to_pwtpoly(x: &PWtPolyZ) -> PWtPoly {
            x.iter().map(|y| (y.0.clone(), From::from(&y.1))).collect()
        }
        let gens = self.gens_normalized1(prec);
        if self.rels.len() > 0 {
            for rel in &self.rels[0] {
                let rel: Vec<_> = rel.iter().map(to_pwtpoly).collect();
                let rel_form = linear_comb(&rel, &gens);
                assert!(rel_form.is_zero());
            }
        }
    }
}

pub fn load_cand(i: u64) -> StrCand {
    let cand_f = File::open(format!("./data/brackets/str{}_cand.sobj", i)).unwrap();
    let monom_f = File::open(format!("./data/brackets/str{}_monoms.sobj", i)).unwrap();
    StrCand::load(i, &cand_f, &monom_f).unwrap()
}

pub fn save_polys_over_z_pickle(xs: &[(PWtPolyZ, Mpz)], f: &mut File) {
    let v: Vec<(Vec<_>, MpzWrapper)> = xs.iter()
        .map(|x| {
            (
                x.0
                    .iter()
                    .map(|elt| (elt.0.idx, Into::<MpzWrapper>::into(&elt.1)))
                    .collect(),
                Into::into(&x.1),
            )
        })
        .collect();
    save_as_pickle(&v, f);
}

/// Corresponds to g2^a * g5^b * f6^c where (a, b, c) = idx.
#[derive(Clone, Debug)]
pub struct MonomFormal {
    pub idx: (usize, usize, usize),
}

impl MonomFormal {
    pub fn into_form(&self, prec: usize) -> HmfGen<Mpz> {
        monom_g2_g5_f6(prec, self.idx)
    }

    pub fn into_form_cached(
        &self,
        prec: usize,
        map_g2: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
        map_g5: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
        map_g6: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    ) -> HmfGen<Mpz> {
        monom_g2_g5_g6_cached(prec, self.idx, map_g2, map_g5, map_g6)
    }

    fn weight(&self) -> usize {
        (2 * self.idx.0 + 5 * self.idx.1 + 6 * self.idx.2)
    }

    pub fn eval<T>(v: &Vec<(MonomFormal, T)>, prec: usize) -> HmfGen<T>
    where
        T: BigNumber + Clone + fmt::Debug,
        for<'b> T: From<&'b Mpz>,
        for<'b> T: AddAssign<&'b T>,
    {
        if v.is_empty() {
            HmfGen::new(prec)
        } else {
            let w = v[0].0.weight();
            let wt = if v.iter().all(|x| x.0.weight() == w) {
                Some((w, w))
            } else {
                None
            };
            let mut res = HmfGen::new(prec);
            res.set(&v[0].0.into_form(prec));
            let mut res: HmfGen<T> = From::from(&res);
            res *= &v[0].1;
            for &(ref monom, ref a) in v.iter().skip(1) {
                let mut tmp: HmfGen<T> = From::from(&monom.into_form(prec));
                tmp *= a;
                res += &tmp;
            }
            res.weight = wt;
            res
        }
    }

    pub fn eval_cached<T>(
        v: &Vec<(MonomFormal, T)>,
        prec: usize,
        map_g2: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
        map_g5: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
        map_g6: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    ) -> HmfGen<T>
    where
        T: BigNumber + Clone + fmt::Debug,
        for<'b> T: From<&'b Mpz>,
        for<'b> T: AddAssign<&'b T>,
    {
        if v.is_empty() {
            HmfGen::new(prec)
        } else {
            let w = v[0].0.weight();
            let wt = if v.iter().all(|x| x.0.weight() == w) {
                Some((w, w))
            } else {
                None
            };
            let mut res = HmfGen::new(prec);
            res.set(&v[0].0.into_form_cached(prec, map_g2, map_g5, map_g6));
            let mut res: HmfGen<T> = From::from(&res);
            res *= &v[0].1;
            for &(ref monom, ref a) in v.iter().skip(1) {
                let mut tmp: HmfGen<T> =
                    From::from(&monom.into_form_cached(prec, map_g2, map_g5, map_g6));
                tmp *= a;
                res += &tmp;
            }
            res.weight = wt;
            res
        }
    }
}


pub fn monoms_of_g2_g5_f6(k: usize) -> Vec<MonomFormal> {
    tpls_of_wt(k)
        .iter()
        .map(|&x| MonomFormal { idx: x })
        .collect()
}

fn linear_comb(coeffs: &[PWtPoly], gens: &[HmfGen<Sqrt5Mpz>]) -> HmfGen<Sqrt5Mpz> {
    let prec = gens[0].prec;
    let mut res = HmfGen::<Sqrt5Mpz>::new(prec);
    for (&ref p, f) in coeffs.iter().zip(gens.iter()) {
        if !p.is_empty() {
            let mut tmp = MonomFormal::eval(p, prec);
            tmp *= f;
            res += &tmp;
        }
    }
    assert!(res.weight.is_some());
    res
}

#[allow(dead_code)]
fn linear_comb_cached(
    coeffs: &[PWtPoly],
    gens: &[HmfGen<Sqrt5Mpz>],
    map_g2: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    map_g5: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
    map_g6: &mut HashMap<(usize, usize), HmfGen<Mpz>>,
) -> HmfGen<Sqrt5Mpz> {
    let prec = gens[0].prec;
    let mut res = HmfGen::<Sqrt5Mpz>::new(prec);
    for (&ref p, f) in coeffs.iter().zip(gens.iter()) {
        if !p.is_empty() {
            let mut tmp = MonomFormal::eval_cached(p, prec, map_g2, map_g5, map_g6);
            tmp *= f;
            res += &tmp;
        }
    }
    assert!(res.weight.is_some());
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
                if !linear_comb(rel, &gens).is_zero() {
                    return false;
                }
            }
            true
        }
    }
}

/// Return a vector of length `len` of mixed weight modular forms.
pub fn mixed_weight_forms(
    df: usize,
    prec: usize,
    len: usize,
) -> Vec<(HmfGen<Sqrt5Mpz>, (usize, usize, usize), (usize, usize, usize))> {
    let mut num = 0;
    let mut res = Vec::new();
    for (i, m) in (2..).flat_map(monoms_of_g2_g5_f6).enumerate() {
        for n in (2..).flat_map(monoms_of_g2_g5_f6).take(if is_even!(df) {
            i + 1
        } else {
            i
        })
        {
            if num >= len {
                return res;
            }
            let f = rankin_cohen_sqrt5(df, &m.into_form(prec), &n.into_form(prec)).unwrap();
            if !f.is_zero() {
                num += 1;
                res.push((f, m.idx, n.idx));
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let v: Vec<_> = (2..).flat_map(monoms_of_g2_g5_f6).take(10).collect();
        println!("{:?}", v);
    }
}

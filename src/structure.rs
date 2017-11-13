use gmp::mpz::Mpz;
use eisenstein::{eisenstein_series, f6_normalized};
use theta_chars::g5_normalized;
use elements::{HmfGen, div_mut, div_mut_with_denom, initial_term};
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
use misc::PowGen;
use rand;
use std::path::Path;
use std::fs::remove_file;
use serde::ser::{Serialize, Serializer};
use serde::{Deserialize, Deserializer};
use libc::{c_ulong, c_long};
use flint::fmpz::Fmpz;
use flint::fmpz_mat::FmpzMat;
use std::convert::From;
use std::ops::AddAssign;
use std::fmt;
use std::collections::HashMap;

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
        let (_, wt_f) = self.free_basis_wts.0;
        let (_, wt_g) = self.free_basis_wts.1;
        let ((ref m0, ref n0), (ref m1, ref n1)) = self.monoms;
        let f = rankin_cohen_sqrt5(self.df as usize, &m0.into_form(prec), &n0.into_form(prec))
            .unwrap();
        let g = rankin_cohen_sqrt5(self.df as usize, &m1.into_form(prec), &n1.into_form(prec))
            .unwrap();
        assert_eq!(f.weight.unwrap().0 as u64, wt_f);
        assert_eq!(g.weight.unwrap().0 as u64, wt_g);
        let free_basis = [f, g];
        self.gens_num
            .iter()
            .map(|nums| {
                let v = [to_pwtpoly(&nums.0), to_pwtpoly(&nums.1)];
                linear_comb(&v, &free_basis)
            })
            .collect()
    }

    pub fn gens(&self, prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
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
        let dnm_form = From::from(&MonomFormal::eval(&self.gens_dnm, prec));
        v.iter()
            .map(|f| {
                let mut res = HmfGen::new(prec);
                div_mut_with_denom(&mut res, f, &dnm_form, true);
                let a = From::from(&res.gcd());
                res /= &a;
                res
            })
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

pub fn three_forms(i: usize, prec: usize) -> Option<Vec<HmfGen<Sqrt5Mpz>>> {
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

pub fn three_forms_rel(i: usize, prec: usize, len: usize) -> Option<[PWtPoly; 3]> {
    let gens = three_forms(i, prec);
    gens.map(|ref x| relation_slow_3gens(x, len))
}

pub fn relation_slow_3gens(gens: &Vec<HmfGen<Sqrt5Mpz>>, len: usize) -> [PWtPoly; 3] {
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
    [v0, v1, v2]
}

pub fn print_3rel(rel: &[PWtPoly]) {
    println!("{:?}", rel[0]);
    println!("{:?}", rel[1]);
    println!("{:?}", rel[2]);
}

impl Structure4 {
    #[allow(dead_code)]
    fn relation_slow() {
        let prec = 10;
        let gens = Self::gens(prec);
        print_3rel(&relation_slow_3gens(&gens, 50));
    }
}

pub struct Structure3;

impl Structure3 {
    // This is not gens.
    pub fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
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
        vec![f2, f0]
    }

    fn relations() -> Option<Vec<Relation>> {
        None
    }
}

pub struct Structure5;

impl Structure for Structure5 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        Self::gens2(prec)
    }

    fn relations() -> Option<Vec<Relation>> {
        let a0 = (MonomFormal { idx: (1, 1, 0) }, Sqrt5Mpz::from_si_g(1));
        let a1 = (
            MonomFormal { idx: (0, 0, 1) },
            Sqrt5Mpz::from_si_g(29937600),
        );
        let a2 = (
            MonomFormal { idx: (0, 1, 0) },
            Sqrt5Mpz::from_si_g(18144000),
        );
        let a3 = (MonomFormal { idx: (1, 0, 0) }, Sqrt5Mpz::from_si_g(-165));

        let b0_1 = (MonomFormal { idx: (3, 0, 0) }, Sqrt5Mpz::from_si_g(-1));
        let b0_2 = (MonomFormal { idx: (0, 0, 1) }, Sqrt5Mpz::from_si_g(1080));
        let b1 = (
            MonomFormal { idx: (0, 1, 0) },
            Sqrt5Mpz::from_si_g(1632960000),
        );
        let b2 = (MonomFormal { idx: (2, 0, 0) }, Sqrt5Mpz::from_si_g(1814400));
        let b3 = (MonomFormal { idx: (0, 0, 0) }, Sqrt5Mpz::from_si_g(0));
        Some(vec![
            vec![vec![a0], vec![a1], vec![a2], vec![a3]],
            vec![vec![b0_1, b0_2], vec![b1], vec![b2], vec![b3]],
        ])

    }
}

impl Structure5 {
    pub fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let g6 = f6_normalized(prec);
        let f7 = rankin_cohen_sqrt5(5, &g2, &g5).unwrap();
        let f7_1 = &Structure2::gens(prec)[0] * &Structure3::gens(prec)[0];
        let f8 = rankin_cohen_sqrt5(5, &g2, &g6).unwrap();
        let res = vec![f7, f7_1, f8];
        for f in res.iter() {
            assert!({
                let (a1, a2) = f.weight.unwrap();
                a2 - a1 == 10
            })
        }
        res
    }

    pub fn gens2(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut gens = Self::gens1(prec);
        let g2 = eisenstein_series(2, prec);
        let mut f6 = HmfGen::new(prec);
        let mut f5 = HmfGen::new(prec);

        {
            let f8 = &gens[2];
            let f7 = &gens[0];
            let f7_1 = &gens[1];

            div_mut(&mut f6, f8, &From::from(&g2));

            let mut f7_2 = f7_1.clone();
            f7_2 *= &Sqrt5Mpz::from_si_g(11880);
            f7_2 += &(f7 * &Sqrt5Mpz::from_si_g(1814400));
            div_mut(&mut f5, &f7_2, &From::from(&g2));
        }
        let f7 = gens.swap_remove(0);
        let gens2 = Structure2::gens(prec);
        let gens3 = Structure3::gens(prec);
        let f10 = &gens2[0] * &gens3[1];
        assert_eq!(f5.weight, Some((5, 15)));
        assert_eq!(f6.weight, Some((6, 16)));
        assert_eq!(f7.weight, Some((7, 17)));
        assert_eq!(f10.weight, Some((10, 20)));
        vec![f5, f6, f7, f10]
    }
}

pub struct Structure6;

impl Structure for Structure6 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let g2 = eisenstein_series(2, prec);
        let mut tmp = HmfGen::new(prec);
        tmp.pow_mut(&g2, 2);
        let mut gens2 = Self::gens2(prec);
        let f2 = gens2.swap_remove(0);
        let f7 = gens2.swap_remove(0);
        let f6 = gens2.swap_remove(0);
        assert_eq!(f2.weight.unwrap().0, 2);
        assert_eq!(f6.weight.unwrap().0, 6);
        assert_eq!(f7.weight.unwrap().0, 7);
        let mut f3 = HmfGen::new(prec);
        div_mut(&mut f3, &f7, &From::from(&tmp));
        vec![f2, f3, f6]
    }

    fn relations() -> Option<Vec<Relation>> {
        let a0 = (MonomFormal { idx: (0, 0, 1) }, Sqrt5Mpz::from_si_g(-1680));
        let a1 = (MonomFormal { idx: (0, 1, 0) }, Sqrt5Mpz::from_si_g(-63000));
        let a2 = (MonomFormal { idx: (1, 0, 0) }, Sqrt5Mpz::from_si_g(1));
        Some(vec![vec![vec![a0], vec![a1], vec![a2]]])
    }
}

impl Structure6 {
    // This is not gens.
    pub fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let i = 6;
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let f4 = rankin_cohen_sqrt5(i, &g2, &g2).unwrap();
        let mut gens = Structure3::gens(prec);
        let mut f6 = gens.remove(0);
        f6.square();
        let f7 = rankin_cohen_sqrt5(i, &g2, &g5).unwrap();
        assert_eq!(f6.weight, Some((6, 18)));
        assert_eq!(f4.weight, Some((4, 16)));
        assert_eq!(f7.weight, Some((7, 19)));
        vec![f4, f6, f7]
    }

    pub fn gens2(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut gens = Self::gens1(prec);
        let mut f2 = HmfGen::new(prec);
        {
            let f4 = &gens[0];
            let g2 = eisenstein_series(2, prec);
            div_mut(&mut f2, f4, &From::from(&g2));
        }
        let f6 = gens.swap_remove(1);
        let mut f7 = gens.swap_remove(1);
        assert_eq!(f6.weight.unwrap().0, 6);
        assert_eq!(f7.weight.unwrap().0, 7);
        let mut tmp = HmfGen::new(prec);
        let g5 = g5_normalized(prec);
        tmp.mul_mut(&f2, &From::from(&g5));
        tmp *= &Sqrt5Mpz::from_ui_g(5);
        f7 *= &Sqrt5Mpz::from_ui_g(8);
        f7 += &tmp;
        vec![f2, f6, f7]
    }
}

pub struct Structure7;

impl Structure for Structure7 {
    fn gens(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        Self::gens3(prec)
    }

    fn relations() -> Option<Vec<Relation>> {
        let from = |x| Sqrt5Mpz::from_si_g(x);
        let m = |x| MonomFormal { idx: x };
        let a0 = from(1);
        let a1 = from(-95256000);
        let a3 = from(1684800);
        let rel0 = vec![
            vec![(m((0, 0, 1)), a0)],
            vec![(m((0, 1, 0)), a1)],
            vec![],
            vec![(m((1, 0, 0)), a3)],
        ];
        let a0 = from(1);
        let a1 = from(-76204800);
        let a2 = from(6264);
        let b2 = from(-9331200);
        let a3 = from(-1244160000);
        let rel1 = vec![
            vec![(m((2, 1, 0)), a0)],
            vec![(m((1, 0, 1)), a1)],
            vec![(m((3, 0, 0)), a2), (m((0, 0, 1)), b2)],
            vec![(m((0, 1, 0)), a3)],
        ];
        Some(vec![rel0, rel1])
    }
}

impl Structure7 {
    pub fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let gens3 = Structure3::gens(prec);
        let gens4 = Structure4::gens(prec);
        let gens2 = Structure2::gens(prec);
        let gens5 = Structure5::gens(prec);
        let f7 = &gens3[0] * &gens4[0];
        let f8 = &gens3[0] * &gens4[1];
        let mut f9 = &gens2[0] * &gens5[0];
        f9 *= 7 as c_ulong;
        let mut g2: HmfGen<Sqrt5Mpz> = From::from(&eisenstein_series(2, prec));
        g2 *= -21384 as c_long;
        f9 += &(&g2 * &f7);
        assert_eq!(f7.weight, Some((7, 21)));
        assert_eq!(f8.weight, Some((8, 22)));
        assert_eq!(f9.weight, Some((9, 23)));
        vec![f7, f8, f9]
    }

    pub fn gens2(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut gens = Self::gens1(prec);
        let mut f4 = HmfGen::new(prec);
        {
            let g5 = From::from(&g5_normalized(prec));
            div_mut(&mut f4, &gens[2], &g5);
        }
        let f7 = gens.swap_remove(0);
        let f8 = gens.swap_remove(1);
        let g2 = eisenstein_series(2, prec);
        let g5 = g5_normalized(prec);
        let f9 = rankin_cohen_sqrt5(7, &(&g2 * &g2), &g5).unwrap();
        assert_eq!(f4.weight, Some((4, 18)));
        assert_eq!(f7.weight, Some((7, 21)));
        assert_eq!(f8.weight, Some((8, 22)));
        vec![f4, f7, f8, f9]
    }

    pub fn gens3(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut gens = Self::gens2(prec);
        let mut f5 = HmfGen::new(prec);
        {
            let f9 = &gens[3];
            let g2 = eisenstein_series(2, prec);
            let g4 = &g2 * &g2;
            let g4: HmfGen<Sqrt5Mpz> = From::from(&g4);
            div_mut(&mut f5, f9, &g4);
        }
        let f4 = gens.remove(0);
        let f7 = gens.remove(0);
        let f8 = gens.remove(0);
        assert_eq!(f4.weight, Some((4, 18)));
        assert_eq!(f5.weight, Some((5, 19)));
        assert_eq!(f7.weight, Some((7, 21)));
        assert_eq!(f8.weight, Some((8, 22)));
        vec![f4, f5, f7, f8]
    }
}


// This does not implement Structrue, but the structure is known and proved
// using Singular.
pub struct Structure8;

impl Structure8 {
    pub fn gens1(prec: usize) -> Vec<HmfGen<Sqrt5Mpz>> {
        let mut gens4 = Structure4::gens(prec);
        let f4 = gens4.remove(0);

        let f8 = &f4 * &f4;

        let mut gens2 = Structure2::gens(prec);
        let mut gens6 = Structure6::gens(prec);
        let f4 = gens2.remove(0);
        let f2 = gens6.remove(0);
        let h3 = gens6.remove(0);

        let f6 = &f2 * &f4;
        let f7 = &h3 * &f4;

        assert_eq!(f6.weight, Some((6, 22)));
        assert_eq!(f7.weight, Some((7, 23)));
        assert_eq!(f8.weight, Some((8, 24)));

        vec![f6, f7, f8]
    }
}

pub struct Structure9;

impl Structure9 {
    pub fn gens1(prec: u64) -> Vec<HmfGen<Sqrt5Mpz>> {
        let gens = three_forms(9, prec as usize).unwrap();
        assert_eq!(gens[0].weight, Some((7, 25)));
        assert_eq!(gens[1].weight, Some((8, 26)));
        assert_eq!(gens[2].weight, Some((11, 29)));
        gens
    }
}

pub struct Structure10;

impl Structure10 {
    pub fn gens1(prec: u64) -> Vec<HmfGen<Sqrt5Mpz>> {
        let gens = three_forms(10, prec as usize).unwrap();
        assert_eq!(gens[0].weight, Some((4, 24)));
        assert_eq!(gens[1].weight, Some((7, 27)));
        assert_eq!(gens[2].weight, Some((11, 31)));
        gens
    }
}

pub fn forms_generated_with_monom(
    k: usize,
    prec: usize,
    gens: &[HmfGen<Sqrt5Mpz>],
) -> Vec<(HmfGen<Sqrt5Mpz>, MonomFormal)> {
    let mut res = Vec::new();
    fn append(
        res: &mut Vec<(HmfGen<Sqrt5Mpz>, MonomFormal)>,
        f: &HmfGen<Sqrt5Mpz>,
        k: usize,
        prec: usize,
    ) {
        let l = f.weight.unwrap().0;
        if k == l {
            let mut tmp = HmfGen::new(prec);
            tmp.set(f);
            res.push((tmp, MonomFormal { idx: (0, 0, 0) }));
        } else if k > l {
            for a in monoms_of_g2_g5_f6(k - f.weight.unwrap().0).into_iter() {
                res.push((f.clone(), a));
            }
        }
    }
    for f in gens.iter() {
        append(&mut res, f, k, prec);
    }
    res
}

pub fn forms_generated_monom(k: usize, gens: &[HmfGen<Sqrt5Mpz>]) -> Vec<Vec<MonomFormal>> {
    fn monoms(f: &HmfGen<Sqrt5Mpz>, k: usize) -> Vec<MonomFormal> {
        let mut v = Vec::new();
        let l = f.weight.unwrap().0;
        if k == l {
            v.push(MonomFormal { idx: (0, 0, 0) });
        } else if k > l {
            for a in monoms_of_g2_g5_f6(k - f.weight.unwrap().0).into_iter() {
                v.push(a);
            }
        }
        v
    }
    gens.iter().map(|f| monoms(f, k)).collect()
}

pub fn forms_generated(k: usize, prec: usize, gens: &[HmfGen<Sqrt5Mpz>]) -> Vec<HmfGen<Sqrt5Mpz>> {
    let forms_w_monm = forms_generated_with_monom(k, prec, gens);
    let mut res = Vec::new();
    for (f, a) in forms_w_monm.into_iter() {
        let g: HmfGen<Sqrt5Mpz> = From::from(&a.into_form(prec));
        res.push(&f * &g);
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

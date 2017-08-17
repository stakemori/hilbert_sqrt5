use gmp::mpz::Mpz;
use eisenstein::{eisenstein_series, f6_normalized};
use theta_chars::g5_normalized;
use elements::HmfGen;
use serde_pickle;
use std::fs::File;
use std::io::Write;
use std::io::Read;
use bignum::RealQuadElement;
use serde;
use std::process::Command;
use std::env;
use bignum::Sqrt5Mpz;


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
pub struct MonomFormal {
    pub idx: (usize, usize, usize),
}

impl MonomFormal {
    pub fn into_form(&self, prec: usize) -> HmfGen<Mpz> {
       monom_g2_g5_f6(prec, self.idx)
    }
    pub fn eval(v: Vec<(Sqrt5Mpz, MonomFormal)>, prec: usize) -> HmfGen<Sqrt5Mpz> {
        let mut res = HmfGen::new(prec);
        for &(ref a, ref monom) in v.iter() {
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
        .map(|&x| {
            MonomFormal {
                idx: x,
            }
        })
        .collect()
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
}

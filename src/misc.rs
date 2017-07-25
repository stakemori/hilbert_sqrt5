use std::ops::{Mul, Add, Sub, Neg, Shr};

/// corresponds to (rt + ir sqrt(5))/2.
#[derive(Debug, Clone)]
pub struct Sqrt5Elt<T> {
    pub rt: T,
    pub ir: T,
}

impl<T> Sqrt5Elt<T> {
    pub fn norm(&self) -> T
    where
        T: Mul<i64, Output = T>
            + Add<Output = T>
            + Mul<T, Output = T>
            + Sub<Output = T>
            + Shr<usize, Output = T>
            + Copy,
    {
        let &Sqrt5Elt { rt: a, ir: b } = self;
        (a * a - b * b * 5) >> 2
    }
}

impl<'a, T> Add<&'a Sqrt5Elt<T>> for Sqrt5Elt<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Sqrt5Elt<T>;
    fn add(self, other: &Sqrt5Elt<T>) -> Sqrt5Elt<T> {
        Sqrt5Elt {
            rt: self.rt + other.rt,
            ir: self.ir + other.ir,
        }
    }
}

impl<'a, T> Sub<&'a Sqrt5Elt<T>> for Sqrt5Elt<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Sqrt5Elt<T>;
    fn sub(self, other: &Sqrt5Elt<T>) -> Sqrt5Elt<T> {
        Sqrt5Elt {
            rt: self.rt - other.rt,
            ir: self.ir - other.ir,
        }
    }
}

impl<'a, T> Neg for &'a Sqrt5Elt<T>
where
    T: Neg<Output = T> + Copy,
{
    type Output = Sqrt5Elt<T>;
    fn neg(self) -> Sqrt5Elt<T> {
        Sqrt5Elt {
            rt: -self.rt,
            ir: -self.ir,
        }
    }
}


impl<'a, 'b, T> Mul<&'a Sqrt5Elt<T>> for &'b Sqrt5Elt<T>
where
    T: Mul<i64, Output = T>
        + Add<Output = T>
        + Mul<T, Output = T>
        + Shr<usize, Output = T>
        + Copy,
{
    type Output = Sqrt5Elt<T>;
    fn mul(self, rhs: &Sqrt5Elt<T>) -> Self::Output {
        let &Sqrt5Elt { rt: a, ir: b } = self;
        let &Sqrt5Elt { rt: c, ir: d } = rhs;
        Sqrt5Elt {
            rt: (a * c + b * d * 5) >> 1,
            ir: (a * d + b * c) >> 1,
        }
    }
}

pub fn prime_sieve(n: usize) -> Vec<usize> {
    let mut vec = Vec::with_capacity(n + 1);
    for _ in 0..(n + 1) {
        vec.push(true);
    }
    let bd = (n as f64).sqrt().floor() as usize;
    for i in 2..(bd + 1) {
        if vec[i] {
            let mut j = i * i;
            while j <= n {
                vec[j] = false;
                j += i;
            }
        }
    }
    let mut res = vec.iter()
        .enumerate()
        .filter(|&(_, &bl)| bl)
        .map(&|(i, _)| i as usize)
        .collect::<Vec<usize>>();
    res.remove(0);
    res.remove(0);
    res
}

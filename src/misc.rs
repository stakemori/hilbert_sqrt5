/// corresponds to (rt + ir sqrt(5))/2.
#[derive(Debug)]
pub struct Sqrt5Elt<T> {
    pub rt: T,
    pub ir: T,
}

pub fn prime_sieve(n: usize) -> Vec<u64> {
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
        .map(&|(i, _)| i as u64)
        .collect::<Vec<u64>>();
    res.remove(0);
    res.remove(0);
    res
}

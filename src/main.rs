use num_bigint::BigInt;
use num_traits::{Zero, One, ToPrimitive};
use num_integer::Integer;  // Add this import

fn main() {
    let b: BigInt = "543134014563836383927251243637".parse().unwrap();
    let m: BigInt = "13637321383345479764947409042367153".parse().unwrap();

    println!("Determine solutions for x^2 ≡ {} (mod {})", b, m);

    match solve_quadratic_congruence(&b, &m) {
        Some(solutions) => {
            if solutions.is_empty() {
                println!("x^2 ≡ {} (mod {}) has no solutions", b, m);
            } else {
                println!("Solutions found: {:?}", solutions);
            }
        }
        None => {
            println!("x^2 ≡ {} (mod {}) has no solutions", b, m);
        }
    }
}

fn solve_quadratic_congruence(b: &BigInt, m: &BigInt) -> Option<Vec<BigInt>> {
    if m.is_even() {
        return solve_mod_power_of_two(b, m);
    }

    if legendre_symbol(b, m) != BigInt::one() {
        return None; // b is not a quadratic residue mod m then no solutions
    }

    // Tonelli-Shanks algorithm for prime moduli
    let mut q: BigInt = m - 1;
    let mut s = BigInt::zero();
    while q.is_even() {
        q /= 2;
        s += 1;
    }

    let mut z = BigInt::from(2);
    while legendre_symbol(&z, m) != -BigInt::one() {
        z += 1;
    }

    let mut c = mod_exp(&z, &q, m);
    let mut r = mod_exp(b, &((q.clone() + 1) / 2), m);
    let mut t = mod_exp(b, &q, m);
    let mut m_s = s.to_usize().unwrap();

    while t != BigInt::one() {
        let mut t2i = t.clone();
        let mut i = 0;
        for j in 1..m_s {
            t2i = (&t2i * &t2i) % m;
            if t2i == BigInt::one() {
                i = j;
                break;
            }
        }

        let b2i = mod_exp(&c, &BigInt::from(2).pow((m_s - i - 1) as u32), m);
        r = (&r * &b2i) % m;
        c = (&b2i * &b2i) % m;
        t = (&t * &c) % m;
        m_s = i;
    }

    Some(vec![r.clone(), m - r])
}

fn solve_mod_power_of_two(b: &BigInt, m: &BigInt) -> Option<Vec<BigInt>> {
    if *m == BigInt::from(2) {
        return Some(vec![b % 2]);
    }

    let k = m.bits() as i64;

    // Solution to x^2 ≡ 1 (mod 2^k)
    if b % 2 == BigInt::one() {
        let mut solutions = vec![BigInt::one(), m - 1];
        for i in 1..k {
            let step = BigInt::from(1) << i;
            let new_solutions: Vec<_> = solutions.iter().map(|x| (x + &step) % m).collect();
            for sol in new_solutions {
                if !solutions.contains(&sol) {
                    solutions.push(sol);
                }
            }
        }
        return Some(solutions);
    }

    // For b even
    None
}

fn legendre_symbol(a: &BigInt, p: &BigInt) -> BigInt {
    let ls = mod_exp(a, &((p - 1) / 2), p);
    if ls == p - 1 {
        -BigInt::one()
    } else {
        ls
    }
}

fn mod_exp(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    let mut result = BigInt::one();
    let mut base = base % modulus;
    let mut exp = exp.clone();

    while exp > BigInt::zero() {
        if &exp % 2 == BigInt::one() {
            result = (&result * &base) % modulus;
        }
        exp >>= 1;
        base = (&base * &base) % modulus;
    }
    result
}





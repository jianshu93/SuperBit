use super::SimHashBits;
use super::sim_hasher::SimHasher;
use std::hash::Hash;
use std::marker::PhantomData;

pub struct SimHash<H, S, const L: usize>
where
    H: SimHasher<T = S>,
    S: SimHashBits,
{
    hasher: H,
    marker: PhantomData<(S, H)>,
}

impl<H, S, const L: usize> SimHash<H, S, L>
where
    H: SimHasher<T = S>,
    S: SimHashBits,
{
    pub fn new(hasher: H) -> Self {
        SimHash {
            hasher,
            marker: PhantomData,
        }
    }

    pub fn create_signature<T, U>(&self, iter: T) -> S
    where
        T: Iterator<Item = U>,
        U: Hash,
    {
        let mut counts = [0i64; L];

        for mut hash in iter.map(|item| self.hasher.hash(&item)) {
            for (_i, count) in counts.iter_mut().enumerate() {
                if hash & S::one() == S::zero() {
                    *count += 1;
                } else {
                    *count -= 1;
                }
                hash >>= 1;
            }
        }

        let mut result = S::zero();
        for i in 0..L {
            if counts[i] > 0 {
                result |= S::one() << i;
            }
        }
        result
    }

    /// Weighted SimHash: each feature contributes ±w instead of ±1.
    /// If the i-th hash bit is 1, add -w; if 0, add +w. Final bit = sign(count).
    pub fn create_signature_weighted<T, U, W>(&self, iter: T) -> S
    where
        T: IntoIterator<Item = (U, W)>,
        U: std::hash::Hash,
        W: Into<f32> + Copy,
    {
        let mut counts = [0f32; L];

        for (item, w0) in iter {
            let mut hash = self.hasher.hash(&item);
            let w = w0.into();

            for i in 0..L {
                if (hash & S::one()) == S::zero() {
                    counts[i] += w;
                } else {
                    counts[i] -= w;
                }
                hash = hash >> 1;
            }
        }

        let mut out = S::zero();
        for i in 0..L {
            if counts[i] > 0.0 {
                out |= S::one() << i;
            }
        }
        out
    }

    pub fn create_centroid<T>(&self, signatures: T) -> S
    where
        T: Iterator<Item = S>,
    {
        let mut counts = [0u64; L];
        let mut len = 0;
        for signature in signatures {
            for i in 0..L {
                if signature >> i & S::one() == S::one() {
                    counts[i] += 1;
                }
            }
            len += 1;
        }
        let mut centroid = S::zero();
        for i in 0..L {
            if counts[i] > len / 2 {
                centroid |= S::one() << i;
            }
        }
        centroid
    }
}

#[cfg(test)]
mod tests {
    use super::super::sim_hasher::{SimSipHasher64, SimSipHasher128, Xxh3Hasher64, Xxh3Hasher128};
    use super::SimHash;
    use super::SimHashBits;

    fn whitespace_split(s: &str) -> impl Iterator<Item = &str> {
        s.split_whitespace()
    }

    static S1: &str = "SimHash is a technique used for detecting near-duplicates or for locality sensitive hashing. It was developed by Moses Charikar and is often used in large-scale applications to reduce the dimensionality of high-dimensional data, making it easier to process";
    static S2: &str = "SimHash is a technique used for detecting near-duplicates or for locality sensitive hashing. It was developed by Moses Charikar and is often utilized in large-scale applications to reduce the dimensionality of high-dimensional data, making it easier to analyze";

    #[test]
    pub fn test_sim_hash_sip64() {
        let sim_hash = SimHash::<SimSipHasher64, u64, 64>::new(SimSipHasher64::new(1, 2));
        let s1 = sim_hash.create_signature(whitespace_split(S1));
        let s2 = sim_hash.create_signature(whitespace_split(S2));
        assert!(s1.hamming_distance(&s2) < 20);
    }

    #[test]
    pub fn test_sim_hash_sip128() {
        let sim_hash = SimHash::<SimSipHasher128, u128, 128>::new(SimSipHasher128::new(1, 2));
        let s1 = sim_hash.create_signature(whitespace_split(S1));
        let s2 = sim_hash.create_signature(whitespace_split(S2));
        assert!(s1.hamming_distance(&s2) < 40);
    }

    #[test]
    pub fn test_sim_hash_xxh3_64() {
        let sim_hash = SimHash::<Xxh3Hasher64, u64, 64>::new(Xxh3Hasher64::new());
        let s1 = sim_hash.create_signature(whitespace_split(S1));
        let s2 = sim_hash.create_signature(whitespace_split(S2));
        assert!(s1.hamming_distance(&s2) < 20);
    }

    #[test]
    pub fn test_sim_hash_xxh3_128() {
        let sim_hash = SimHash::<Xxh3Hasher128, u128, 128>::new(Xxh3Hasher128::new());
        let s1 = sim_hash.create_signature(whitespace_split(S1));
        let s2 = sim_hash.create_signature(whitespace_split(S2));
        assert!(s1.hamming_distance(&s2) < 40);
    }

    // cargo test --release simhash_large -- --nocapture
    #[test]
    fn simhash_large() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        const N: usize = 100_000;
        let mut rng = StdRng::seed_from_u64(42);

        // Make two vectors that differ in ~25% of bits.
        let data1: Vec<u8> = (0..N).map(|_| rng.gen_range(0..2)).collect();
        let mut data2 = data1.clone();
        for i in (0..N).step_by(4) {
            data2[i] = 1 - data2[i];
        }

        let sim_hash = SimHash::<Xxh3Hasher128, u128, 128>::new(Xxh3Hasher128::new());
        let h1 = sim_hash.create_signature((0..N).map(|i| (i as u64, data1[i])));
        let h2 = sim_hash.create_signature((0..N).map(|i| (i as u64, data2[i])));
        let hd = h1.hamming_distance(&h2);

        // ~12–13% of 128 ≈ 15–17 in expectation; allow slack
        assert!(hd < 40);
    }
}

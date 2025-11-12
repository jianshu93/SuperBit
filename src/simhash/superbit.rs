// Super-Bit:
use core::marker::PhantomData;
use std::hash::Hash;
use std::hash::Hasher;
use super::sim_hasher::SimHasher;
use super::SimHashBits;
use xxhash_rust::xxh3::Xxh3;

pub struct SuperBitSimHash<H, S, const L: usize>
where
    H: SimHasher<T = u64>,
    S: SimHashBits,
{
    hasher: H,
    r: usize,                 // superbit depth (block size)
    m: usize,                 // number of blocks = L / r
    q_blocks: Vec<Vec<f32>>,  // each is r×r row-major orthonormal
    seed: u64,
    _phantom: PhantomData<S>, // keep S "used" at the type level
}

impl<H, S, const L: usize> SuperBitSimHash<H, S, L>
where
    H: SimHasher<T = u64>,
    S: SimHashBits,
{
    pub fn new(hasher: H, r: usize, seed: u64) -> Self {
        assert!(r > 0 && L % r == 0, "r must divide L");
        let m = L / r;
        let mut q_blocks = Vec::with_capacity(m);
        for b in 0..m {
            q_blocks.push(Self::make_orthonormal_block(
                r,
                seed ^ ((b as u64) * 0x9E37_79B9),
            ));
        }
        Self {
            hasher,
            r,
            m,
            q_blocks,
            seed,
            _phantom: PhantomData,
        }
    }

    // Classical (stable enough for small r) Gram–Schmidt on a random r×r
    fn make_orthonormal_block(r: usize, seed: u64) -> Vec<f32> {
        let mut mat = vec![0f32; r * r];
        // fill with pseudo-random N(0,1)-ish via hashed Box–Muller
        let mut k = 0u64;
        for i in 0..r {
            for j in 0..r {
                let u1 = Self::u01(seed, k ^ ((i as u64) << 16) ^ (j as u64));
                let u2 = Self::u01(seed, k.wrapping_mul(0x9E37_79B97F4A7C15) ^ 0xBF58_476D);
                k = k.wrapping_add(1);
                let r2 = (-2.0f32 * u1.ln()).sqrt();
                let th = 2.0f32 * std::f32::consts::PI * u2;
                mat[i * r + j] = r2 * th.cos(); // one normal; good enough for r<=32
            }
        }
        // Gram–Schmidt orthogonal
        for j in 0..r {
            // subtract projections
            for p in 0..j {
                let mut dot = 0f32;
                for i in 0..r {
                    dot += mat[i * r + j] * mat[i * r + p];
                }
                for i in 0..r {
                    mat[i * r + j] -= dot * mat[i * r + p];
                }
            }
            // normalize
            let mut n = 0f32;
            for i in 0..r {
                n += mat[i * r + j] * mat[i * r + j];
            }
            let n = n.sqrt().max(1e-12);
            for i in 0..r {
                mat[i * r + j] /= n;
            }
        }
        mat
    }

    #[inline]
    fn u01(seed: u64, x: u64) -> f32 {
        // xxhash-rust API: use `with_seed`
        let mut h = Xxh3::with_seed(seed);
        h.update(&x.to_le_bytes());
        let v = h.finish();
        // map 53-bit mantissa to (0,1); avoid exact 0/1
        ((v >> 11) as f32 + 0.5) * (1.0 / ((1u64 << 53) as f32))
    }

    /// Unweighted items: treat each item as weight 1.0 (like classic SimHash).
    pub fn create_signature<U>(&self, iter: impl Iterator<Item = U>) -> S
    where
        U: Hash,
    {
        self.create_signature_weighted(iter.map(|u| (u, 1.0f32)))
    }

    /// Weighted variant (useful for TF/IDF or MS intensities after L2 norm).
        /// Weighted variant (useful for TF/IDF or MS intensities after L2 norm).
    pub fn create_signature_weighted<U>(&self, iter: impl Iterator<Item = (U, f32)>) -> S
    where
        U: Hash,
    {
        let mut counts = vec![0f32; L];
        // Reuse a single buffer for g to avoid per-block allocations.
        let mut g = vec![0f32; self.r];

        for (item, w) in iter {
            if w == 0.0 { continue; }

            // Stable per-item 64-bit base using the generic hasher H
            let base: u64 = self.hasher.hash(&item);

            // for each block, build g in {+1,-1}^r and accumulate Q_b * g
            for b in 0..self.m {
                let qb = &self.q_blocks[b]; // row-major r×r

                // Seed a tiny PRNG (SplitMix64) once per (item, block)
                let mut s = self.seed
                    ^ base
                    ^ ((b as u64) << 32)
                    ^ 0x9E37_79B9_7F4A_7C15;

                // Fill g[j] ∈ {+1,-1} using SplitMix64
                for j in 0..self.r {
                    s = Self::splitmix64(s);
                    g[j] = if (s >> 63) == 0 { 1.0 } else { -1.0 };
                }

                // u = Q_b * g (r dot-products, row-major is cache-friendly here)
                let off = b * self.r;
                for row in 0..self.r {
                    let mut acc = 0f32;
                    let row_off = row * self.r;
                    // dot(row, g)
                    for col in 0..self.r {
                        acc += qb[row_off + col] * g[col];
                    }
                    counts[off + row] += w * acc;
                }
            }
        }

        // threshold -> bits
        let mut out = S::zero();
        for i in 0..L {
            if counts[i] > 0.0 {
                out |= S::one() << i;
            }
        }
        out
    }

    fn splitmix64(mut x: u64) -> u64 {
        x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = x;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{BitArray, SimHashBits};
    use crate::simhash::sim_hasher::Xxh3Hasher64;
    use rand::rngs::StdRng;
    use rand::Rng;
    use rand::SeedableRng;

    // cargo test --release superbit_simhash_bitarray -- --nocapture
    #[test]
    fn superbit_simhash_bitarray() {
        type Bits = BitArray<16>; // 16×64 = 1024 bits
        const L: usize = 1024;
        const R: usize = 32;      // block size; L % R == 0
        const N: usize = 100_000;

        let mut rng = StdRng::seed_from_u64(12345);
        let data1: Vec<u8> = (0..N).map(|_| rng.gen_range(0..=1)).collect();
        let mut data2 = data1.clone();
        for i in (0..N).step_by(4) {
            data2[i] ^= 1;
        }

        // Ground-truth cosine/angle from the two bit-vectors treated as {0,1}^N
        let (mut dot, mut n1, mut n2) = (0f64, 0f64, 0f64);
        for i in 0..N {
            let x = data1[i] as f64;
            let y = data2[i] as f64;
            dot += x * y;
            n1 += x * x;
            n2 += y * y;
        }
        let cosine = (dot / (n1.sqrt() * n2.sqrt())).clamp(-1.0, 1.0);
        let theta = cosine.acos();
        let p_bit = theta / std::f64::consts::PI; // P(bit differs) for SRP

        let mean = p_bit * L as f64;
        let sigma = (L as f64 * p_bit * (1.0 - p_bit)).sqrt();
        let low = (mean - 3.0 * sigma).floor().max(0.0) as usize;   // 3σ band for robustness
        let high = (mean + 3.0 * sigma).ceil().min(L as f64) as usize;

        let sb = SuperBitSimHash::<Xxh3Hasher64, Bits, L>::new(Xxh3Hasher64::new(), R, 0xDEAD_BEEF);
        let h1 = sb.create_signature_weighted((0..N).map(|i| (i as u64, data1[i] as f32)));
        let h2 = sb.create_signature_weighted((0..N).map(|i| (i as u64, data2[i] as f32)));
        let hd = h1.hamming_distance(&h2);

        eprintln!(
            "SuperBit SimHash (L={L}, R={R}, N={N}): HD = {}, expected ≈ {}–{} (p_bit ≈ {:.3}, mean {:.1}, σ {:.1})",
            hd, low, high, p_bit, mean, sigma
        );
        assert!((low..=high).contains(&hd));
    }
}
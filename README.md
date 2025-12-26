# Batch-Orthogonal Locality-Sensitive Hashing for Angular Similarity

Signed Random Projection Locality Sensitive Hashing (SRP-LSH) or simply SimHash was invented in 2002 by Moses Charikar (1) to quickly approximate Angular similarity and has been the standard since then.

Based on exisiting SRP-LSH, two new SimHash algorithms were proposed (2,3). The key idea is to orthogonalize the random projection vectors in batches. Both the Super-Bit paper and BOLSH paper followed this idea. This crate was created to implement the Super-Bit and BOLSH with optimal performance. Now the bottleneck is pure memory bandwidth because bit Hamming is the limiting step. Several hash libraries, e.g., SipHash and XxHash3 can be used as the underlying hashing algorithm. 

A bit hacked version of original SimHash, called fast SimHash was also provided, which use packed counters to update many bits at once and can be 10 times faster than original SimHash (recommended).


## Key features

	•	SimHash (classic),
	•	FastSimHash (Dynatrace bit-hack port),
	•	SuperBitSimHash (orthogonalized SRP),
	•	and BitArray container + hashers.

## Install
```
# Cargo.toml of your project
[dependencies]
superbit = { git = "https://github.com/jianshu93/SuperBit" }

```

## Usage
Classic SimHash:

```rust
use superbit::{SimHash, SimSipHasher64, SimSipHasher128, Xxh3Hasher64, Xxh3Hasher128};

// Tokenize however you like; any Hash-able item works (e.g., &str, u64, etc.)
fn tokens(s: &str) -> impl Iterator<Item = &str> {
    s.split_whitespace()
}

fn main() {
    // 64-bit SimHash with SipHash-2-4 keyed hasher
    let sh64 = SimHash::<SimSipHasher64, u64, 64>::new(SimSipHasher64::new(0xDEAD_BEEF, 0xC0FF_EE));
    let a64 = sh64.create_signature(tokens("foo bar baz foo"));
    let b64 = sh64.create_signature(tokens("foo bar qux"));
    let hd64 = a64.hamming_distance(&b64);
    println!("HD(64) = {}", hd64);

    // 128-bit SimHash with XXH3
    let sh128 = SimHash::<Xxh3Hasher128, u128, 128>::new(superbit::Xxh3Hasher128::new());
    let a128 = sh128.create_signature(tokens("lorem ipsum dolor sit amet"));
    let b128 = sh128.create_signature(tokens("lorem ipsum dolor sit amet amet"));
    let hd128 = a128.hamming_distance(&b128);
    println!("HD(128) = {}", hd128);

    // (Optional) Estimate cosine from HD:
    // For sign random projections, P[bit differs] = θ/π, where θ = arccos(cosine).
    let l = 128.0;
    let theta = (hd128 as f64 / l) * std::f64::consts::PI;
    let cosine_est = theta.cos();
    println!("cosine ≈ {:.4}", cosine_est);
}

```

Fast SimHash;

```rust
use superbit::{SimHash, SimSipHasher64, SimSipHasher128, Xxh3Hasher64, Xxh3Hasher128};

// Tokenize however you like; any Hash-able item works (e.g., &str, u64, etc.)
fn tokens(s: &str) -> impl Iterator<Item = &str> {
    s.split_whitespace()
}

fn main() {
    // 64-bit SimHash with SipHash-2-4 keyed hasher
    let sh64 = SimHash::<SimSipHasher64, u64, 64>::new(SimSipHasher64::new(0xDEAD_BEEF, 0xC0FF_EE));
    let a64 = sh64.create_signature(tokens("foo bar baz foo"));
    let b64 = sh64.create_signature(tokens("foo bar qux"));
    let hd64 = a64.hamming_distance(&b64);
    println!("HD(64) = {}", hd64);

    // 128-bit SimHash with XXH3
    let sh128 = SimHash::<Xxh3Hasher128, u128, 128>::new(superbit::Xxh3Hasher128::new());
    let a128 = sh128.create_signature(tokens("lorem ipsum dolor sit amet"));
    let b128 = sh128.create_signature(tokens("lorem ipsum dolor sit amet amet"));
    let hd128 = a128.hamming_distance(&b128);
    println!("HD(128) = {}", hd128);

    // (Optional) Estimate cosine from HD:
    // For sign random projections, P[bit differs] = θ/π, where θ = arccos(cosine).
    let l = 128.0;
    let theta = (hd128 as f64 / l) * std::f64::consts::PI;
    let cosine_est = theta.cos();
    println!("cosine ≈ {:.4}", cosine_est);
}

```


SuperBit/BOLSH:

```rust
use superbit::{BitArray, SuperBitSimHash, Xxh3Hasher64};
use std::f64::consts::PI;

fn main() {
    type Bits = BitArray<16>; // 16 × 64 = 1024 bits
    const L: usize = 1024;
    const R: usize = 32;      // block size; must divide L

    // Seed controls the orthonormal basis; fix for reproducibility
    let sb = SuperBitSimHash::<Xxh3Hasher64, Bits, L>::new(Xxh3Hasher64::new(), R, 0xDEAD_BEEF);

    // Unweighted example: stream of tokens (any Hash-able)
    let s1 = "the quick brown fox jumps over the lazy dog";
    let s2 = "the quick brown fox jumps over the very lazy dog";

    let h1 = sb.create_signature(s1.split_whitespace());
    let h2 = sb.create_signature(s2.split_whitespace());

    let hd = h1.hamming_distance(&h2);
    println!("HD(superbit, 1024) = {}", hd);

    // Angle / cosine estimate as with SimHash:
    let theta = (hd as f64 / L as f64) * PI;
    let cosine_est = theta.cos();
    println!("cosine ≈ {:.4}", cosine_est);

    // Weighted example (e.g., TF/IDF):
    // Pass (item, weight) pairs to create_signature_weighted.
    let weighted = vec![
        ("alpha", 0.7f32), ("beta", 0.2), ("gamma", 1.3), ("beta", 0.2),
    ];
    let hw = sb.create_signature_weighted(weighted.into_iter());
    let _ = hw; // use as needed
}

```

This crate was inspired by [gaoya](https://crates.io/crates/gaoya). 



# References
1.Charikar, M.S., 2002, May. Similarity estimation techniques from rounding algorithms. In Proceedings of the thiry-fourth annual ACM symposium on Theory of computing (pp. 380-388).

2.Ji, J., Li, J., Yan, S., Zhang, B. and Tian, Q., 2012. Super-bit locality-sensitive hashing. Advances in neural information processing systems, 25.

3.Ji, J., Yan, S., Li, J., Gao, G., Tian, Q. and Zhang, B., 2014. Batch-orthogonal locality-sensitive hashing for angular similarity. IEEE transactions on pattern analysis and machine intelligence, 36(10), pp.1963-1974.

4.https://www.dynatrace.com/engineering/blog/speeding-up-simhash-by-10x-using-a-bit-hack/

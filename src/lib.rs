
pub mod simhash;

// Re-export common types so users can `use superbit::prelude::*;`
pub mod prelude {
    pub use crate::simhash::{
        BitArray, FastSimHash, SimHash, SimHashBits, SuperBitSimHash,
        SimSipHasher64, SimSipHasher128, Xxh3Hasher64, Xxh3Hasher128,
    };
}

// (Optional) also export at crate root:
pub use simhash::{
    BitArray, FastSimHash, SimHash, SimHashBits, SuperBitSimHash,
    SimSipHasher64, SimSipHasher128, Xxh3Hasher64, Xxh3Hasher128,
};

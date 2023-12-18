use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use core::fmt;

#[derive(Debug, Copy, Clone)]
pub struct SecondSeq( pub u64, pub u8);

impl PartialEq for SecondSeq {
    fn eq(&self, other: &Self) -> bool {
        let mask: u64;

        if other.1 > self.1{
            mask = (1 << (self.1 as u64) *2 ) - 1;
        }else {
            mask = (1 << (other.1 as u64) *2 ) - 1;
        }

        (self.0 & mask) == (other.0 & mask)
    }
}

impl Eq for SecondSeq {}

impl Hash for SecondSeq {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

// Implementing Display trait for SecondSeq
impl fmt::Display for SecondSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(u64: {:b}, u8: {})", self.0, self.1)
    }
}


impl SecondSeq {
    /// calculate the bit flips between two u64 sequences
    pub fn hamming_distance(self, other: &SecondSeq) -> u32 {
        let mask:u64;
        if other.1 > self.1{
            mask = (1 << (self.1 as u64) *2 ) - 1;
        }else {
            mask = (1 << (other.1 as u64) *2 ) - 1;
        }
        ((self.0 & mask) ^ (other.0 & mask)).count_ones()
    }
    pub fn print_second_seq(seq: &SecondSeq) {
        println!("Contents of SecondSeq:");
        println!("u64: {:b}", seq.0);
        println!("{} sig bp", seq.1);
    }

    pub fn fuzzy_match(&self, other:&SecondSeq, max_dist:u32 ) -> bool {

        return self.hamming_distance( other ) <= max_dist.try_into().unwrap()
    }
    pub fn to_le_bytes(&self) -> [u8; std::mem::size_of::<u64>() + std::mem::size_of::<u8>()] {
        let mut bytes = [0; std::mem::size_of::<u64>() + std::mem::size_of::<u8>()];
        bytes[..std::mem::size_of::<u64>()].copy_from_slice(&self.0.to_le_bytes());
        bytes[std::mem::size_of::<u64>()..].copy_from_slice(&self.1.to_le_bytes());
        bytes
    }
    pub fn from_le_bytes(bytes: [u8; 9]) -> Option<Self> {
        if bytes.len() >= std::mem::size_of::<u64>() + std::mem::size_of::<u8>() {
            let mut u64_bytes = [0; std::mem::size_of::<u64>()];
            let mut u8_bytes = [0; std::mem::size_of::<u8>()];

            u64_bytes.copy_from_slice(&bytes[..std::mem::size_of::<u64>()]);
            u8_bytes.copy_from_slice(&bytes[std::mem::size_of::<u64>()..]);

            let u64_val = u64::from_le_bytes(u64_bytes);
            let u8_val = u8::from_le_bytes(u8_bytes);

            Some(SecondSeq(u64_val, u8_val))
        } else {
            None
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_second_seq_equality() {
        let seq1 = SecondSeq(0b10101010, 4);
        let seq2 = SecondSeq(0b10101010, 4);
        let seq3 = SecondSeq(0b11110000, 4);
        let seq4 = SecondSeq(0b1010101011, 6);

        assert_eq!(seq1, seq2, "same sequences"); // Test equal sequences
        assert_ne!(seq1, seq3, "different sequences"); // Test unequal sequences
        assert_eq!(seq1, seq4, "same, but one is longer");
    }

    #[test]
    fn test_second_seq_hashing() {
        let mut map = std::collections::HashMap::new();
        let seq1 = SecondSeq(0b10101010, 4);
        let seq2 = SecondSeq(0b11110000, 4);

        map.insert(seq1, "Value for seq1");

        assert_eq!(map.get(&seq1), Some(&"Value for seq1")); // Test retrieval by equal key
        assert_eq!(map.get(&seq2), None); // Test retrieval by different key
    }
}
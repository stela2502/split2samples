
use std::hash::{Hash, Hasher};
use core::fmt;

#[derive(Debug, Copy, Clone)]
pub struct SecondSeq( pub u64, pub u8);

impl PartialEq for SecondSeq {
    fn eq(&self, other: &Self) -> bool {
        let mask: u64;
        let length:u8;
        if other.1 > self.1{
            length =  self.1
        }else {
            length = other.1;
        }
        if length < 20 {
            return false
        }
        mask = (1 << (length as u64) *2 ) - 1;
        //eprintln!("I'll compare {:b} to {:b}", (self.0 & mask) , (other.0 & mask));
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

    /// the == takes a mininmal matching region into a account
    /// This same function soes not do that.
    pub fn same(&self, other: &Self ) -> bool {
        self.0 == other.0
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
    pub fn table(&self) -> std::collections::HashMap<char, u32> {
        let mut a_cnt = 0;
        let mut c_cnt = 0;
        let mut g_cnt = 0;
        let mut t_cnt = 0;
        let sequence = self.0;

        for i in 0..self.1 {
            let pair = (sequence >> (i * 2)) & 0b11;
            match pair {
                0b00 => a_cnt += 1,
                0b01 => c_cnt += 1,
                0b10 => g_cnt += 1,
                0b11 => t_cnt += 1,
                _ => {} // Handle invalid pairs if needed
            }
        }

        let mut counts = std::collections::HashMap::new();
        counts.insert('A', a_cnt);
        counts.insert('C', c_cnt);
        counts.insert('G', g_cnt);
        counts.insert('T', t_cnt);

        counts
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_second_seq_equality() {
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq3 = SecondSeq(0b1111000010101010101010101010101010101010, 20);
        let seq4 = SecondSeq(0b101010101010101010101010101010101010101010101010, 24);

        assert_eq!(seq1, seq2, "same sequences"); // Test equal sequences
        assert_ne!(seq1, seq3, "different sequences"); // Test unequal sequences
        assert_eq!(seq1, seq4,  "same, but one is longer");
    }

    #[test]
    fn test_second_table() {
        let seq1 = SecondSeq(0b1010101011010101, 32);
        let mut exp = std::collections::HashMap::new();
        exp.insert('A', 24);
        exp.insert('C', 3);
        exp.insert('G', 4);
        exp.insert('T', 1);

        assert_eq!(seq1.table() , exp,  "table did return the right counts");
    }

    #[test]
    fn test_second_seq_hashing() {
        let mut map = std::collections::HashMap::new();
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1111000010101010101010101010101010101010, 20);

        map.insert(seq1, "Value for seq1");

        assert_eq!(map.get(&seq1), Some(&"Value for seq1")); // Test retrieval by equal key
        assert_eq!(map.get(&seq2), None); // Test retrieval by different key
    }
}
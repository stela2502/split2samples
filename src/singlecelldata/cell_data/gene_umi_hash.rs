use std::hash::{Hash, Hasher};
use core::fmt;
use std::cmp::Ord;

#[derive(Debug, Copy, Clone)]
pub struct GeneUmiHash( pub usize, pub u64);
impl PartialEq for GeneUmiHash {
    fn eq(&self, other: &Self) -> bool {
    	(self.0 == other.0) && (self.1 == other.1)
    }
}


impl Ord for GeneUmiHash {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare based on the first field (usize) first
        match self.0.cmp(&other.0) {
            std::cmp::Ordering::Equal => {
                // If the first fields are equal, compare based on the second field (u64)
                self.1.cmp(&other.1)
            }
            ordering => ordering,
        }
    }
}

impl PartialOrd for GeneUmiHash {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for GeneUmiHash {}

impl Hash for GeneUmiHash {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
        self.1.hash(state);
    }
}

// Implementing Display trait for GeneUmiHash
impl fmt::Display for GeneUmiHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GeneUmiHash (gene_id: {}, umi: {})", self.0, self.1 )
    }
}



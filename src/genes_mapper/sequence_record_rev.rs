use std::default::Default;
use core::ops::Add;
use std::ops::AddAssign;
use core::fmt;

#[derive(Debug, Clone)]
pub struct SeqRec<'a> {
    id: &'a [u8],
    seq: &'a [u8],
    qual: &'a [u8],
}

// Implementing Display trait for SeqRec
impl<'a> fmt::Display for SeqRec<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "@{}\n{}\n+\n{}",
            String::from_utf8_lossy(self.id),
            self.as_dna_string(),
            String::from_utf8_lossy(self.qual)
        )
    }
}

impl<'a> PartialEq for SeqRec<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id && self.seq == other.seq && self.qual == other.qual
    }
}

impl<'a> Add for SeqRec<'a> {
    type Output = SeqRec<'a>;

    fn add(self, other: Self) -> Self::Output {
        let mut id = self.id.to_vec();
        id.extend_from_slice(other.id);

        let mut seq = self.seq.to_vec();
        seq.extend_from_slice(other.seq);

        let mut qual = self.qual.to_vec();
        qual.extend_from_slice(other.qual);

        SeqRec {
            id: Box::leak(id.into_boxed_slice()),
            seq: Box::leak(seq.into_boxed_slice()),
            qual: Box::leak(qual.into_boxed_slice()),
        }
    }
}

impl<'a> AddAssign for SeqRec<'a> {
    fn add_assign(&mut self, other: Self) {
        let new_id = [self.id, other.id].concat();
        let new_seq = [self.seq, other.seq].concat();
        let new_qual = [self.qual, other.qual].concat();

        self.id = Box::leak(new_id.into_boxed_slice());
        self.seq = Box::leak(new_seq.into_boxed_slice());
        self.qual = Box::leak(new_qual.into_boxed_slice());
    }
}

impl<'a> Default for SeqRec<'a> {
    fn default() -> Self {
        SeqRec {
            id: &[],
            seq: &[],
            qual: &[],
        }
    }
}

impl<'a> SeqRec<'a> {
    pub fn new(id: &'a [u8], seq: &'a [u8], qual: &'a [u8]) -> Self {
        if seq.len() != qual.len() {
            panic!(
                "I need to get three arrays of the same length! And I got arrays of different lengths here!:\n{:?}, \n{:?},\n{:?}",
                id, seq, qual
            )
        }
        Self { id, seq, qual }
    }

    pub fn id(&self) -> &[u8] {
        self.id
    }

    pub fn seq(&self) -> &[u8] {
        self.seq
    }

    pub fn qual(&self) -> &[u8] {
        self.qual
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn slice(&self, start: usize, len: usize) -> Option<Self> {
        let real_start = if start > self.len() { 0 } else { start };

        let used_len = if real_start + len > self.seq.len() {
            #[cfg(debug_assertions)]
            eprintln!(
                "You try a length of {} from start {} but the sequence is only {} bp long\nI'll give you the max length.\n{}",
                len, real_start, self.seq.len(), String::from_utf8_lossy(self.seq)
            );
            self.seq.len() - start
        } else {
            len
        };

        let end = real_start + used_len;

        let id = format!("{}+{}", String::from_utf8_lossy(self.id), real_start + used_len).into_bytes();
        Some(Self::new(&id, &self.seq[real_start..end], &self.qual[real_start..end]))
    }

    pub fn as_dna_string(&self) -> String {
        String::from_utf8_lossy(self.seq).to_string()
    }

    // Convert the first 32 seq data to a u64
    pub fn to_u64(&self) -> u64 {
        let mut ret = 0_u64;
        for byte in self.seq[0..32.min(self.seq.len())].iter().rev() {
            ret <<= 2;
            ret |= self.enc::<u64>(byte);
        }
        ret
    }

    // Convert the first 16 seq data to a u32
    pub fn to_u32(&self) -> u32 {
        let mut ret = 0_u32;
        for byte in self.seq[0..16.min(self.seq.len())].iter().rev() {
            ret <<= 2;
            ret |= self.enc::<u32>(byte);
        }
        ret
    }

    // Convert the first 8 seq data to a u16
    pub fn to_u16(&self) -> u16 {
        let mut ret = 0_u16;
        for byte in self.seq[0..8.min(self.seq.len())].iter().rev() {
            ret <<= 2;
            ret |= self.enc::<u16>(byte);
        }
        ret
    }

    fn enc<T: From<u8>>(&self, base: &u8) -> T {
        match *base {
            b'A' | b'a' => T::from(0_u8),
            b'C' | b'c' => T::from(1_u8),
            b'G' | b'g' => T::from(2_u8),
            b'T' | b't' => T::from(3_u8),
            b'N' | b'n' => T::from(0_u8), // You might want to handle 'N' differently
            _ => panic!("cannot encode base into 2-bit encoding"),
        }
    }
}

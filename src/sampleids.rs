/// Sampleids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use kmers::naive_impl::Kmer;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;


// and here the data
pub struct SampleIds<'a>{    
    kmers: BTreeMap<u64, u32>,
    kmer_size: usize,
    max_value: u32,
    reads: Vec<u32>
}

// here the functions
impl SampleIds<'_>{
    pub fn new(kmer_size: usize)-> Self {
        let mut kmers = BTreeMap::<u64, u32>::new();
        let max_value:u32 = 0;
        Self {
            kmers,
            kmer_size.unwrap_or(3),
            max_value
        }
    }

    pub fn add(&self, seq: &[u8] ){
        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size ) {
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            self.kmers.insert(km, self.max_value);
        }
        self.max_value += 1;
    }

    pub fn prepare(&self){
        self.reads = vec![0 as u32, self.max_value ];
    }

    pub fn get(&self, seq: &[u8], jump: usize) -> Result< u32, &str>{
        let jump_loc = jump.unwrap_or(1);
        let res = vec![0; self.max_value];
        let mut max_value = 0;
        let mut ret:u32 = 0;
        kmer_vec = needletail::kmer::Kmers::new(seq, self.kmer_size );
        while id < kmer_vec.len() {
            let km = Kmer::from(kmer_vec[id]).into_u64();
            match self.kmers.get(&km){
                Some(c1) => {
                        //println!("to_cellid the c1 {}", c1 );
                        res[c1] += 1;
                    },
                None => (), 
            }
            id += jump_loc;
        }

        for i in 0..self.max_value{
            if res[i] > max_value{
                ret = i;
                max_value = res[i];
            }
        }

        let mut z = 0;
        if max_value > 2 {
            for i in 0..res.len(){
                if res[i] == max_value {
                    ret = i as u32;
                    z += 1;
                }
            }
        }else {
            return Err::<u32, &str>( "NoMatch");
        }
        if z == 1{
            println!("I got a match to entry {} with {} matches", ret, max_value);
        }else {
            println!("Error: I got {} ids with {} matches", z);
            return Err::<u32, &str>( "MultiMatch");
        }
        self.reads[ret] += 1;
        Ok(ret);
    }
}


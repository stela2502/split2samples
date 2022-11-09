/// Sampleids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use kmers::naive_impl::Kmer;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;

#[derive(Debug)]
pub struct Reads {
        pub total: i32,
        pub part: i32, 
}


impl Reads{
        pub fn new( )-> Self {
            let total = 0;
            let part = 0;
            Self {
                total,
                part,
            }
        }
        pub fn add_total(&mut self) {
            self.total += 1;
        }
        pub fn add_part( &mut self ){
            self.part +=1;
        }
}

// and here the data
pub struct SampleIds{    
    kmers: BTreeMap<u64, u32>,
    kmer_size: usize,
    max_value: u32,
    pub read: BTreeMap<u32, Reads>
}

// here the functions
impl SampleIds{
    pub fn new(kmer_size: usize)-> Self {
        let kmers = BTreeMap::<u64, u32>::new();
        let max_value:u32 = 0;
        let read = BTreeMap::<u32, Reads>::new();
        Self {
            kmers,
            kmer_size: kmer_size,
            max_value,
            read
        }
    }

    pub fn read_get(&self, i: u32) -> Option<&Reads> {
        return self.read.get( &i )
    }

    // pub fn read_get_mut(&mut self, i: u32) -> Option<&mut Reads> {
    //     return self.read.get_mut( &i )
    // }

    pub fn add(&mut self, seq: &[u8] ){
        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            self.kmers.insert(km, self.max_value);
        }
        self.read.insert( self.max_value, Reads::new() );
        self.max_value += 1;
    }


    pub fn get(&mut self, seq: &[u8], jump: usize, min_value: usize, min_z: usize ) -> Result< u32, &str>{
        
        for (_key, value) in self.read.iter_mut() {
            value.part = 0;
        }
        let mut max_value = 0;
        let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut id = 0;
        while id < kmer_vec.len() {
            let km = kmer_vec[id];
            //println!("SampleIds::get - checking this sequence: {} and the at pos {}", km, id );
            match self.kmers.get_mut(&km){
                Some(c1) => {
                        //println!("to_cellid the c1 {}", c1 );
                        match self.read.get_mut(c1){
                            Some(read) => {
                                read.add_part();
                            }
                            None => ()
                        };
                    },
                None => (), 
            }
            id += jump;
        }

        for i in 0..self.max_value{
            match self.read.get(&i){
                Some(value) => {
                    if value.part > max_value{
                        max_value = value.part;
                    }
                },
                None => (),
            };
        }

        let mut z = 0;
        if max_value > min_value as i32 {
            for i in 0..self.max_value{
                match self.read.get(&i){
                    Some(value) => {
                        if value.part == max_value{
                            ret = i;
                            z += 1;
                        }
                    },
                    None => (),
                };
            }
        }else {
            return Err::<u32, &str>( "Samples NoMatch - likely data read");
        }

        if z > min_z as i32{
            println!("Error: I got {} ids with {} matches", z, max_value);
            return Err::<u32, &str>( "MultiMatch");
        }
        match self.read.get_mut( &ret ) {
            Some(val) => val.add_total(),
            None => (),
        };
        Ok(ret)
    }
}


fn fill_kmer_vec<'a>(seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
   kmer_vec.clear();
   let mut bad = 0;
   for km in seq {
        // I would like to add a try catch here - possibly the '?' works?
        // if this can not be converted it does not even make sense to keep this info
        for nuc in km{
            if *nuc ==b'N'{
                bad = 1;
            }
        }
        if bad == 0{
            // let s = str::from_utf8(km);
            // println!( "this is the lib: {:?}",  s );
            kmer_vec.push(Kmer::from(km).into_u64());
        }
   }
}


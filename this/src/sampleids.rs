/// Sampleids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use kmers::naive_impl::Kmer;
use std::collections::HashSet;

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
    pub read: BTreeMap<u32, Reads>,
    bad_entries: HashSet<u64>
}

// here the functions
impl SampleIds{
    pub fn new(kmer_size: usize)-> Self {
        let kmers = BTreeMap::<u64, u32>::new();
        let max_value:u32 = 0;
        let read = BTreeMap::<u32, Reads>::new();
        let bad_entries = HashSet::<u64>::new();

        Self {
            kmers,
            kmer_size: kmer_size,
            max_value,
            read,
            bad_entries
        }
    }

    /// adds the BD Rhapsody sample primers to the object
    pub fn init_rhapsody(&mut self, specie:&str ){

        if  specie.eq("human") {
            // get all the human sample IDs into this.
            self.add( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG" );
            self.add( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG" );
            self.add( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT" );
            self.add( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT" );
            self.add( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG" );
            self.add( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC" );
            self.add( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT" );
            self.add( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT" );
            self.add( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG" );
            self.add( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC" );
            self.add( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC" );
            self.add( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG" );

        }
        else if specie.eq("mouse") {
            // and the mouse ones
            self.add( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG" );
            self.add( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC" );
            self.add( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC" );
            self.add( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT" );
            self.add( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT" );
            self.add( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC" );
            self.add( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC" );
            self.add( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG" );
            self.add( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC" );
            self.add( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC" );
            self.add( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT" );
            self.add( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC" );

        } else {
            println!("Sorry, but I have no primers for species {}", specie);
            std::process::exit(1)
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
            for nuc in kmer {
                if *nuc ==b'N'{
                    continue;
                }
            }
            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            if self.bad_entries.contains( &km ){
                continue
            }
            if self.kmers.contains_key ( &km ){
                self.bad_entries.insert( km.clone() );
                self.kmers.remove( &km );
            }else {
                self.kmers.insert(km, self.max_value);
            }        
        }
        self.read.insert( self.max_value, Reads::new() );
        self.max_value += 1;
    }


    pub fn get(&mut self, seq: &[u8], jump: usize, start: usize ) -> Result< u32, &str>{
        
        for (_key, value) in self.read.iter_mut() {
            value.part = 0;
        }
        let min_value = 2;
        let min_z = 1;
        let mut max_value = 0;
        let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut id = start;
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

// #[cfg(test)]
// mod tests {
//     #[test]    
//     fn test_genes_names_ids () {
//         let mut genes = super::parse_bc_map( "testData/HTOs.csv", 9 );

//         let exp = vec![0,1,2,3,4,5 ];
//         let mut data = Vec::<usize>::with_capacity(7);
//         for ( _name, id ) in &genes.names{
//             eprintln!( "{}", id);
//             data.push(*id);
//         }
//         assert_eq!( exp, data);
//     }
//     fn test_genes_get (){

//         let mut genes = ;
//         // Hope I get the correct id:
//         let mut data = b"CTTGCCGCATGTCAT";
//         let mut val = genes.get( data );
//         assert_eq!( Some(2), val );
//         data = b"ACCCACCAGTAAGAC";
//         val = genes.get( data );
//         assert_eq!( Some(0), val );
//         data = b"NNNGCCGCATGTCAN" ;
//         val = genes.get( data );
//         assert_eq!( Some(2), val );

//         val = genes.get( b"NNNGCCNCATGTCAN" );
//         let val2:Option<usize> = None;
//         assert_eq!( val2, val );
//     }
// }
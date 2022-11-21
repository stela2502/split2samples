/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use std::collections::HashSet;

use kmers::naive_impl::Kmer;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;

// #[derive(Debug)]
// pub struct Info {
// //    pub cells: CellIds10x,
//     pub id: u64,
//     pub name: std::string::String
// }


// impl  Info{
//         pub fn new( id: u64, name: std::string::String )-> Self {
//             let loc_name = name.clone();
//             Self {
//                 id,
//                 name: loc_name,
//             }
//         }
// }

/// GeneIds harbors the antibody tags
/// these sequences have (to my knowmledge) 15 bp os length
/// but that can be different, too. 
/// Hence we store here 
/// kmers       : the search object
/// seq_len     : the length of the sequences (10x oversequences them)
/// kmer_size   : the length of the kmers
/// names       : a hashset for the gene names
/// bad_entries : a hash to save bad entries (repetetive ones)
pub struct GeneIds{    
    pub kmers: BTreeMap<u64, usize>, // the search map with kamer u64 reps.
    pub seq_len: usize, // size of the sequence that has been split into kmers
    kmer_size: usize, // size of the kmers
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    bad_entries: HashSet<u64>, // non unique u64 values that will not be recoreded.
    max_id: usize // hope I get the ids right this way...
}

// here the functions
impl GeneIds{
    pub fn new(kmer_size: usize )-> Self {
        let kmers = BTreeMap::<u64, usize>::new();
        let names = BTreeMap::<std::string::String, usize>::new();
        let bad_entries = HashSet::<u64>::new();
        let seq_len = 0;
        let max_id = 0;
        Self {
            kmers,
            seq_len,
            kmer_size: kmer_size,
            names,
            bad_entries,
            max_id
        }
    }

    pub fn add(&mut self, seq: &[u8], name: std::string::String ){
        
        if seq.len() > self.seq_len{
            self.seq_len = seq.len() 
        }
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
                //let info = Info::new(km, name.clone() );
                if ! self.names.contains_key( &name ){
                    self.names.insert( name.clone(), self.max_id );
                    self.max_id += 1;
                }
                //println!("I insert a kmer for id {}", self.max_id-1 );
                self.kmers.insert(km, self.max_id-1 );
            }
        }
    }

    pub fn get(&mut self, seq: &[u8] ) -> Option< usize >{
        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut ret:Option<usize> = None;

        if kmer_vec.len() == 0 {
            //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
            return ret
        }  
        let mut sums = vec![0 ;self.names.len()];
        let mut max = 0;

        for km in kmer_vec{
            //println!( "searching for kmer {}", km);
            match self.kmers.get(&km){
                Some(c1) => {
                    //println!("And got a match: {}", c1);
                    sums[*c1] += 1;
                    if max < sums[*c1]{
                        //println!("the new max is {}", max);
                        max =  sums[*c1];
                        if max > 1{
                            break // 2 unique hits should be enough
                        }
                    };
                }
                None => ()
            };
        }
        for i in 0..sums.len(){
            if sums[i] == max{
                //println!("Now the ret hould have the value {} resp {:?}", i, Some(i));
                ret = Some(i);
                break;
            }
        }
        //println!("return geneid {:?}", ret);
        return ret
    }

    // pub fn to_ids( &self,  ret:&mut Vec<Info> )  {
    //     ret.clear();
    //     for (_i, obj) in &self.kmers {
    //         ret.push(*obj );
    //     }
    // }

    pub fn to_header( &self ) -> std::string::String {
        let mut ret= Vec::<std::string::String>::with_capacity( self.names.len() +2 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for (obj, _id) in &self.names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( format!( "{}", obj)) ;
        }
        ret.push("Most likely name".to_string());
        ret.push("Faction total".to_string());
        return "CellID\t".to_owned()+&ret.join("\t")
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
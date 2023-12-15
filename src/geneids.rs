/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use std::collections::HashSet;

use kmers::naive_impl::Kmer;

use crate::ofiles::Ofilesr;
use crate::ifiles::Ifilesr;

use std::path::Path;
use std::io::Write;
use std::fs;

use std::io::BufRead;
use std::io::Read;

pub use crate::traits::Index;

/// GeneIds harbors the antibody tags, the mRNA tags and whatever tags more you search the R2 for
/// but that can be different, too. 
/// Hence we store here 
/// kmers       : the search object
/// seq_len     : the length of the sequences (10x oversequences them)
/// kmer_size   : the length of the kmers
/// names       : a hashset for the gene names
/// bad_entries : a hash to save bad entries (repetetive ones)
/// and a lot of private ease of live id to name or vice versa BTreeMaps

#[derive(Debug,PartialEq)]
pub struct GeneIds{    
    pub kmers: BTreeMap<u64, usize>, // the search map with kamer u64 reps.
    pub seq_len: usize, // size of the sequence that has been split into kmers
    kmer_size: usize, // size of the kmers
    names_store: BTreeMap< usize, std::string::String>,
    kmer_store: BTreeMap< u64, std::string::String>,
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    bad_entries: HashSet<u64>, // non unique u64 values that will not be recoreded.
    pub good_entries: HashSet<u64>, // upon export as sparse matrix it has to be checked if a gene has a value
    pub max_id: usize,// hope I get the ids right this way...
    unchecked: bool, //if add_unchecked was used results should always be as get_unchecked,
    pub last_count: usize,
    pub last_kmers: Vec<String>,
}


impl Index for GeneIds {

    fn add(&mut self, seq: &Vec<u8>, name: String, _class_ids: Vec<String>) -> usize {
        
        if seq.len() > self.seq_len{
            self.seq_len = seq.len() 
        }

        let mut checker = BTreeMap::<u8, usize>::new();
        let mut total = 0;

        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            checker.clear();
            //check for too simple kmers
            for nuc in kmer {

                match checker.get_mut( nuc ){
                    Some( map ) => *map += 1,
                    None => {
                        if checker.insert( *nuc, 1).is_some(){};
                    }
                };
                if *nuc ==b'N'{
                    continue;
                }
            }
            if checker.len() < 3{
                //println!( "kmer for gene {} is too simple/not enough diff nucs): {:?}", name, std::str::from_utf8(kmer)  );
                continue;
            }
            for ( _key, value ) in checker.iter(){
                if *value as f32 / self.kmer_size as f32 > 0.6 {
                    //println!( "kmer for gene {} is too simple/too many nucs same: {:?}", name, std::str::from_utf8(kmer)  );
                    continue;
                } 
            }

            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();

            if self.bad_entries.contains( &km ){
                continue
            }

            if let std::collections::btree_map::Entry::Vacant(e) = self.kmers.entry(km) {
                //let info = Info::new(km, name.clone() );
                if ! self.names.contains_key( &name ){
                    self.names.insert( name.clone(), self.max_id );
                    self.names_store.insert( self.max_id, name.clone() );
                    self.max_id += 1;
                }
                //println!("I insert a kmer {} for name {}", std::str::from_utf8(kmer).unwrap().to_string(), name );
                self.kmer_store.insert( km , std::str::from_utf8(kmer).unwrap().to_string() );
                e.insert(self.max_id-1);
                total +=1;
            } 
            else {
                self.bad_entries.insert( km );
                self.kmers.remove( &km );
            }
        }
        if total < 3{
            eprintln!( "Sequence for gene {name} has too little OK kmers to match!");
            return 0;
        }
        total
        //println!( "{} kmers for gene {}", total, name );

    }

    fn get(&self, seq: &[u8]) -> Option<usize> {        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        //let report = self.get_id("2810417H13Rik".to_string());

        let mut ret:Option<usize> = None;

        if kmer_vec.is_empty() {
            //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
            return ret
        }  
        let mut sums = vec![0 ;self.names.len()];
        let mut max = 0;
        //self.last_kmers.clear();

        for km in kmer_vec{
            //println!( "searching for kmer {}", km);
            if let Some(c1) = self.kmers.get(&km) {
                sums[*c1] += 1;
                if max < sums[*c1]{
                    max =  sums[*c1];
                }
                //self.last_kmers.push(  self.kmer_store.get( &km).unwrap().to_string() );
            }

        }
        //self.last_count = max;
        if max >2 {
            for (i, sum) in sums.iter().enumerate() {
                if *sum == max{
                    //println!("Now the ret hould have the value {} resp {:?}", i, Some(i));
                    if ret.is_some() {
                        //eprintln!("I have two genes matching with max value of {}: {:?} and {}", max, ret, i);
                        ret =None;
                        return ret
                    }
                    ret = Some(i);
                }
            }
        }
        //println!("return geneid {:?}", ret);
        ret
    }

    fn get_id(&self, name: String) -> usize {
        let id = match self.names.get( &name ) {
            Some( id ) => id,
            None => panic!("Gene {name} not defined in the GeneID object"),
        };
        *id
    }

    fn print(&self) {
        println!("I have {} kmers and {} genes", self.kmers.len(), self.names.len() );
    }

    fn change_start_id(&mut self, new_start: usize) {
        self.max_id = new_start;
    }

    fn names_len(&self) -> usize {
        self.names.len()
    }

    fn names(&self) -> Vec<String> {
        self.names.keys().cloned().collect()
    }

    fn to_header_n(&self, names: &Vec<String>) -> String {
        let mut ret= Vec::<std::string::String>::with_capacity( names.len() +2 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for name in names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( name.to_string() ) ;
        }
        ret.push("AsignedSampleName".to_string());
        ret.push("FractionTotal".to_string());
        "CellID\t".to_owned()+&ret.join("\t")
    }

    fn max_id(&self) -> usize {
        self.max_id
    }

    fn write_index( &mut self, path: String ) -> Result< (), &str>{
        let rs = Path::new( &path ).exists();

        if ! rs {
            match fs::create_dir ( path.clone() ){
                Ok(_file) => (),
                Err(err) => {
                     eprintln!("Error?: {err:#?}");
                 }
            };
        }

        // remove the old files if they exist:
        let fpath = Path::new(&path ) ;
        if fs::remove_file(fpath.join("index.1.Index.gz") ).is_ok(){};
        if fs::remove_file(fpath.join("index.1.gene.txt.gz") ).is_ok(){};

        let mut ofile = Ofilesr::new( 1, "index", "Index", "gene.txt",  &path );

        let seq_len:u64 = self.kmer_size as u64;


        match ofile.buff1.write( &seq_len.to_le_bytes() ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("seq_len could not be written"),
        };

        for ( kmer, gene_id) in &self.kmers {
            match ofile.buff1.write( &kmer.to_le_bytes() ){
                Ok(_) => println!("{} -> {:?} -> {:?}",kmer, &kmer.to_le_bytes(),  &u64::from_le_bytes(kmer.to_le_bytes())  ) ,
                Err(_err) => return Err::<(), &str>("kmer could not be written"),
            };

            match self.names_store.get( &gene_id ){
                Some(name) => {
                    write!(ofile.buff2,"{name}\n" ).unwrap();
                },
                None => panic!("Gene id {gene_id} not found in this object!"),
            };
        }


        ofile.close();
        Ok(())
    }

    fn load_index( &mut self, path: String ) -> Result< (), &str>{

        let f1 = Path::new( &path ).join("index.1.Index");
        let f2 = Path::new( &path ).join("index.1.gene.txt");
        if ! f1.exists() {
            panic!("Index file {} does not exist", f1.to_str().unwrap() );
        }
        if ! f2.exists() {
            panic!("kmer names file {} does not exist", f2.to_str().unwrap() );
        }

        let mut ifile =Ifilesr::new( 1, "index", "Index", "gene.txt",  &path  );
        let mut id = 0;

        let mut buff = [0_u8 ;8 ];

        ifile.buff1.read_exact(&mut buff).unwrap();
        //println! ("I got the buff {buff:?}");
        self.kmer_size = u64::from_le_bytes( buff ) as usize;

        for name in ifile.buff2.lines() {
            //println!("I read this name: {name:?}");
            //let mut buff = [0_u8 ;8 ];
            ifile.buff1.read_exact(&mut buff).unwrap();

            let km = &u64::from_le_bytes( buff );
            
            //println!("I try to insert km {km} and gene name '{name}' ({len:?}) ");
            id += 1;
            let mut data:String = "".to_string() ; // needs to be a constant...
            match name{
                Ok(na) => {
                    if let std::collections::btree_map::Entry::Vacant(e) = self.kmers.entry(*km) {
                        //let info = Info::new(km, name.clone() );
                        if ! self.names.contains_key( &na ){
                            self.names.insert( na.clone(), self.max_id );
                            self.names_store.insert( self.max_id, na.clone() );
                            self.max_id += 1;
                        }
                        let gene_id = self.names.get( &na ).unwrap();
                        //self.kmer_store.insert( *km , "not recovered".to_string() );
                        //println!("These are the bytes I need to work on: {:?}",&km.to_le_bytes());
                        Self::u64_to_str(self.kmer_size, km, &mut data );
                        println!( "({})Can I get that {km} to string? {:?}", self.seq_len, data );
                        
                        self.kmer_store.insert( *km , data.to_string() );
                        e.insert( *gene_id );
                    } 
                    else {
                        self.bad_entries.insert( *km );
                        self.kmers.remove( km );
                    }
                },
                Err(err) => panic!("This is not working: {err:?}"),
            };

            //println!(" {:?} -> {:?}", buff,  &u64::from_le_bytes(buff)  ) ;
        }
        if id == 0 {
            panic!("No data read!");
        }

        //self.print();
        Ok(())
    }
}
// here the functions

impl GeneIds{

    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new(kmer_len: usize, _allocate:usize )-> Self {
        let kmers = BTreeMap::<u64, usize>::new();
        let names = BTreeMap::<std::string::String, usize>::new();
        let names_store = BTreeMap::< usize, std::string::String>::new();
        let kmer_store = BTreeMap::< u64, std::string::String>::new();
        let names4sparse = BTreeMap::<std::string::String, usize>::new();
        let bad_entries = HashSet::<u64>::new();
        let good_entries = HashSet::<u64>::new();
        let last_kmers = Vec::with_capacity(60);
        Self {
            kmers,
            seq_len: 0, // Example initial value, change as needed
            kmer_size: kmer_len, // Using the parameter to set the kmer_size field
            names_store,
            kmer_store,
            names,
            names4sparse,
            bad_entries,
            good_entries,
            max_id: 0, // Example initial value, change as needed
            unchecked: false, // Example initial value, change as needed
            last_count: 0, // Example initial value, change as needed
            last_kmers,
        }
    }

    pub fn add_unchecked(&mut self, seq: &[u8], name: std::string::String ){
        
        if seq.len() > self.seq_len{
            self.seq_len = seq.len() 
        }
        self.unchecked = true;
        let mut checker = BTreeMap::<u8, usize>::new();
        let mut total = 0;

        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            checker.clear();
            for nuc in kmer {
                match checker.get_mut( nuc ){
                    Some( map ) => *map += 1,
                    None => {
                        checker.insert( *nuc, 1);
                    }
                }
                if *nuc ==b'N'{
                    continue;
                }
            }
            if checker.len() < 3{
                //println!( "kmer for gene {} is too simple/not enough diff nucs): {:?}", name, std::str::from_utf8(kmer)  );
                continue;
            }
            for ( _key, value ) in checker.iter(){
                if *value as f32 / self.kmer_size as f32 > 0.6 {
                    //println!( "kmer for gene {} is too simple/too many nucs same: {:?}", name, std::str::from_utf8(kmer)  );
                    continue;
                } 
            }

            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            if self.bad_entries.contains( &km ){
                continue
            }

            if let std::collections::btree_map::Entry::Vacant(e) = self.kmers.entry(km) {
                //let info = Info::new(km, name.clone() );
                if ! self.names.contains_key( &name ){
                    self.names.insert( name.clone(), self.max_id );
                    self.names_store.insert( self.max_id, name.clone() );
                    self.max_id += 1;
                }
                //println!("I insert a kmer {} for name {}", std::str::from_utf8(kmer).unwrap().to_string(), name );
                self.kmer_store.insert( km , std::str::from_utf8(kmer).unwrap().to_string() );
                e.insert(self.max_id-1);
                total +=1;
            } 
            else {
                self.bad_entries.insert( km );
                self.kmers.remove( &km );
            }
            //println!("I insert a kmer for id {}", self.max_id-1 );
        }
        if total < 3{
            eprintln!( "Sequence for gene {name} has too little OK kmers to match!");
        }
        //println!( "{} kmers for gene {}", total, name );

    }

    

    pub fn get_name( &self, id:usize) -> String{
        let name = match self.names_store.get( &id ){
            Some( na ) => na,
            None => panic!("GeneID {id} not defined in the GeneID object"),
        };
        name.to_string()
    }

    

    pub fn get_unchecked(&mut self, seq: &[u8] ) -> Option< usize >{
        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut ret:Option<usize> = None;

        if kmer_vec.is_empty() {
            //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
            return ret
        }  
        let mut sums = vec![0 ;self.names.len()];
        let mut max = 0;
        self.last_kmers.clear();

        for km in kmer_vec{
            //println!( "searching for kmer {}", km);
            if let Some(c1) = self.kmers.get(&km) {
                sums[*c1] += 1;
                if max < sums[*c1]{
                    max =  sums[*c1];
                }
                self.last_kmers.push(  self.kmer_store.get( &km).unwrap().to_string() );
            }
        }
        self.last_count = max;
        if max >3 {
            for (i, sum) in sums.iter().enumerate(){
                if *sum == max{
                    if ret.is_some() {
                        //eprintln!("I have two genes matching with max value of {}: {:?} and {}", max, ret, i);
                        ret =None;
                        return ret
                    }
                    //println!("max = {} -> now the ret hould have the value {:?}", max, Some(i));
                    ret = Some(i);
                    //break;
                }
            }
        }
        // else {
        //     eprintln!("The max was {} => no (good) gene found", max);
        // }
        //println!("return geneid {:?}", ret);
        ret
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
        for obj in self.names.keys() {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push(  obj.to_string() ) ;
        }
        ret.push("Most likely name".to_string());
        ret.push("Faction total".to_string());
        "CellID\t".to_owned()+&ret.join("\t")
    }

    

    

    /*
    if that works here is the CatGPG discussion  that taught me that:
    rust how to access bits of a u8 0,1 and 2,3 and so forth?

    To access bits of a u8, use the bitwise operators & (AND) and | (OR). For example, to access the 0th and 1st bits of a u8, you can use the expression:

    u8 & 0x03

    What would be the result of this operation?

    The result of this operation would be the 0th and 1st bits of the u8.

    Can that be applioed multiple times to get the next 2 ans so on?

    Yes, you can use this operation multiple times to get the next 2 bits and so on. For example, to access the 2nd and 3rd bits of a u8, you could use the expression:

    u8 & 0x0C

    Can I get that for the rest, too?

    Yes, you can use this operation to get the rest of the bits as well. For example, to access the 4th and 5th bits of a u8, you could use the expression:

    u8 & 0x30

    and one more?

    Yes, to access the 6th and 7th bits of a u8, you could use the expression:

    u8 & 0xC0

    That is it then - there are no more bits in an u8 - or?

    Yes, that is it - there are no more bits in an u8.

    EPIC! Thank you!

    You're welcome!
    */
    fn u64_to_str (kmer_size:usize, km:&u64,  data:&mut String ){

        let mut i = 0;
        for u8_rep in km.to_le_bytes(){
            if i < kmer_size{
                match u8_rep & 0x03{
                    0b00000000 => *data +="A",
                    0b00000001 => *data +="C",
                    0b00000010 => *data +="G",
                    0b00000011 => *data +="T",
                    _ => *data +="N",
                };
                i += 1;
            }
            if i < kmer_size{
                match u8_rep & 0x0C{
                    0b00000000 => *data +="A",
                    0b00000100 => *data +="C",
                    0b00001000 => *data +="G",
                    0b00001100 => *data +="T",
                    _ => *data +="N",
                };
                i += 1;
            }
            if i < kmer_size{
                match u8_rep & 0x30{
                    0b00000000 => *data +="A",
                    0b00010000 => *data +="C",
                    0b00100000 => *data +="G",
                    0b00110000 => *data +="T",
                    _ => *data +="N",
                };
                i += 1;
            }
            if i < kmer_size{
                match u8_rep & 0xC0{
                    0b00000000 => *data +="A",
                    0b01000000 => *data +="C",
                    0b10000000 => *data +="G",
                    0b11000000 => *data +="T",
                    _ => *data +="N",
                };
                i += 1;
            }
            if i == kmer_size{
                break;
            }

        //println!( "({})Can I get that {u8_rep} to u2 ({i})? {:?}, {:?}, {:?}, {:?}, {:?}", self.seq_len, u8_rep & 0x03,  u8_rep & 0x0C, u8_rep & 0x30 , u8_rep & 0xC0, data );
        }
    }

}




fn fill_kmer_vec(seq: needletail::kmer::Kmers, kmer_vec: &mut Vec<u64>) {
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


#[cfg(test)]
mod tests {

    use crate::geneids::GeneIds;
    use std::path::Path;
    use crate::traits::Index;

     #[test]
    fn check_geneids() {
        let mut genes = GeneIds::new( 7, 1 );

        let mut geneid = 0;
        
        genes.add( &b"AGCTGCTAGCCGATATT".to_vec(), "Gene1".to_string(), Vec::<String>::new() );
        genes.names4sparse.insert( "Gene1".to_string(), geneid );
        genes.add( &b"CTGTGTAGATACTATAGATAA".to_vec(), "Gene2".to_string(), Vec::<String>::new() );
        genes.names4sparse.insert( "Gene1".to_string(), geneid );
        // the next two should not be in the output
        genes.add( &b"CGCGATCGGATAGCTAGATAGG".to_vec(), "Gene3".to_string(), Vec::<String>::new() );
        genes.add( &b"CATACAACTACGATCGAATCG".to_vec(), "Gene4".to_string(), Vec::<String>::new() );

        geneid = genes.get_id( "Gene1".to_string() );
        assert_eq!( geneid,  0 ); 

        geneid = genes.get_id( "Gene3".to_string() );
        assert_eq!( geneid,  2 ); 

        genes.write_index( Path::new("testData/output_index_test").to_str().unwrap().to_string() ).unwrap();

        let mut genes2 = GeneIds::new( 7, 1 );

        genes2.load_index( Path::new("testData/output_index_test").to_str().unwrap().to_string() ).unwrap();

        assert_eq!( Vec::from_iter( genes.kmers.keys().cloned() ),  Vec::from_iter( genes2.kmers.keys().cloned()) );
        assert_eq!( Vec::from_iter( genes.kmer_store.keys().cloned()),  Vec::from_iter( genes2.kmer_store.keys().cloned()) );
        assert_eq!( genes.kmer_size, genes2.kmer_size )
    }

    #[test]
    fn test_u64_to_str(){

        let num:u64 = 15561;
        let mut data:String = "".to_string();
        GeneIds::u64_to_str( 7, &num, &mut data );
        assert_eq!( data, "CGATATT".to_string() )
    } 
}

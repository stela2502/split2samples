/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use needletail::kmer::Kmers;
use kmers::naive_impl::Kmer;
use std::collections::BTreeMap;

use crate::ofiles::Ofilesr;
use crate::ifiles::Ifilesr;
use crate::geneids::GeneIds;

use std::path::Path;
use std::io::Write;
use std::fs;

use std::io::BufRead;
use std::io::Read;


/// GeneIds harbors the antibody tags, the mRNA tags and whatever tags more you search the R2 for
/// but that can be different, too. 
/// Hence we store here 
/// mapper      : the search vector of GeneIDs
/// names       : a hashset for the gene names
/// names4sparse: and internal BtreeMap for the data export
/// last_count  : internal value for data export

/// and a lot of private ease of live id to name or vice versa BTreeMaps

#[derive(Debug,PartialEq)]
pub struct FastMapper<'a>{
    pub mapper: Vec<&'a mut GeneIds>, // the position is the sequence!
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names_store: Vec<String>, //sore gene names as vector
    pub max_id: usize,// hope I get the ids right this way...
    pub last_count: usize,
    pub last_kmers: Vec<String>,
    pub with_data: usize,
}

// here the functions

/* Seams I need to implement a Clone function for my GeneIds class. What a pain!

stefan.lang@med.lu.se
Hi and good Morning. How do I create a to_owned functionality for my own Rust class?

Good morning! To create a to_owned functionality for your own Rust class, you can implement the Clone trait for your class. The Clone trait defines a clone method that returns a new instance of the object with the same data as the original instance.

Here's an example implementation:

rust

#[derive(Debug, Clone)]
struct MyStruct {
    // fields of your struct
}

impl MyStruct {
    fn to_owned(&self) -> Self {
        self.clone()
    }
}

In this example, we define a MyStruct struct and implement the Clone trait for it. The to_owned method simply calls the clone method to create a new instance of the struct with the same data as the original instance. We can then use to_owned to create a new owned instance of MyStruct from a reference to an existing instance.

Note that we also derive the Debug trait for MyStruct to allow for debugging output using println!("{:?}", my_struct).

Here's an example usage:

rust

fn main() {
    let my_struct = MyStruct { /* initialize fields */ };
    let owned_struct = my_struct.to_owned();
    // `owned_struct` is a new owned instance of `MyStruct` with the same data as `my_struct`
}

I hope this helps! Let me know if you have any more questions.
stefan.lang@med.lu.se
:-D And the clone function - just re-run the new function?

Yes, to implement the Clone trait, you need to define a clone method that returns a new instance of the struct with the same data as the original instance. You can do this by simply copying the values of the fields from the original instance to the new instance.

Here's an example implementation:

rust

#[derive(Debug, Clone)]
struct MyStruct {
    // fields of your struct
}

impl Clone for MyStruct {
    fn clone(&self) -> Self {
        // Create a new instance of the struct with the same data as `self`
        MyStruct {
            // Copy the values of the fields from `self` to the new instance
            // For example:
            // field1: self.field1.clone(),
            // field2: self.field2.clone(),
            // ...
        }
    }
}

In this example, we implement the Clone trait for MyStruct by defining a clone method that creates a new instance of the struct with the same data as self.

To copy the values of the fields from self to the new instance, you can call the clone method on each field that implements the Clone trait. If a field does not implement the Clone trait, you may need to implement it yourself by defining a clone method for that field.

Once you've implemented the Clone trait for your struct, you can call the clone method on any instance of the struct to create a new instance with the same data.

I hope this answers your question! Let me know if you have any more questions.
*/

impl <'a> FastMapper<'a>{
    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new( )-> Self {
        let mapper = Vec::with_capacity(u16::MAX as usize);
        let names = BTreeMap::<std::string::String, usize>::new();
        let names4sparse = BTreeMap::<std::string::String, usize>::new();
        let names_store = Vec::with_capacity(1000);
        let max_id = 0;
        let last_count = 0;
        let last_kmers = Vec::with_capacity(60);
        let with_data = 0;
        Self {
            mapper,
            names,
            names4sparse,
            names_store,
            max_id,
            last_count,
            last_kmers,
            with_data,
        }
    }

    pub fn add(&mut self, seq: &[u8], name: std::string::String ){
        
        let mut total = 0;
        let mut add = 0;

        for i in 0..3{
            // add three kmers at position 0+4 4+4 and 8+4
            let kmer = needletail::kmer::Kmers::new(&seq[(1*i)..(8+8*i)], 8 as u8 ).next().unwrap();
            for value in Kmer::from(kmer).k.to_le_bytes(){
                let index = value as usize;
                let mut gene_ids = self.mapper.get_mut(index)
                    .unwrap_or_else(|| {
                        self.mapper.insert(index, &mut Box::new(&mut GeneIds::new(32) ) );
                        self.mapper.get_mut( index).unwrap()
                        });
                gene_ids.add(&seq[(8 + 8 * i)..(36 + 8 * i)], name.to_string() );
            }
            total +=1;
        }
        if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.max_id );
            self.names_store.push( name.clone() );
            self.max_id += 1;
        }

        println!( "{} kmers for gene {}", total, name );

    }



    pub fn print( &self ){
        println!("I have {} kmers and {} genes", self.with_data, self.names.len() );
    }

    // pub fn get_id( &mut self, name: String ) -> usize{
    //     let id = match self.names.get( &name ) {
    //         Some( id ) => id,
    //         None => panic!("Gene {name} not defined in the GeneID object"),
    //     };
    //     *id
    // }

    // pub fn get_name( &self, id:usize) -> String{
    //     let name = match self.names_store.get( &id ){
    //         Some( na ) => na,
    //         None => panic!("GeneID {id} not defined in the GeneID object"),
    //     };
    //     name.to_string()
    // }

    // pub fn get(&mut self, seq: &[u8] ) -> Option< usize >{
        
    //     // let min_value = 2;
    //     // let min_z = 1;
    //     // let mut max_value = 0;
    //     // let mut ret:u32 = 0;
    //     if self.unchecked{
    //         return self.get_unchecked( seq );
    //     }
    //     let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
    //     let mut kmer_vec = Vec::<u64>::with_capacity(60);

    //     fill_kmer_vec(kmers, &mut kmer_vec);

    //     //let report = self.get_id("2810417H13Rik".to_string());

    //     let mut ret:Option<usize> = None;

    //     if kmer_vec.is_empty() {
    //         //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
    //         return ret
    //     }  
    //     let mut sums = vec![0 ;self.names.len()];
    //     let mut max = 0;
    //     self.last_kmers.clear();

    //     for km in kmer_vec{
    //         //println!( "searching for kmer {}", km);
    //         if let Some(c1) = self.kmers.get(&km) {
    //             sums[*c1] += 1;
    //             if max < sums[*c1]{
    //                 max =  sums[*c1];
    //             }
    //             self.last_kmers.push(  self.kmer_store.get( &km).unwrap().to_string() );
    //         }

    //     }
    //     self.last_count = max;
    //     if max >2 {
    //         for (i, sum) in sums.iter().enumerate() {
    //             if *sum == max{
    //                 //println!("Now the ret hould have the value {} resp {:?}", i, Some(i));
    //                 if ret.is_some() {
    //                     //eprintln!("I have two genes matching with max value of {}: {:?} and {}", max, ret, i);
    //                     ret =None;
    //                     return ret
    //                 }
    //                 ret = Some(i);
    //             }
    //         }
    //     }
    //     //println!("return geneid {:?}", ret);
    //     ret
    // }

    // pub fn get_unchecked(&mut self, seq: &[u8] ) -> Option< usize >{
        
    //     // let min_value = 2;
    //     // let min_z = 1;
    //     // let mut max_value = 0;
    //     // let mut ret:u32 = 0;
    //     let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
    //     let mut kmer_vec = Vec::<u64>::with_capacity(60);

    //     fill_kmer_vec(kmers, &mut kmer_vec);

    //     let mut ret:Option<usize> = None;

    //     if kmer_vec.is_empty() {
    //         //eprintln!( "bad sequence: {:?}", std::str::from_utf8( seq ) );
    //         return ret
    //     }  
    //     let mut sums = vec![0 ;self.names.len()];
    //     let mut max = 0;
    //     self.last_kmers.clear();

    //     for km in kmer_vec{
    //         //println!( "searching for kmer {}", km);
    //         if let Some(c1) = self.kmers.get(&km) {
    //             sums[*c1] += 1;
    //             if max < sums[*c1]{
    //                 max =  sums[*c1];
    //             }
    //             self.last_kmers.push(  self.kmer_store.get( &km).unwrap().to_string() );
    //         }
    //     }
    //     self.last_count = max;
    //     if max >3 {
    //         for (i, sum) in sums.iter().enumerate(){
    //             if *sum == max{
    //                 if ret.is_some() {
    //                     //eprintln!("I have two genes matching with max value of {}: {:?} and {}", max, ret, i);
    //                     ret =None;
    //                     return ret
    //                 }
    //                 //println!("max = {} -> now the ret hould have the value {:?}", max, Some(i));
    //                 ret = Some(i);
    //                 //break;
    //             }
    //         }
    //     }
    //     // else {
    //     //     eprintln!("The max was {} => no (good) gene found", max);
    //     // }
    //     //println!("return geneid {:?}", ret);
    //     ret
    // }
    // // pub fn to_ids( &self,  ret:&mut Vec<Info> )  {
    // //     ret.clear();
    // //     for (_i, obj) in &self.kmers {
    // //         ret.push(*obj );
    // //     }
    // // }

    // pub fn to_header( &self ) -> std::string::String {
    //     let mut ret= Vec::<std::string::String>::with_capacity( self.names.len() +2 );
    //     //println!( "I get try to push into a {} sized vector", self.names.len());
    //     for obj in self.names.keys() {
    //         //println!( "Pushing {} -> {}", obj, *id-1);
    //         ret.push(  obj.to_string() ) ;
    //     }
    //     ret.push("Most likely name".to_string());
    //     ret.push("Faction total".to_string());
    //     "CellID\t".to_owned()+&ret.join("\t")
    // }

    // pub fn to_header_n( &self, names: &Vec<String> ) -> std::string::String {
    //     let mut ret= Vec::<std::string::String>::with_capacity( names.len() +2 );
    //     //println!( "I get try to push into a {} sized vector", self.names.len());
    //     for name in names {
    //         //println!( "Pushing {} -> {}", obj, *id-1);
    //         ret.push( name.to_string() ) ;
    //     }
    //     ret.push("AsignedSampleName".to_string());
    //     ret.push("FractionTotal".to_string());
    //     "CellID\t".to_owned()+&ret.join("\t")
    // }

    // pub fn write_index( &mut self, path: String ) -> Result< (), &str>{
    //     let rs = Path::new( &path ).exists();

    //     if ! rs {
    //         match fs::create_dir ( path.clone() ){
    //             Ok(_file) => (),
    //             Err(err) => {
    //                  eprintln!("Error?: {err:#?}");
    //              }
    //         };
    //     }

    //     // remove the old files if they exist:
    //     let fpath = Path::new(&path ) ;
    //     if fs::remove_file(fpath.join("index.1.Index.gz") ).is_ok(){};
    //     if fs::remove_file(fpath.join("index.1.gene.txt.gz") ).is_ok(){};

    //     let mut ofile = Ofilesr::new( 1, "index", "Index", "gene.txt",  &path );

    //     let seq_len:u64 = self.kmer_size as u64;


    //     match ofile.buff1.write( &seq_len.to_le_bytes() ){
    //         Ok(_) => (),
    //         Err(_err) => return Err::<(), &str>("seq_len could not be written"),
    //     };

    //     for ( kmer, gene_id) in &self.kmers {
    //         match ofile.buff1.write( &kmer.to_le_bytes() ){
    //             Ok(_) => println!("{} -> {:?} -> {:?}",kmer, &kmer.to_le_bytes(),  &u64::from_le_bytes(kmer.to_le_bytes())  ) ,
    //             Err(_err) => return Err::<(), &str>("kmer could not be written"),
    //         };

    //         match self.names_store.get( &gene_id ){
    //             Some(name) => {
    //                 write!(ofile.buff2,"{name}\n" ).unwrap();
    //             },
    //             None => panic!("Gene id {gene_id} not found in this object!"),
    //         };
    //     }


    //     ofile.close();
    //     Ok(())
    // }

    // pub fn load_index( &mut self, path: String ) -> Result< (), &str>{

    //     let f1 = Path::new( &path ).join("index.1.Index");
    //     let f2 = Path::new( &path ).join("index.1.gene.txt");
    //     if ! f1.exists() {
    //         panic!("Index file {} does not exist", f1.to_str().unwrap() );
    //     }
    //     if ! f2.exists() {
    //         panic!("kmer names file {} does not exist", f2.to_str().unwrap() );
    //     }

    //     let mut ifile =Ifilesr::new( 1, "index", "Index", "gene.txt",  &path  );
    //     let mut id = 0;

    //     let mut buff = [0_u8 ;8 ];

    //     ifile.buff1.read_exact(&mut buff).unwrap();
    //     //println! ("I got the buff {buff:?}");
    //     self.kmer_size = u64::from_le_bytes( buff ) as usize;

    //     for name in ifile.buff2.lines() {
    //         //println!("I read this name: {name:?}");
    //         //let mut buff = [0_u8 ;8 ];
    //         ifile.buff1.read_exact(&mut buff).unwrap();

    //         let km = &u64::from_le_bytes( buff );
            
    //         //println!("I try to insert km {km} and gene name '{name}' ({len:?}) ");
    //         id += 1;
    //         let mut data:String = "".to_string() ; // needs to be a constant...
    //         match name{
    //             Ok(na) => {
    //                 if let std::collections::btree_map::Entry::Vacant(e) = self.kmers.entry(*km) {
    //                     //let info = Info::new(km, name.clone() );
    //                     if ! self.names.contains_key( &na ){
    //                         self.names.insert( na.clone(), self.max_id );
    //                         self.names_store.insert( self.max_id, na.clone() );
    //                         self.max_id += 1;
    //                     }
    //                     let gene_id = self.names.get( &na ).unwrap();
    //                     //self.kmer_store.insert( *km , "not recovered".to_string() );
    //                     //println!("These are the bytes I need to work on: {:?}",&km.to_le_bytes());
    //                     Self::u64_to_str(self.kmer_size, km, &mut data );
    //                     println!( "({})Can I get that {km} to string? {:?}", self.seq_len, data );
                        
    //                     self.kmer_store.insert( *km , data.to_string() );
    //                     e.insert( *gene_id );
    //                 } 
    //                 else {
    //                     self.bad_entries.insert( *km );
    //                     self.kmers.remove( km );
    //                 }
    //             },
    //             Err(err) => panic!("This is not working: {err:?}"),
    //         };

    //         //println!(" {:?} -> {:?}", buff,  &u64::from_le_bytes(buff)  ) ;
    //     }
    //     if id == 0 {
    //         panic!("No data read!");
    //     }

    //     //self.print();
    //     Ok(())
    // }

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
    use crate::geneids::FastMapper;
    use std::path::Path;

     #[test]
    fn check_geneids() {
        let mut mapper = FastMapper::new();

        let mut geneid = 0;
        
        mapper.add( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC", "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", "Gene2".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 6 );
        mapper.print();

    }

    #[test]
    fn test_u64_to_str(){

        let num:u64 = 15561;
        let mut data:String = "".to_string();
        GeneIds::u64_to_str( 7, &num, &mut data );
        assert_eq!( data, "CGATATT".to_string() )
    } 
}

/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

//use needletail::kmer::Kmers;
use kmers::naive_impl::Kmer;
use std::collections::BTreeMap;
//use needletail::bitkmer::BitKmer;

use crate::ofiles::Ofilesr;
use crate::ifiles::Ifilesr;
use crate::mapper_entries::MapperEntry;

use std::path::Path;
use std::io::Write;
use std::fs::{self, DirBuilder};

use std::io::BufRead;
use std::io::Read;

//use crate::geneids::GeneIds;


/// GeneIds harbors the antibody tags, the mRNA tags and whatever tags more you search the R2 for
/// but that can be different, too. 
/// Hence we store here 
/// mapper      : the search vector of GeneIDs
/// names       : a hashset for the gene names
/// names4sparse: and internal BtreeMap for the data export
/// last_count  : internal value for data export

/// and a lot of private ease of live id to name or vice versa BTreeMaps

#[derive(Debug,PartialEq)]
pub struct FastMapper{
    pub kmer_len:usize, // how long should the second entries be (up to 32 bp + 8bp initial map) 
    spacer:usize,
    pub mapper: Vec<Box<MapperEntry>>, // the position is the sequence!
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names_store: Vec<String>, //sore gene names as vector
    pub max_id: usize,// hope I get the ids right this way...
    pub last_count: usize,
    pub last_kmers: Vec<String>,
    pub with_data: usize,
}

// here the functions

impl  FastMapper{
    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new( kmer_len:usize )-> Self {
        println!("I would add {} new entries here:", u16::MAX);
        let mut mapper: Vec<Box< MapperEntry>> = Vec::with_capacity(u16::MAX as usize);
        for _i in 0..u16::MAX{
            let b = Box::new( MapperEntry::new( ) );
            mapper.push( b );
        }
        let names = BTreeMap::<std::string::String, usize>::new();
        let names4sparse = BTreeMap::<std::string::String, usize>::new();
        let names_store = Vec::with_capacity(1000);
        let max_id = 0;
        let last_count = 0;
        let last_kmers = Vec::with_capacity(60);
        let with_data = 0;
        let spacer = 3;
        Self {
            kmer_len,
            spacer,
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



    pub fn add(&mut self, seq: &Vec<u8>, name: std::string::String ){
        
        let mut total = 0;

        if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.max_id );
            self.names_store.push( name.clone() );
            self.max_id += 1;
        }

        let gene_id = self.get_id( name.to_string() );


        'this: for i in 0..10{
            // add three kmers at position 0+4 4+4 and 8+4
            if seq.len() < self.kmer_len + 3 + 8 * i{
                break;
            }
            let mut longer: u64 = 0;
            for kmer_l in needletail::kmer::Kmers::new(&seq[(8 + self.spacer * i)..( self.kmer_len +8 + self.spacer * i)], self.kmer_len as u8){
                //println!("And I try to get a longer from kmer_l {kmer_l:?}:");
                for nuc in kmer_l {
                    if *nuc ==b'N'{
                        continue 'this;
                    }
                }
                longer = Kmer::from(kmer_l).into_u64();
                break;
            }


            //let kmer = needletail::kmer::Kmers::new(&seq[(0 + self.spacer *i)..(8+self.spacer*i)], 8 as u8 ).next().unwrap();

            let mut index : usize = 0;
            //println!("Is this a [u8] ({}, {})?: {:?}", (0 + 8 * i), ( 8 + 8 * i ), &seq[(8 + 8 * i)..(8 + 8 * i)]);
            for kmer in needletail::kmer::Kmers::new(&seq[(0+ self.spacer * i)..(8 +self.spacer * i)], 8 as u8){
                //println!("And I try to get a index from kmer {kmer:?}:");
                for nuc in kmer {
                    if *nuc ==b'N'{
                        continue 'this;
                    }
                }
                index = Kmer::from(kmer).into_u64() as usize;
                break;
            }

            //println!("This is what I would add here: {index} and {longer}");


            match self.mapper.get_mut(index){
                Some( genebox ) => {
                    if genebox.map.len() == 0{
                        self.with_data +=1;
                    }
                    genebox.add( longer , gene_id );
                    //println!("Cool I got a gene box! {genebox:?}");
                },
                None => {
                    let mut b = Box::new( MapperEntry::new( ) );
                    b.add( longer , gene_id  );
                    //println!("Cool I got a new gene box!? {b:?}");
                    self.mapper.insert(index,  b ); 
                }
            }
            total +=1;
        }
        
        if total == 0 {
            eprintln!("gene {} does not get an entry in the fast_mapper object! - too short!!", &name.as_str() );
        }
        //println!( "{} kmers for gene {} ({})", total, &name.to_string(), gene_id );

    }

    pub fn print( &self ){
        println!("I have {} kmers and {} genes", self.with_data, self.names.len() );
    }

    pub fn get_id( &mut self, name: String ) -> usize{
        let id = match self.names.get( &name ) {
            Some( id ) => id,
            None => panic!("Gene {name} not defined in the GeneID object"),
        };
        *id
    }

    pub fn get_name( &self, id:usize) -> String{
        let name = match self.names_store.get( id ){
            Some( na ) => na,
            None => panic!("GeneID {id} not defined in the GeneID object"),
        };
        name.to_string()
    }

    pub fn get(&mut self, seq: &[u8] ) -> Option< usize >{
        
        let mut id = 0;

        for kmer in needletail::kmer::Kmers::new(&seq, 8 as u8){
            //println!("And I try to get a index from kmer {kmer:?} at pos {id}:");
            let mut ok = true; 
            if seq.len() < 40 + 3 * id{
                break
            }
            for n in &seq[(id)..(self.kmer_len +8 + id)]{
                if  n == &b'N' {
                    ok = false;
                }
            }
            if ! ok{
                continue;
            }
            
            let index = Kmer::from(kmer).into_u64() as usize;
            let mapper = match self.mapper.get_mut(index){
                Some( map) => map,
                None => return None
            };
            //println!("got the index {index} and the mapper {mapper:?} searching for");
            //let mapper:MapperEntry = &*box_;
            if mapper.has_data() {
                let mut longer: u64;
                for kmer_l in needletail::kmer::Kmers::new(&seq[(8 + id)..(self.kmer_len +8 + id)], self.kmer_len as u8){
                    //println!("And I try to get a longer from kmer_l {kmer_l:?}:");
                    longer = Kmer::from(kmer_l).into_u64();
                    //println!("got the index {index} and the mapper {mapper:?} searching for longer {longer}");
                    match mapper.get(&longer){
                        Some( gene_id ) => {
                            //println!("Got one: {gene_id}");
                            return Some(gene_id);
                        },
                        None => (),
                    }
                }
            }
            id +=1;
        }
        None
    }

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

    pub fn to_header_n( &self, names: &Vec<String> ) -> std::string::String {
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

    pub fn write_index( &mut self, path: String ) -> Result< (), &str>{
        let rs = Path::new( &path ).exists();

        if ! rs {
            let mut dir_builder = DirBuilder::new();
            dir_builder.recursive(true);

            match dir_builder.create ( path.clone() ){
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

        
        let mut count:usize;
        let mut idx:usize;

        let seq_len:u64 = self.kmer_len as u64;
        match ofile.buff1.write( &seq_len.to_le_bytes() ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("seq_len could not be written"),
        };

        let with_data:u64 = self.with_data as u64;
        match ofile.buff1.write( &with_data.to_le_bytes() ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("with_data could not be written"),
        };


        for i in 0..u16::MAX{
            // write the 8bp kmer as int
            idx = i as usize;
            if self.mapper[idx].has_data(){
                match ofile.buff1.write( &i.to_le_bytes() ){
                    Ok(_) => println!("i: {} -> {:?}",i, &i.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("i could not be written"),
                };
                //write the amount of downstream entries
                count = self.mapper[idx].map.len();
                match ofile.buff1.write( &count.to_le_bytes() ){
                    Ok(_) => println!("count: {} -> {:?}",count, &count.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("count could not be written"),
                };
                for (key, value) in self.mapper[idx].map.iter(){
                    match ofile.buff1.write( &key.to_le_bytes() ){
                        Ok(_) => println!("key: {} -> {:?}",key, &key.to_le_bytes()  ) ,
                        Err(_err) => return Err::<(), &str>("key could not be written"),
                    };
                    match ofile.buff1.write( &value.to_le_bytes() ){
                        Ok(_) => println!("value: {} -> {:?}",value, &value.to_le_bytes()  ) ,
                        Err(_err) => return Err::<(), &str>("value could not be written"),
                    };
                }
            }
        }

        let string_to_write = self.names_store.join("\n");
        match ofile.buff2.write( string_to_write.as_bytes() ){
            Ok(_) => println!("I wrote this: {string_to_write}",),
            Err(_err) => return Err::<(), &str>("gene names could not be written"),
        };

        ofile.close();
        Ok(())
    }

    pub fn load_index( &mut self, path: String ) -> Result< (), &str>{

        let f1 = Path::new( &path ).join("index.1.Index");
        let f2 = Path::new( &path ).join("index.1.gene.txt");
        if ! f1.exists() {
            panic!("Index file {} does not exist", f1.to_str().unwrap() );
        }
        if ! f2.exists() {
            panic!("kmer names file {} does not exist", f2.to_str().unwrap() );
        }

        let mut ifile =Ifilesr::new( 1, "index", "Index", "gene.txt",  &path  );

        let mut buff_u64 = [0_u8 ;8 ];
        let mut buff_u16 = [0_u8 ;2 ];

        let mut kmer:u64;
        let mut gene_id:usize;

        // first 16 bytes are the kmer_len and the with_data usize values.
        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        //println! ("I got the buff {buff:?}");
        self.kmer_len = u64::from_le_bytes( buff_u64 ) as usize;

        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        self.with_data = u64::from_le_bytes( buff_u64 ) as usize;

        for _i in 0..self.with_data {
            // read in the single mapper elements 
            // first 2 - kmer position
            ifile.buff1.read_exact(&mut buff_u16).unwrap();
            let idx = u16::from_le_bytes ( buff_u16) as usize;
            if self.mapper[idx].has_data(){
                panic!("Loading an mapper_entry that has alreadyx been filled!");
            }
            // next 8 amount of [u64, gene id] following
            ifile.buff1.read_exact(&mut buff_u64).unwrap();
            let i = u64::from_le_bytes( buff_u64 ) as usize;
            for _i in 0..i{
                // next 8 u64 kmer
                ifile.buff1.read_exact(&mut buff_u64).unwrap();
                kmer = u64::from_le_bytes( buff_u64 );
                // next 8 u64 gene id
                ifile.buff1.read_exact(&mut buff_u64).unwrap();
                gene_id = usize::from_le_bytes( buff_u64 );
                // save in object
                self.mapper[idx].map.insert( kmer, gene_id );
            }
        }


        for name in ifile.buff2.lines() {
            //println!("I read this name: {name:?}");
            //let mut buff = [0_u8 ;8 ];
            match name {
                Ok(n) => {
                    if ! self.names.contains_key( &n ){
                        self.names.insert( n.clone(), self.max_id );
                        self.names_store.push( n.clone() );
                        self.max_id += 1;
                    }
                },
                Err(e) => panic!("I could not read from the genes table {e:?}"),
            };
        }
            
        if self.with_data == 0 {
            panic!("No data read!");
        }

        //self.print();
        Ok(())
    }

}



#[cfg(test)]
mod tests {

    use crate::fast_mapper::FastMapper;
    use std::path::Path;
    use crate::fast_mapper::Kmer;

    #[test]
    fn check_geneids() {
        let mut mapper = FastMapper::new( 16 );

        let mut geneid = 0;
        
        mapper.add( b"ATCCCATCCTTCATTGTTCGCCTGGA", "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", "Gene2".to_string() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 9 );

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC" ), Some(0));
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC" ), Some(1));

        mapper.add( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", "Gene3".to_string() );

        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC" ), Some(2));

        let mut gnames = Vec::<String>::with_capacity(3);
        gnames.push( "Gene1".to_string() );
        gnames.push( "Gene2".to_string() );
        gnames.push( "Gene3".to_string() );
        assert_eq!( mapper.names_store, gnames );
        mapper.print();

    }

    #[test]
    fn check_write_index() {
        let mut mapper = FastMapper::new( 16 );

        let mut geneid = 0;
        
        mapper.add( b"ATCCCATCCTTCATTGTTCGCCTGGA", "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", "Gene2".to_string() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 9 );

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC" ), Some(0));
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC" ), Some(1));

        let opath = "testData/output_index_test";
        mapper.write_index( opath.to_string() );
        mapper.print();

        let mut mapper2 = FastMapper::new( 16 );
        mapper2.load_index( opath.to_string() );

        assert_eq!(  mapper.with_data, mapper2.with_data );

        assert_eq!(  mapper.kmer_len, mapper2.kmer_len );

        assert_eq!(  mapper.kmer_len, mapper2.kmer_len );

        let key = b"ATCCCATC";
        let kmer = needletail::kmer::Kmers::new(key, 8 as u8).next();
        let idx = Kmer::from(kmer.unwrap()).into_u64() as usize;

        assert_eq!(  mapper.names_store, mapper2.names_store );

    }

}

/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

//use needletail::kmer::Kmers;
use kmers::naive_impl::Kmer;
use std::collections::BTreeMap;
//use needletail::bitkmer::BitKmer;
use crate::int_to_str::IntToStr;
use std::collections::HashMap;

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
    pub mapper: Vec<MapperEntry>, // the position is the sequence!
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names_store: Vec<String>, //store gene names as vector
    pub max_id: usize,// hope I get the ids right this way...
    pub last_count: usize,
    pub last_kmers: Vec<String>,
    pub with_data: usize,
    pub version: usize, //the version of this
    pub tool: IntToStr,

}

// here the functions

impl  FastMapper{
    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new( kmer_len:usize )-> Self {
        //println!("I would add {} new entries here:", u16::MAX);
        let mut mapper: Vec<MapperEntry> = Vec::with_capacity(u16::MAX as usize);
        for _i in 0..u16::MAX{
            let b = MapperEntry::new( );
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
        let tool = IntToStr::new(b"".to_vec(), kmer_len );
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
            version: 2,
            tool,
        }
    }

    /// returns true if the sequence passes, fasle if there is a problem like N nucleotides
    /// or too low sequences variability in the first 8bp
    fn seq_ok(&self, seq:&Vec<u8>, i:usize, checker:&mut BTreeMap::<u8, usize>) -> Option<bool> {

        let start = (i+1)*8;

        if seq.len() < (i+1)*8 +10 {
            return None
        }
        let mut to = (i+1)*8+self.kmer_len;
        if seq.len()< to{
            to = seq.len();
        }
        // check the initial 8 bp as they will create a key
        for nuc in &seq[i*8..(i+1)*8 ] {
            match checker.get_mut( nuc ){
                Some( count ) => *count += 1,
                None => {
                    if checker.insert( *nuc, 1).is_some(){};
                }
            };
            if *nuc ==b'N'{
                return Some(false);
            }
        }
        if checker.len() < 3{
            //println!( "kmer is too simple/not enough diff nucs" );
            return Some(false);
        }
        
        for ( key, value ) in checker.iter(){
            //println!( "sequence from {start} to {to} is too simple/too many nucs same" );

            if *value as f32 / (to-start) as f32 > 0.6 {
                //println!( "sequence from {start} to {to} is too simple/too many nucs same" );
                return Some(false);
            } 
        }
        
        // and now check the rest, too
        checker.clear();

        for nuc in &seq[(i+1)*8..to ] {
            match checker.get_mut( nuc ){
                Some( count ) => *count += 1,
                None => {
                    if checker.insert( *nuc, 1).is_some(){};
                }
            };
            if *nuc ==b'N'{
                return Some(false);
            }
        }
        if checker.len() < 3{
            //println!( "kmer is too simple/not enough diff nucs" );
            return Some(false);
        }
        
        for ( key, value ) in checker.iter(){
            //println!( "sequence from {start} to {to} is too simple/too many nucs same" );

            if *value as f32 / (to-start) as f32 > 0.6 {
                //println!( "sequence from {start} to {to} is too simple/too many nucs same" );
                return Some(false);
            } 
        }




        return Some(true)
    }



    pub fn add(&mut self, seq: &Vec<u8>, name: std::string::String ){
        
        let mut total = 0;

        println!("fast_mapper adds this genes: {}", name);

        if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.max_id );
            self.names_store.push( name.clone() );
            self.max_id += 1;
        }
        
        //println!("Adding gene {name} to the index:");

        let gene_id = self.get_id( name.to_string() );

        //let mut long:u64;
        //let mut short:u16;
        let mut checker = BTreeMap::<u8, usize>::new();
        self.tool.from_vec_u8( seq.to_vec() );
        let mut tmp = "".to_string();
        self.tool.to_string( seq.len(), &mut tmp );
        self.tool.print();
        //println!("Adding gene to the index: \n\t>{name}\n\t{tmp}");
        let mut i = 0;

        self.tool.shifted = 3; // I will only add one set of kmers!

        let mut item = self.tool.next();
        while let Some(entries) = item{

            // checker.clear();
            // match self.seq_ok( &seq.to_vec(), i, &mut checker){
            //     Some(false) => {
            //         i+=1;
            //         continue 'main
            //     },
            //     Some(true) => (),
            //     None => break 'main,

            // }
            println!( "short {:b} long: {:b}",entries.0, entries.1);

            // comment out if not debugging
            let mut tmp = "".to_string();
            self.tool.u8_array_to_str( 8, entries.0.clone().to_le_bytes().to_vec(), &mut tmp );
            print!("\t\t\tindex\t{} {} ", &tmp, &entries.0 );
            print!("[ ");
            for en in entries.0.to_le_bytes().to_vec(){
                print!(", {en:b}");
            }
            println!(" ]");
            //

            // comment out if not debugging
            tmp.clear();
            self.tool.u8_array_to_str( 32, entries.1.to_le_bytes().to_vec(), &mut tmp );
            println!("\t\t\tlong\t{tmp}, {}",entries.1);
            //

            if ! self.mapper[entries.0 as usize].has_data() {
                // will add in the next step so
                self.with_data +=1;
            } 
            self.mapper[entries.0 as usize].add( entries.1, gene_id );
            i+=1;
            item = self.tool.next();
        }
        
        //println!("I have now {i} mappers for the gene {name}");
        if i == 0 {
            panic!("gene {} does not get an entry in the fast_mapper object! - too short!!", &name.as_str() );
        }

        self.with_data += total;
        //println!( "{} kmers for gene {} ({})", total, &name.to_string(), gene_id );

    }
    /// [entries with values, total second level entries, total single gene first/second pairs, total >1 first/second level pairs]
    pub fn info( &self ) -> [usize;4]{
        let mut ret: [usize;4] = [0,0,0,0];
        let mut second: [usize;3];
        for entry in &self.mapper {
            if entry.has_data() {
                ret[0] +=1;
                second = entry.info();
                for id in 0..3{
                    ret[id+1] += second[id];
                }
            }
        }
        ret
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

    pub fn get(&mut self, seq: &[u8] ) -> Option<usize>{
        
        let mut id = 0;
        let mut genes:HashMap::<usize, usize>= HashMap::new();

        let mut checker = BTreeMap::<u8, usize>::new();
        let mut total = 0;
        let mut i =0;

        let mut long:u64;
        let mut short:u16;
        let mut stop: bool;
        let mut tmp = "".to_string();
        //self.tool.shifted = 3; // I will only add one set of kmers!

        let mut item = self.tool.next();
        'main :while let Some(( short, long )) = item{

        //'main: while i < 30{
            //println!("And I try to get a index from kmer {kmer:?} at pos {id}:");

            // checker.clear();
            // match self.seq_ok( &seq.to_vec(), i, &mut checker){
            //     Some(false) => {
            //         i+=1;
            //         continue 'main
            //     },
            //     Some(true) => (),
            //     None => break 'main,

            // }

            tmp.clear();
            self.tool.u8_array_to_str( 8, short.to_le_bytes().to_vec(), &mut tmp );
            println!("\t\t\tindex\t{tmp} {short}" );
            tmp.clear();
            self.tool.u8_array_to_str( 32, long.to_le_bytes().to_vec(), &mut tmp );
            println!("\t\t\tlong\t{tmp} {long}" );


            //println!("got the index {index} and the mapper {mapper:?} searching for");
            //let mapper:MapperEntry = &*box_;
            stop = false;
            if self.mapper[short as usize].has_data() {
                match self.mapper[short as usize].get(&long){
                    Some( gene_id ) => {
                        //println!("Got one: {gene_id:?}");
                        //return ( gene_id );
                        for gid in gene_id{
                            match genes.get_mut(&gid) {
                                Some(mut gene_count) => {
                                    *gene_count +=1;
                                    if *gene_count == 4 {
                                        stop = true;
                                    }
                                },
                                None => {
                                    genes.insert( gid, 1);
                                },
                            };
                        }
                    },
                    None => (),
                }
            }
            i +=1;
            if stop{
                break 'main;
            }
        }
        let mut max = 0;
        id = 0;
        if genes.len() > 0{
            for (gene_id, count) in &genes{
                if count > &max{
                    max = *count;
                    id = *gene_id;
                    i = 0;
                }else if count == &max {
                    i+=1;
                }
            }
            if i == 0{
                // we have exactly one gene with the max count
                return Some(id);
            }else {
                // we have more than one gene as a max count - that is inacceptable
                // for the debug I need to know which genes I get here!
                let mut gene_ids= Vec::<usize>::with_capacity(i);
                for (gene_id, count) in &genes{
                    if count == &max {
                        gene_ids.push(*gene_id);
                    }
                }
                if gene_ids.len() == 2 && ( gene_ids[0] +1 == gene_ids[1] || gene_ids[0] == gene_ids[1] +1 ) {
                    // this is the gene and it's _int transcript not having a clear difference.
                    // needs to be fixed in a later version of the index...
                    return Some(gene_ids[0] )
                }
                println!("For the sequence {}", std::str::from_utf8(&seq).expect("Invalid UTF-8") );
                println!("I got a multi mapper ({i}): {:?} -> returning None", gene_ids );
                return None
            }
            
        }
        // we have not found any gene here!
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

        let version:u64 = self.version as u64;
        match ofile.buff1.write( &version.to_le_bytes() ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("version could not be written"),
        };

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
                    Ok(_) => (), 
                    //Ok(_) => println!("8bp mapper: {:b} binary -> {:?} bytes",i, &i.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("i could not be written"),
                };
                //write the amount of downstream entries
                count = self.mapper[idx].map.len();
                match ofile.buff1.write( &count.to_le_bytes() ){
                    Ok(_) => (), 
                    //Ok(_) => println!("amount of 32 bp second kmers to the first 8bp: {} -> {:?} bytes",count, &count.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("count could not be written"),
                };
                for tuple  in self.mapper[idx].map.iter(){
                    match ofile.buff1.write( &tuple.0.to_le_bytes() ){
                        Ok(_) => (), 
                        //Ok(_) => println!("the u64 kmer: {} or {:b} binary -> {:?} bytes",tuple.0, tuple.0, &tuple.0.to_le_bytes()  ) ,
                        Err(_err) => return Err::<(), &str>("key could not be written"),
                    };
                    // now we need the NameEntry.len() to to_le_bytes()
                    match ofile.buff1.write( &tuple.1.data.len().to_le_bytes() ){
                        Ok(_) => (),
                        //Ok( len ) => println!("length of gene list attached to that 32bp kmer: {:?}",len   ) ,
                        Err(_err) => return Err::<(), &str>("value could not be written"),
                    };
                    for id in &tuple.1.data{
                        match ofile.buff1.write( &id.to_le_bytes() ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                    }
                    
                }
            }
        }

        let string_to_write = self.names_store.join("\n");
        match ofile.buff2.write( string_to_write.as_bytes() ){
            Ok(_) => (),//eprintln!("Genes indexed: {string_to_write}",),
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
        let mut len:usize;


        // the real first 8 bytes are the version
        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        //println! ("I got the buff {buff:?}");
        if ! self.version == u64::from_le_bytes( buff_u64 ) as usize{
            panic!("Sorry - this version of the index {} does not match the requirements {}.", u64::from_le_bytes( buff_u64 ), self.version);
        }

        self.kmer_len = u64::from_le_bytes( buff_u64 ) as usize;
        // next 8bytes are the kmer_len
        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        //println! ("I got the buff {buff:?}");
        self.kmer_len = u64::from_le_bytes( buff_u64 ) as usize;
        // and finally (preable) the with_data usize value
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
                len = usize::from_le_bytes( buff_u64 );
                for _id in 0..len{
                    ifile.buff1.read_exact(&mut buff_u64).unwrap();
                    gene_id = usize::from_le_bytes( buff_u64 );
                    self.mapper[idx].add( kmer, gene_id );
                }
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




/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

//use needletail::kmer::Kmers;
//use kmers::naive_impl::Kmer;
use std::collections::BTreeMap;
//use needletail::bitkmer::BitKmer;
use crate::int_to_str::IntToStr;
use std::collections::HashMap;
//use std::collections::HashSet;

use regex::Regex;

use crate::ofiles::Ofilesr;
use crate::ifiles::Ifilesr;
use crate::mapper_entries::MapperEntry;

use std::path::Path;
use std::io::Write;
use std::fs::{self, DirBuilder};

use std::io::BufRead;
use std::io::Read;
//use std::error::Error;
//use crate::geneids::GeneIds;
use crate::mapping_info::MappingInfo;

/// GeneIds harbors the antibody tags, the mRNA tags and whatever tags more you search the R2 for
/// but that can be different, too. 
/// Hence we store here 
/// mapper      : the search vector of GeneIDs
/// names       : a hashset for the gene names
/// names4sparse: and internal BtreeMap for the data export
/// last_count  : internal value for data export

/// and a lot of private ease of live id to name or vice versa BTreeMaps

static EMPTY_VEC: Vec<usize> = Vec::new();



#[derive(Debug,PartialEq)]
pub struct FastMapper{
    pub kmer_len:usize, // how long should the second entries be (up to 32 bp + 8bp initial map) 
    spacer:usize,
    pub mapper: Vec<MapperEntry>, // the position is the sequence!
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names_store: Vec<String>, // store gene names as vector
    pub max_id: usize,// hope I get the ids right this way...
    pub last_count: usize,
    pub with_data: usize,
    pub version: usize, //the version of this
    pub tool: IntToStr,
    mask:u64, //a mask to fill matching sequences to match the index's kmer_len
    pub pos: usize, // when creating the index - this is the amount of times we added a new matching sequence
    pub neg: usize, // when creating the index - this is the amount of times we failed to add a new matching sequence (as we already had this connection)
}

// here the functions

impl  FastMapper{
    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new( kmer_len:usize, allocate:usize )-> Self {
        //println!("I would add {} new entries here:", u16::MAX);
        let mut mapper: Vec<MapperEntry> = Vec::with_capacity(u16::MAX as usize);
        for _i in 0..u16::MAX{
            let b = MapperEntry::new( 10 );
            mapper.push( b );
        }
        let names = BTreeMap::<std::string::String, usize>::new();
        let names4sparse = BTreeMap::<std::string::String, usize>::new();
        let names_store = Vec::with_capacity( allocate );
        let max_id = 0;
        let last_count = 0;
        let with_data = 0;
        let spacer = 3;
        let tool = IntToStr::new(b"".to_vec(), kmer_len );
        //let mask:u64 = 0b11 << (2* 32- kmer_len);
        let mask: u64 = (1 << (2 * kmer_len)) - 1;
        let neg = 0;
        let pos = 0;
        Self {
            kmer_len,
            spacer,
            mapper,
            names,
            names4sparse,
            names_store,
            max_id,
            last_count,
            with_data,
            version: 5,
            tool,
            mask,
            pos,
            neg,
        }
    }

    /// For the multicore processing I need this object to be mergable
    /// All genes from the other object need to be incorporated and all primaryid - secondaryid compbos need to be collected.
    /// All other classes have to be copied, too. So we need to use our own incorporate_match_combo function.
    pub fn merge( &mut self, other: FastMapper ) {

        for (primary_id, mapper_object) in other.mapper.iter().enumerate(){
            for ( secondary_match, name_entry) in mapper_object.map.iter(){
                for idx in 0..name_entry.classes.len(){
                    // Now I need to copy the gene names from the other object to my own.
                    let gene_names = other.gene_names_for_ids( &name_entry.classes[idx] );
                    // And get my own ids for these names
                    //let gene_ids = self.ids_for_gene_names( &gene_names );
                    let _ = self.incorporate_match_combo( primary_id, *secondary_match as usize, gene_names[0].clone(), gene_names );
                }           
            }
        }

    }

    pub fn gene_names_for_ids( &self, ids:&Vec<usize> ) -> Vec<String> {
        let mut ret = Vec::<String>::with_capacity( ids.len() );
        for id in ids.iter(){
            ret.push( self.names_store[*id].clone())
        }
        ret
    }

    pub fn ids_for_gene_names( &mut self, ids:&Vec<String> ) -> Vec<usize> {
        let mut ret = Vec::<usize>::with_capacity( ids.len() );
        for name in ids.iter(){
            if ! self.names.contains_key( name ){
                self.names.insert( name.clone(), self.max_id );
                self.names_store.push( name.clone() );
                self.max_id += 1;
            }
            let id = self.get_id( name.to_string() ).clone();
            ret.push( id);
        }
        ret
    }

    pub fn incorporate_match_combo( &mut self, initial_match:usize, secondary_match:usize, name:String, ids: Vec<String> )  ->Result<(), &str>{

        if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.max_id );
            self.names_store.push( name.clone() );
            self.max_id += 1;
        }
        let gene_id = self.get_id( name.to_string() ).clone();
        let classes =  self.ids_for_gene_names( &ids );
        if ! self.mapper[initial_match].has_data() {
            // will add in the next step so
            self.with_data +=1;
        }
        match self.mapper[initial_match].add( secondary_match.try_into().unwrap(), ( gene_id, 0), classes ){
            Ok(_) => self.pos += 1,
            Err(_) => {
                //eprintln!("incorporate_match_combo - gene seq was already registered {}", name);
                self.neg +=1
            },
        };
        Ok( () )
    }

    pub fn fill_sequence_with_a( self, seq:&u64 ) -> u64{
        let filled_sed = seq | (!self.mask & (0b11 << (2 * self.kmer_len )));
        return filled_sed
    }


    pub fn add(&mut self, seq: &Vec<u8>, name: std::string::String, class_ids: Vec<String> ) -> usize{
        
        //println!("fast_mapper adds this genes: {} of length {}", name, seq.len() );

        if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.max_id );
            self.names_store.push( name.clone() );
            self.max_id += 1;
        }

        let classes =  self.ids_for_gene_names( &class_ids );
        let gene_id = self.get_id( name.to_string() );

        //let mut long:u64;
        //let mut short:u16;

        self.tool.from_vec_u8( seq.to_vec() );
        let mut tmp = "".to_string();
        self.tool.to_string( seq.len(), &mut tmp );
        //self.tool.print();
        //println!("Adding gene to the index: \n\t>{name}\n\t{tmp}");

        // self.tool.shift(); // I will only add one set of kmers!
        // self.tool.shift();
        // self.tool.shift();

        //let mut item = self.tool.next();
        let mut i =0;
        //let mut short = String::from("");
        //let mut long = String::from("");
        let mut n =20;

        while let Some(entries) = self.tool.next(){
            n -=1;
            if n == 0{
                break;
            }
            // index_read, longer_read, longer length
            //short.clear();
            //long.clear();
            //self.tool.u16_to_str( 8, &entries.0, &mut short);
            //self.tool.u64_to_str( entries.2, &entries.1, &mut long);
            //println!("I got some entry [short {short}, long {long}, length {}", entries.2 );

            if entries.2 < 5{
                // too short second entry
                break;
            }

            if ! self.mapper[entries.0 as usize].has_data() {
                // will add in the next step so
                self.with_data +=1;
            }
            // if entries.2  < self.kmer_len{
            //     println!("adding a smaller than expected sequence for gene {name}/{gene_id} ({})", entries.2);
            // }
            match self.mapper[entries.0 as usize].add( entries.1, (gene_id, 0), classes.clone() ){
                Ok(_) => self.pos += 1,
                Err(_) => self.neg +=1,
            };
            i+=1;
        }
        
        //println!("I have now {i} mappers for the gene {name}");
        if i == 0 {
            eprint!(".");
            //eprintln!("gene {} ({:?}) does not get an entry in the fast_mapper object! - too short? {} bp", &name.as_str(), classes, seq.len() );
        }
        /*else {
            println!("I added {i} mappper entries for gene {name}");
        }*/

        //self.with_data += i;
        //println!( "{} kmers for gene {} ({})", total, &name.to_string(), gene_id );
        gene_id
    }


    /// The index is assumed to have been build using a gene family object.
    /// Therefore we have a lot of almost sequences indexed here.
    /// It is to be assumed that the index contains many not informative 40 bp areas.
    /// This function tries to get the maximum out of the index by chaning the 40 index elements associated gene lists tp onyl contain the
    /// least diverse class of the matching elements (the new classes objects in this class).
    /// The class object is mean to have the same kind of entry in the same position of the vector - say
    /// [NameEntry.name, "gene_name", "family_name", "class_name"] - each position will be summed for 40bp mapping elements
    /// and only the contents of the least class will be used. The original gene name will also be used - that does not need to be part of the ids.
    pub fn make_index_te_ready( &mut self ) {

        //println!("Before make_index_te_ready:");
        //self.print();
        for mapper_entry in self.mapper.iter_mut() {
            if mapper_entry.has_data() {
                mapper_entry.collapse_classes();
            }
        }

        self.with_data = 0;
        for mapper_entry in &self.mapper{
            if mapper_entry.has_data() {
                self.with_data += 1;
            }
        }

        //println!("\nAfter:");
        //self.print();
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
        println!("I have {} kmers and {} genes and {}% duplicate entries", self.with_data, self.names.len(), self.neg as f32 / (self.pos + self.neg) as f32 );
    }

    pub fn eprint( &self ){
        eprintln!("I have {} kmers and {} genesand {}% duplicate entries", self.with_data, self.names.len(), self.neg as f32 / (self.pos + self.neg) as f32 );
    }

    pub fn get_id( &self, name: String ) -> usize{
        let id = match self.names.get( &name ) {
            Some( id ) => id,
            None => panic!("Gene {name} not defined in the GeneID object"),
        };
        id.clone()
    }

    pub fn get_name( &self, id:usize) -> String{
        let name = match self.names_store.get( id ){
            Some( na ) => na,
            None => panic!("GeneID {id} not defined in the GeneID object"),
        };
        name.to_string()
    }

    fn get_best_gene( &self, genes:&HashMap::<usize, usize>, 
        possible_gene_levels:&HashMap::<usize, usize>, ret: &mut Vec::<usize> ) -> bool{
        ret.clear();
        //eprintln!("I suppose we have a problem with the level being 1 and not 0 here? {genes:?}");
        if genes.len() > 0{
            let mut max = 0;
            let mut id= usize::MAX;
            let mut i = usize::MAX;

            let mut min_level = usize::MAX;
            for level in possible_gene_levels{
                //eprintln!("The optional gene level is {:?} < {min_level}?", level);
                if min_level > *level.1{
                    min_level = *level.1
                }
            }
            //eprintln!("The min gene level is {min_level}");
            for gene_id in genes.keys() {
                // Access the values in both HashMaps using the key
                if let Some(level) = possible_gene_levels.get(gene_id) {
                    if level > &min_level {
                        continue;
                    }
                    if let Some(count) = genes.get(gene_id) {
                        if count > &max{
                            max = *count;
                            id = *gene_id;
                            i = 0;
                        }else if count == &max {
                            i+=1;
                        }
                    }
                }
            }

            //eprintln!("We have obtained the max gene count of {max}");

            if max > 1{
                // there needs to be some reproducibility here!
                return true
            }
            if i == 0{
                // we have exactly one gene with the max count
                //eprintln!("I found one gene: {}", self.names_store[id] );
                ret.push(id);
                return true
            }else {
                // we have more than one gene as a max count - that is inacceptable
                // for the debug I need to know which genes I get here!
                let mut gene_ids= Vec::<usize>::with_capacity(i);
                for (gene_id, count) in genes{
                    if count == &max {
                        gene_ids.push(*gene_id);
                    }
                }
                if gene_ids.len() == 2{
                    if self.names_store[gene_ids[0]] == self.names_store[gene_ids[1]].to_string()+"_int"{
                        ret.push( gene_ids[0] )
                    }
                    if self.names_store[gene_ids[1]] == self.names_store[gene_ids[0]].to_string()+"_int" {
                        ret.push( gene_ids[1] )
                    }
                }
                
                let mut genes = Vec::<String>::with_capacity( gene_ids.len() );
                for gid in gene_ids{
                    ret.push(gid);
                    genes.push( self.names_store[gid].to_string())
                }
                return false;
                //println!("I got a multi mapper ({i}): {:?} -> returning None", genes );
                //multimapper += 1;
            }
            
        }
        return false

    }


    pub fn get(&self, seq: &[u8], report:&mut MappingInfo ) -> Option< usize >{ // gene_id, gene_level
        
        //let mut id:usize;
        let mut genes:HashMap::<usize, usize>= HashMap::new();

        let mut stop: bool;

        //let mut possible_genes = HashMap::<usize, usize>::with_capacity(10);
        let mut possible_gene_levels = HashMap::<usize, usize>::with_capacity(10);

        let mut tool = IntToStr::new(seq.to_vec(), self.tool.kmer_size);
        let mut i = 0;
        tool.from_vec_u8( seq.to_vec() );

        //let mut item = self.tool.next();

        let mut matching_geneids = Vec::<usize>::with_capacity(10);

        // entries is a Option<(u16, u64, usize)
        stop = false;
        'main :while let Some(entries) = tool.next(){

            if self.mapper[entries.0 as usize].has_data() {
                // the 8bp bit is a match
                i +=1;
                //eprintln!("We are at iteration {i}");

                // if matching_geneids.len() == 1 && i > 3 {
                //     //eprintln!("we have one best gene identified! {}",matching_geneids[0]);
                //     return Some( matching_geneids[0] )
                //     //break;
                // }

                //eprintln!("Ill create a new tool here based on seq {}", entries.1);
                let seq_u64 = tool.mask_u64( &entries.1 );

                match &self.mapper[entries.0 as usize].get( &seq_u64 ){
                    Some( gene_id ) => {
                        //eprintln!("Got one: {gene_id:?}");
                        //return ( gene_id );
                        for gid in &gene_id.data{
                            match genes.get_mut( &gid.0) {
                                Some(gene_count) => {
                                    *gene_count +=1;
                                    if *gene_count == 4 {
                                        stop = true;
                                    }
                                },
                                None => {
                                    //eprintln!( "Adding a new gene {} with count 1 here!", gid.0);
                                    genes.insert( gid.0, 1);
                                    possible_gene_levels.insert( gid.0, gid.1 );
                                },
                            };

                        }
                    },
                    None => {
                        //eprintln!("Got one no  gene id in the first run:");
                        match self.mapper[entries.0 as usize].find(&seq_u64 ){
                            Some( gene_id ) => {
                                //eprintln!("But in the second I got one: {gene_id:?}");
                                //return ( gene_id );
                                for gid in &gene_id.data{
                                    match genes.get_mut(&gid.0) {
                                        Some(gene_count) => {
                                            *gene_count +=1;
                                            if *gene_count == 4 {
                                                stop = true;
                                            }
                                        },
                                        None => {
                                            //eprintln!("This is a new gene - I'll insert {} and {}",gid.0, gid.1 );
                                            genes.insert( gid.0, 1);
                                            possible_gene_levels.insert( gid.0, gid.1 );
                                            //eprintln!("I have finished with the indert");
                                        },
                                    };
                                }
                            },
                            None => {        
                                //no_32bp_match +=1;
                            },
                        }
                    },
                }
            }
            // if self.get_best_gene( &genes, &possible_gene_levels, &mut matching_geneids ){
            //     //eprintln!("We found a best gene!");
            //     break 'main;
            // }
            // if stop{
            //     //eprintln!("We did not find a best gene!");
            //     break 'main;
            // }
        }
        // check if there is only one gene //
        let _ = self.get_best_gene( &genes, &possible_gene_levels, &mut matching_geneids );
        let names = self.gene_names_for_ids( &matching_geneids );
        if names.len() > 1 {
            eprintln!("I have multiple genes for this read: {:?}", names);
        }
        //eprintln!("I have {} matching genes", matching_geneids.len());
        if matching_geneids.len() == 0 {
            return None
        }
        // we have not found any gene here!
        // only if all the 8bp have only linked to a single gene...
        if matching_geneids.len() == 1{
            return Some( matching_geneids[0] )
        }



        //let mut max = 0;
        //let mut gid = 0;
        //let mut n =0;
        let mut no_int = Vec::<usize>::with_capacity(5);
        let re = Regex::new(r".*_int").unwrap();
        
        for id in &matching_geneids{
            if ! re.is_match(&self.names_store[*id]) { 
                no_int.push(*id);
            };
        }
        if no_int.len() ==1 {
            //eprintln!("I selected the only RNA seqence {}",self.names_store[ no_int[0] ]);

            return Some(no_int[0])
        }

        let gnames: Vec<String> = matching_geneids.iter().map(|&index| self.names_store[index].to_string()).collect();
        report.write_to_log( format!("Multimapping sequence to {} genes: {:?}", matching_geneids.len(), gnames ));
        report.write_to_log( format!("For the sequence {}", std::str::from_utf8(&seq).expect("Invalid UTF-8") ));
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
        ret.push("n".to_string());
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
        eprintln!("Writing index version {}", version);
        match ofile.buff1.write( &version.to_le_bytes() ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("version could not be written"),
        };

        let seq_len:u64 = self.kmer_len as u64;
        eprintln!("with kmer_len {}", seq_len);
        match ofile.buff1.write( &seq_len.to_le_bytes() ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("seq_len could not be written"),
        };

        let with_data:u64 = self.with_data as u64;
        eprintln!("And a total of {} data entries", with_data);
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
                    for id in 0..tuple.1.data.len(){
                        match ofile.buff1.write( &tuple.1.data[id].0.to_le_bytes() ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                        match ofile.buff1.write( &tuple.1.data[id].1.to_le_bytes() ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                        // match ofile.buff1.write( &tuple.1.significant_bp[id].to_le_bytes() ){
                        //     Ok(_) => (),
                        //     //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                        //     Err(_err) => return Err::<(), &str>("value could not be written"),
                        // };
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
        let mut gene_level:usize;
        let mut len:usize;


        // the real first 8 bytes are the version
        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        //println! ("I got the buff {buff:?}");
        if ! self.version == u64::from_le_bytes( buff_u64 ) as usize{
            panic!("Sorry - this version of the index {} does not match the requirements {}.", u64::from_le_bytes( buff_u64 ), self.version);
        }
        eprintln!("index version {}", self.version);

        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        self.kmer_len = u64::from_le_bytes( buff_u64 ) as usize;
        eprintln!("Index kmer_len: {}", self.kmer_len); // next 8bytes are the kmer_len
        
        ifile.buff1.read_exact(&mut buff_u64).unwrap();
        self.with_data = u64::from_le_bytes( buff_u64 ) as usize;
        eprintln!("reading {} entries",self.with_data );
        
        for i in 0..self.with_data {
            // read in the single mapper elements 
            // first 2 - kmer position
            match ifile.buff1.read_exact(&mut buff_u16){
                Ok(_) => (),
                Err(e) => {
                    eprintln!("On entry {i} we hit this error: {e:?}");
                    break
                }
            };
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
                    ifile.buff1.read_exact(&mut buff_u64).unwrap();
                    gene_level = usize::from_le_bytes( buff_u64 );
                    // ifile.buff1.read_exact(&mut buff_u64).unwrap();
                    // sig_bits= usize::from_le_bytes( buff_u64 );
                    // this should never throw an error.
                    let _ = self.mapper[idx].add( kmer, (gene_id, gene_level), EMPTY_VEC.clone() );

                }
            }
        }

        eprintln!("Reading the gene names\n");
        for name in ifile.buff2.lines() {
            //println!("I read this name: {name:?}");
            //let mut buff = [0_u8 ;8 ];
            match name {
                Ok(n) => {
                    //eprintln!("the name of the gene: {n}\n");
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

    pub fn write_index_txt( &mut self, path: String ) -> Result< (), &str>{
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
        if fs::remove_file(fpath.join("index.1.Index.txt.gz") ).is_ok(){};
        if fs::remove_file(fpath.join("index.1.gene.txt.gz") ).is_ok(){};

        let mut ofile = Ofilesr::new( 1, "index", "Index.txt", "gene.txt",  &path );

        
        let mut count:usize;
        let mut idx:usize;

        let version:u64 = self.version as u64;
        match write!( ofile.buff1, "version: {}\n", version ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("version could not be written"),
        };

        let seq_len:u64 = self.kmer_len as u64;
        match write!( ofile.buff1, "seq_len: {}\n", seq_len ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("seq_len could not be written"),
        };

        let with_data:u64 = self.with_data as u64;
        match write!( ofile.buff1, "with_data: {}\n", with_data ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("with_data could not be written"),
        };

        let mut nucl= String::from ("");
        let mut tuple_idx: usize;

        for i in 0..u16::MAX{
            // write the 8bp kmer as int
            idx = i as usize;
            nucl.clear();
            self.tool.u64_to_str( 8, &(i as u64), &mut nucl ); 
            if self.mapper[idx].has_data(){
                match write!(ofile.buff1, "\ni: {} {} ", i, nucl ){
                    Ok(_) => (), 
                    //Ok(_) => println!("8bp mapper: {:b} binary -> {:?} bytes",i, &i.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("i could not be written"),
                };
                //write the amount of downstream entries
                count = self.mapper[idx].map.len();
                match write!(ofile.buff1, "count {} ", count ){
                    Ok(_) => (), 
                    //Ok(_) => println!("amount of 32 bp second kmers to the first 8bp: {} -> {:?} bytes",count, &count.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("count could not be written"),
                };
                tuple_idx= 0;
                for tuple  in self.mapper[idx].map.iter(){
                    nucl.clear();
                    self.tool.u64_to_str( 32, &tuple.0, &mut nucl ); 
                    match write!(ofile.buff1, "\n\t\t{tuple_idx}: {} resp. {} ",&tuple.0, nucl){
                        Ok(_) => (), 
                        //Ok(_) => println!("the u64 kmer: {} or {:b} binary -> {:?} bytes",tuple.0, tuple.0, &tuple.0.to_le_bytes()  ) ,
                        Err(_err) => return Err::<(), &str>("key could not be written"),
                    };
                    // now we need the NameEntry.len() to to_le_bytes()
                    match write!(ofile.buff1, "n genes: {} ", tuple.1.data.len() ){
                        Ok(_) => (),
                        //Ok( len ) => println!("length of gene list attached to that 32bp kmer: {:?}",len   ) ,
                        Err(_err) => return Err::<(), &str>("value could not be written"),
                    };
                    for id in 0..tuple.1.data.len(){
                        match write!(ofile.buff1, "gene id: {} resp {} ", tuple.1.data[id].0, self.names_store[tuple.1.data[id].0] ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                        match write!(ofile.buff1, "gene level: {}  ", tuple.1.data[id].1 ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                        /* not needed here as we already used this info
                        match write!(ofile.buff1, &tuple.1.significant_bp[id].to_le_bytes() ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };*/
                    }
                    tuple_idx += 1;
                    
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

}





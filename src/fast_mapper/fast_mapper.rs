/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

//use needletail::kmer::Kmers;
//use kmers::naive_impl::Kmer;
use std::collections::BTreeMap;
//use needletail::bitkmer::BitKmer;
use crate::int_to_str::IntToStr;
use std::collections::HashMap;
//use std::collections::HashSet;
use crate::errors::MappingError; // not really errors but enums!
use std::mem;

use rayon::prelude::*; // for the make_index_te_ready

use crate::fast_mapper::mapper_entries::second_seq::SecondSeq;

use crate::ofiles::Ofilesr;
use crate::ifiles::Ifilesr;
use crate::fast_mapper::mapper_entries::MapperEntry;
use crate::fast_mapper::mapper_entries::NameEntry;

use std::path::Path;
use std::io::Write;
use std::fs::{self, DirBuilder};

use std::io::BufRead;
use std::io::Read;
//use std::error::Error;
//use crate::geneids::GeneIds;

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
    debug: bool, //write detailed info for a mapping SLOW
    pub mapper: Vec<MapperEntry>, // the position is the sequence!
    names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names4sparse:  BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names_store: Vec<String>, // store gene names as vector
    pub names_count: Vec<usize>, // store the number of mappers for this gene
    pub max_id: usize,// THIS IS THE ITERATOR FOR THE names_store
    pub last_count: usize, // THIS IS THE ITERATOR FOR THE main names!
    pub with_data: usize,
    pub version: usize, //the version of this
    pub tool: IntToStr,
    mask:u64, // a mask to fill matching sequences to match the index's kmer_len
    pub pos: usize, // when creating the index - this is the amount of times we added a new matching sequence
    pub neg: usize, // when creating the index - this is the amount of times we failed to add a new matching sequence (as we already had this connection)
    /*offset is used to make a fixed index report 'wrong' gene_ids.
    Here is an example: You have two indices and want to collect the data from both into one SingleCellData matrix.
    As these indeices have been created with a basis of 0 you can not collect both data into the same object.
    But using this offset you can give the second index an offset of the max id of the other index.
    Therefore the second index will report ids shifted by that value.
    But it will be totally transperent to the ouside world.
     */
    offset: usize, // the offset is to make a index report 'wrong' ids - say you want
    pub report4: Option<usize>, // report for a gene?
}

// here the functions

impl FastMapper{
    /// kmer_size: how long should the single kmers to search in the sequences be (rec. 9)
    pub fn new( kmer_len:usize, allocate:usize, offset: usize )-> Self {
        //println!("I would add {} new entries here:", u16::MAX);
        let mut mapper: Vec<MapperEntry> = Vec::with_capacity(u16::MAX as usize);
        
        for _i in 0..u16::MAX{
            let b = MapperEntry::new();
            mapper.push( b );
        }
        let names = BTreeMap::<std::string::String, usize>::new();
        let names4sparse = BTreeMap::<std::string::String, usize>::new();
        let names_store = Vec::with_capacity( allocate );
        let names_count = Vec::with_capacity( allocate );
        let max_id = 0;
        let last_count = 0;
        let with_data = 0;
        let spacer = 3;
        let tool = IntToStr::new(b"".to_vec(), kmer_len );
        // this mask will nerver mask ALL bits.
        // But that would also be an error in this whole approach!
        let size = match kmer_len <= 31{
            true => kmer_len,
            false => 31,
        };
        let mask: u64 = (1 << (2 * size)) - 1;
        let neg = 0;
        let pos = 0;

        Self {
            kmer_len,
            spacer,
            debug: false,
            mapper,
            names,
            names4sparse,
            names_store,
            names_count,
            max_id,
            last_count,
            with_data,
            version: 6,
            tool,
            mask,
            pos,
            neg,
            offset, // just to be on the save side here.
            report4:None,
        }
    }

    // Method to calculate memory size
    pub fn memory_size(&self) -> usize {
        let size = mem::size_of::<usize>() * 12 // Sizes of usize fields
                 + mem::size_of::<u64>()        // Size of u64 field
                 + mem::size_of::<Vec<String>>() + self.names_store.iter().map(|s| mem::size_of::<String>() + s.capacity()).sum::<usize>() // Size of Vec<String> and its contents
                 + mem::size_of::<Vec<usize>>() + self.names_count.len() * mem::size_of::<usize>() // Size of Vec<usize> and its contents
                 + mem::size_of::<Vec<MapperEntry>>() + self.mapper.iter().map(|entry| entry.memory_size()).sum::<usize>(); // Size of Vec<MapperEntry> and its contents
        size
    }

    pub fn change_start_id ( &mut self, new_start :usize ){
        if self.offset == 0 {
            self.offset = new_start;
        }else {
            panic!("You try to change the start id twice - not supported!");
        }
        /*self.last_count = new_start;
        for _i in 0..new_start{
            self.names_store.push("PLACEHOLDER".to_string());
            self.names_count.push(0);
        }*/
    }

    /// For the multicore processing I need this object to be mergable
    /// All genes from the other object need to be incorporated and all primaryid - secondaryid compbos need to be collected.
    /// All other classes have to be copied, too. So we need to use our own incorporate_match_combo function.
    pub fn merge( &mut self, other: FastMapper ) {

        for (primary_id, mapper_object) in other.mapper.iter().enumerate(){
            for (key, name_entry) in mapper_object.map.iter(){
                if let Some(entry) = self.mapper[primary_id].map.get_mut( key ){
                    if (! entry.keep) || (! name_entry.keep){
                        entry.keep = false; // nothing more needed
                        continue
                    }
                } else if ! name_entry.keep {
                    let mut this = NameEntry::new( *key );
                    this.keep=false;
                    self.mapper[primary_id].map.insert(*key, this );
                    continue
                }

                for idx in 0..name_entry.classes.len(){   
                    // Now I need to copy the gene names from the other object to my own.
                    let gene_names = other.gene_names_for_ids( &name_entry.classes[idx] );
                    let gene_name = other.gene_names_for_ids( &vec![ name_entry.data[idx].0 ]);
                    // And get my own ids for these names
                    //let gene_ids = self.ids_for_gene_names_mut( &gene_names );
                    let _ = self.incorporate_match_combo( primary_id, name_entry.key , gene_name[0].clone(), gene_names );
                }           
            }
        }

    }

    pub fn debug(&mut self, dbg: Option<bool>) -> bool{
        if let Some(val) = dbg { 
            self.debug = val 
        }
        self.debug
    }


    /// You get a Vector og gene names for a vector of IDs.
    /// This function takes into account that we might have an offset here.
    pub fn gene_names_for_ids( &self, ids:&Vec<usize> ) -> Vec<String> {
        let mut ret = Vec::<String>::with_capacity( ids.len() );
        for id in ids.iter(){
            let fixed = *id - self.offset;
            ret.push( self.names_store[fixed].clone())
        }
        ret
    }

    /// how many genes does this index describe?
    pub fn get_gene_count( &self) -> usize {
        self.names.len()
    }

    /// Returns all possible gene names this index object describes
    pub fn get_all_gene_names (&self) -> Vec::<String> {
        self.names.keys().map(|name| name.to_string()).collect()
    }
    
/*
    fn gene_names_for_intern_ids( &self, ids:&Vec<usize> ) -> Vec<String> {
        let mut ret = Vec::<String>::with_capacity( ids.len() );
        for id in ids.iter(){
            let fixed = *id;
            ret.push( self.names_store[fixed].clone())
        }
        ret
    }
*/

    pub fn ids_for_gene_names( &self, ids:&Vec<String> ) -> Vec<usize> {
        let mut ret = Vec::<usize>::with_capacity( ids.len() );
        for name in ids.iter(){
            let id = self.get_id( name.to_string() );

            ret.push( id );
        }
        ret
    }


    pub fn ids_for_gene_names_mut( &mut self, ids:&Vec<String> ) -> Vec<usize> {
        let mut ret = Vec::<usize>::with_capacity( ids.len() );
        for name in ids.iter(){
            if ! self.names.contains_key( name ){
                self.names.insert( name.clone(), self.last_count );
                self.names_store.push( name.clone() );
                self.names_count.push( 0 );
                self.last_count += 1;
            }
            let id = self.get_id( name.to_string() );
            ret.push( id  );
        }
        ret
    }

    fn intern_ids_for_gene_names( &mut self, ids:&Vec<String> ) -> Vec<usize> {
        let mut ret = Vec::<usize>::with_capacity( ids.len() );
        for name in ids.iter(){
            if ! self.names.contains_key( name ){
                self.names.insert( name.clone(), self.last_count );
                self.names_store.push( name.clone() );
                self.names_count.push( 0 );
                self.last_count += 1;
            }
            let id = self.get_intern_id( name.to_string() );
            ret.push( id );
        }
        ret
    }

    fn get_extern_id_from_intern( &self, id: usize ) -> usize{
        id + self.offset
    }

    pub fn get_id( &self, name: String ) -> usize{
        let id = match self.names.get( &name ) {
            Some( id ) => id + self.offset,
            None => panic!("Gene {name} not defined in the GeneID object"),
        };
        id
    }

    pub fn extern_id_for_gname ( &self, name:&str ) -> Option<usize> {
        self.names.get( name ).map(|id| id + self.offset)
    }

    fn get_intern_id( &self, name: String ) -> usize{
        let id = match self.names.get( &name ) {
            Some( id ) => id ,
            None => panic!("Gene {name} not defined in the GeneID object"),
        };
        *id
    }

    pub fn get_name( &self, id:usize) -> String{
        let name = match self.names_store.get( id - self.offset ){
            Some( na ) => na,
            None => panic!("GeneID {id} not defined in the GeneID object"),
        };
        name.to_string()
    }

    /*fn local_gene_id (&self, id:usize) -> usize {
        id - self.offset
    }*/

    pub fn incorporate_match_combo( &mut self, initial_match:usize, secondary_match:SecondSeq, name:String, ids: Vec<String> )  ->Result<(), &str>{
        // this is cheked beforehand in the merge function
        /*let tmp = NameEntry::new( secondary_match.clone() );
        if let Some(entry) = self.mapper[initial_match].get( &tmp ){
            if ! entry.pass {
                // useless to add something else here
                return Ok ( () )
            }
        }*/
        if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.last_count );
            self.names_store.push( name.clone() );
            self.names_count.push( 0 );
            self.last_count += 1;

            //println!("fast_mapper: I have {} -> {}",  name, self.names.len()-1);

        }
        let gene_id = self.get_id( name.to_string() );
        let classes =  self.ids_for_gene_names_mut( &ids );

        if ! self.mapper[initial_match].has_data() {
            // will add in the next step so
            self.with_data +=1;
        }
        if self.mapper[initial_match].add( secondary_match, ( gene_id, 0), classes ){
            self.pos += 1
        }


        Ok( () )
    }

    pub fn fill_sequence_with_a( self, seq:&u64 ) -> u64{
        
        seq | (!self.mask & (0b11 << (2 * self.kmer_len )))
    }


    /*pub fn add_small(&mut self, seq: &Vec<u8>, name: std::string::String, class_ids: Vec<String> ) -> usize{


        let mut classes_vec = vec![name.clone()];
        classes_vec.extend( class_ids);

        let classes =  self.intern_ids_for_gene_names( &classes_vec );
        let gene_id = classes[0];
        //eprintln!("fast_mapper: I have {} -> {}", name, gene_id);

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
        let mut n = 50;

        while let Some(entries) = self.tool.next_small(){
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

            if entries.1.1 < 5_u8{
                // too short second entry
                break;
            }
            if entries.0 > 65534{
                continue; // this would be TTTTTTTT and hence not use that!
            }
            if ! self.mapper[entries.0 as usize].has_data() {
                // will add in the next step so
                self.with_data +=1;
            }
            if entries.1.1  < self.kmer_len as u8 {
                //eprintln!("Too small sequence");
                continue;
            }
            if self.mapper[entries.0 as usize].add( entries.1, (gene_id, 0), classes.clone() ){
                //eprintln!("I have added a sequence! - getting the gene_id {gene_id}");
                self.pos += 1;
                self.names_count[gene_id] +=1;
                i+=1;
            }else {
                //eprintln!("NOT - I have added a sequence!");
                self.neg +=1
            }
            
        }

        //eprintln!("fast_mapper adds this genes: {} of length {} ({}, {})", name, seq.len(), self.pos, self.neg );
        
        //println!("I have now {i} mappers for the gene {name}");
        if i == 0 {
            //eprint!(".");
            println!("gene {} ({:?}) does not get an entry in the fast_mapper object! - too short? {} bp", &name.as_str(), classes, seq.len() );
        }
        /*else {
            println!("I added {i} mappper entries for gene {name}");
        }*/

        //self.with_data += i;
        //println!( "{} kmers for gene {} ({})", total, &name.to_string(), gene_id );
        gene_id
    }*/


    /// The add function returns the intern id. 
    pub fn add(&mut self, seq: &Vec<u8>, name: std::string::String, class_ids: Vec<String> ) -> usize{


        let mut classes_vec = vec![name.clone()];
        classes_vec.extend( class_ids);

        let classes =  self.intern_ids_for_gene_names( &classes_vec );
        let gene_id = classes[0];

        

        let mut i =0;
        //let mut n = 20;

        // let mut short = String::from("");
        // let mut long = String::from("");

        self.tool.from_vec_u8( seq.to_vec() );
        for entries in self.tool.by_ref(){

            //println!("add -> I got the u16 {} and the u64 as {}",entries.0, entries.1);
            //tmp = tmp + format!("it {i} ");
            //i +=1;
            //n -=1;
            //if n == 0{
            //    break;
            //}
            //self.tool.print_second_seq( entries.0, entries.1 );

            // index_read, longer_read, longer length
            // short.clear();
            // long.clear();
            // self.tool.u16_to_str( 8, &entries.0, &mut short);
            // self.tool.u64_to_str( entries.1.1.into(), &entries.1.0, &mut long);
            // println!("I got some entry [short {short}, long {long}, length {}", entries.1.1 );

            if entries.1.1  < 10{
                //eprintln!("Too small sequence {}",entries.1.1 );
                continue;
            }
            if entries.0 > 65534{
                continue; // this would be TTTTTTTT and hence not use that!
            }

            //println!("Added?!");

            // And here check if the longer and shorter entry together are complex enough


            if ! self.mapper[entries.0 as usize].has_data() {
                // will add in the next step so
                self.with_data +=1;
            }
            
            if self.mapper[entries.0 as usize].add( entries.1, (gene_id, 0), classes.clone() ){
                //println!("Addin the id {} to this index!", self.tool.u64_to_string( 8,&(entries.0 as u64) ) );
                self.pos += 1;
                self.names_count[ gene_id] +=1;
                if i < self.names_count[gene_id]{
                    i = self.names_count[gene_id];
                }
                i+=1;
            }else {
                //println!("NOT added this sequence! {:#b}+{} -> {gene_id} & {classes:?} ",entries.0, entries.1);
                self.neg +=1
            }
            
        }

        //eprintln!("fast_mapper adds this genes: {} of length {} ({}, {})", name, seq.len(), self.pos, self.neg );
        
        //println!("I have now {i} mappers for the gene {name}");
        if i == 0 {
            //eprint!(".");
            println!("add -> gene {} ({:?}) does not get an entry in the fast_mapper object! - too short or a duplicate? {} bp", &name.as_str(), classes, seq.len() );
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
    pub fn make_index_te_ready_single( &mut self ) {

        eprintln!("Before make_index_te_ready:");
        self.eprint();
        for mapper_entry in self.mapper.iter_mut() {
            if mapper_entry.has_data() {
                //eprintln!("I collapse the data for the mapper entry {:?}", mapper_entry);
                mapper_entry.collapse_classes();
                //eprintln!("is it collapes now\n\n {:?}\n\n?", mapper_entry);
            }
        }

        self.with_data = 0;
        for mapper_entry in &self.mapper{
            if mapper_entry.has_data() {
                self.with_data += 1;
            }
        }

        eprintln!("\nAfter:");
        self.eprint();
    }


    pub fn make_index_te_ready(&mut self) {
        eprintln!("Before make_index_te_ready:");
        self.eprint();

        self.mapper.par_iter_mut().for_each(|mapper_entry| {
            if mapper_entry.has_data() {
                mapper_entry.collapse_classes();
            }
        });

        // Recalculate self.with_data
        self.with_data = self.mapper.par_iter().filter(|entry| entry.has_data()).count();

        eprintln!("\nAfter:");
        self.eprint();
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
        if !self.names.is_empty() {
            println!("I have {} kmers for {} genes with {}% duplicate entries", self.with_data, self.names.len(), self.neg as f32 / (self.pos + self.neg) as f32 );
            println!("gene names like '{}'", self.names_store[self.names_store.len()-1]);
            println!("gene_ids range from {:?} to {:?}\n", self.names.values().min(), self.names.values().max());
        }else {
            println!("This index is empty\n");
        }    
    }

    pub fn eprint( &self ){
        if !self.names.is_empty() {
            eprintln!("I have {} kmers for {} genes with {}% duplicate entries", self.with_data, self.names.len(), self.neg as f32 / (self.pos + self.neg) as f32 );
            eprintln!("gene names like '{}'", self.names_store[self.names_store.len()-1]);
            eprintln!("gene_ids range from {:?} to {:?}\n", self.names.values().min(), self.names.values().max());
        }else {
            eprintln!("This index is empty\n");
        }
    }

    fn mean( numbers: &[f32] ) -> Option<f32>{
        if numbers.is_empty() {
            None
        }else {
            Some(numbers.iter().sum::<f32>() / numbers.len() as f32)
        }
    }

    fn get_best_gene( &self, genes:&HashMap::<(usize, usize), Vec<f32> >, ret: &mut Vec::<usize> ) -> bool{
        ret.clear();
        if genes.is_empty(){
            return false
        }

        /*let mut report = false;
        if genes.len() > 2{
            report = true;
            eprintln!("Lots of genes matched here?: {genes:?}");
            for  (gene_id, _level) in genes.keys(){
                eprint!(" {}", self.names_store[*gene_id] );
            }
            eprint!("\n");
        }*/
        let mut min_genes = Vec::new();
        let mut mean_value: Option<f32> = None;
        let mut most_matches = usize::MAX;

        for ((gene, _class), vec) in genes {
            if let Some(mean_val) = Self::mean(vec) {
                if mean_val < mean_value.unwrap_or(f32::INFINITY) {
                    //println!("I got a new mean value: {mean_val} and gene {gene}");
                    min_genes.clear(); // Clear the previous min_genes
                    min_genes.push(*gene);
                    mean_value = Some(mean_val);
                    most_matches = vec.len();
                } else if mean_val == mean_value.unwrap_or(f32::INFINITY) {
                    //println!("more gene with the mean values!: {mean_val} and gene {gene}");
                    min_genes.push(*gene);
                }
            }
        }

        //println!("I have some best gene(s): {min_genes:?} with a minimum mapping value of {mean_value:?} and a total of {most_matches} mapping events");
        /*if genes.len() > 1 {
            eprintln!("I got this as most matches {most_matches} for that data {genes:?}");
        }*/
        
        if most_matches < 3 {
            return false
        }

        if let Some(value) = mean_value{
            if value > 0.4 {
                return false
            }
        }else {
            return false
        }
                // Process min_genes as needed
        for gene_id in min_genes {
            //println!("I push the gene {gene_id}");
            ret.push(gene_id);
            
        }
        true
    }

    pub fn get_strict(&self, seq: &[u8], tool: &mut IntToStr ) -> Result< Vec<usize>, MappingError >{ // gene_id, gene_level
        
        //let mut id:usize;
        let mut genes:HashMap::<(usize, usize), Vec<f32>>= HashMap::new();

        //let mut possible_genes = HashMap::<usize, usize>::with_capacity(10);
        //let mut tool = IntToStr::new(seq.to_vec(), self.tool.kmer_size);
        let mut i = 0;
        tool.from_vec_u8( seq.to_vec() );
        let mut na = 0;
        //let mut item = self.tool.next();

        let mut matching_geneids = Vec::<usize>::with_capacity(10);

        // entries is a Option<(u16, u64, usize)
        for entries in tool.by_ref(){
            if entries.0 == 65535{
                continue;
            }

            if na == 10 && genes.is_empty(){
                match self.debug {
                    true => println!("get - I have no match after 10 oligo pairs - giving up"),
                    false => (),
                };
                return Err(MappingError::NoMatch)
            }
            if self.mapper[entries.0 as usize].has_data() {
                // the 8bp bit is a match
                i +=1;
                //println!("We are at iteration {i}");

                match &self.mapper[entries.0 as usize].get( &entries.1 ){
                    Some( (gene_ids, nw_value) ) => {
                        //println!("And we got a match! ({gene_id:?})");
                        for name_entry in gene_ids {
                            for gid in &name_entry.data{
                                let ext_gid = self.get_extern_id_from_intern( gid.0 );
                                genes.entry( (ext_gid, gid.1) ).or_insert(Vec::with_capacity(10)).push(*nw_value);
                            }
                        }
                    },
                    None => {
                        na +=1;
                    },
                }
            }

        }
        match self.debug {
            true => {
                let names: Vec<String> = genes.keys()
                    .map(|&(id, _)| self.get_name(id) )
                    .collect();
                println!("get - testing {i} oligos I collected these genes: {genes:?} with these names: {names:?}")
            },
            false => (),
        };
        if genes.is_empty() {
            return Err(MappingError::NoMatch)
        }
        // check if there is only one gene //
        if self.get_best_gene( &genes, &mut matching_geneids ){
            //println!("And I got a match! ({matching_geneids:?})");
            if matching_geneids.len() == 1 {
                return Ok( matching_geneids )
            }else {
                return Err(MappingError::MultiMatch)
            }
        }
        Err(MappingError::NoMatch)
    }


    pub fn get(&self, seq: &[u8], tool: &mut IntToStr  ) -> Result< Vec<usize>, MappingError >{ // gene_id, gene_level
        
        //let mut id:usize;

        if let Ok(gene_id) = self.get_strict( seq, tool ) { 
            return Ok(gene_id) 
        }

        //println!("\n\nStart fast_mapper::get()");

        let mut genes:HashMap::<(usize, usize), Vec<f32>>= HashMap::new();

        //let mut possible_genes = HashMap::<usize, usize>::with_capacity(10);

        tool.from_vec_u8( seq.to_vec() );

        //let mut item = self.tool.next();

        let mut matching_geneids = Vec::< usize>::with_capacity(10);
        let mut na = 0;
        // entries is a Option<(u16, u64, usize)

        let mut i = 0;
        //let mut tmp = "".to_string();
        /*
        let mut small_seq :String = Default::default();
        let mut large_seq: String = Default::default();*/

        for entries in tool.by_ref(){

            //tmp = tmp + &format!("iteration {i} ");
            //i +=1;

            if na == 10 && matching_geneids.is_empty() {
                // This is obviousely not exisitingn in the index.
                //println!("{tmp} - stopped after 10 failed matches\n");
                match self.debug {
                    true => println!("get - I have no match after 10 oligo pairs - giving up"),
                    false => (),
                };
                return Err(MappingError::NoMatch)
            }
            if entries.0 == 65535{
                //tmp += "useless u16;\n";
                continue;
            }

            if self.mapper[entries.0 as usize].has_data() {
                // the 8bp bit is a match
                //println!("The first id {} is known to this index!", tool.u64_to_string( 8,&(entries.0 as u64) ) );

                //eprintln!("We are at iteration {i}");
                if entries.1.1 < 20_u8{
                    //tmp += "too short u64;\n";
                    //will never map as the mapper fails every sequence below 20 bp - too many false positives!
                    continue
                }
                i +=1;

                if (i > 12) && genes.is_empty() {
                    //tmp += "failing after 4 failed matches and no match at all;\n";
                    //println!("First four matches have not given anything! breaking!");
                    break;
                }
                /*small_seq.clear();
                large_seq.clear();
                tool.u16_to_str( 8, &entries.0, &mut small_seq);
                tool.u64_to_str( entries.1.1.into(), &entries.1.0, &mut large_seq);
                println!("We are on iteration {i} with seq {}-{} ( {}-{})", small_seq, large_seq, &entries.0, &entries.1);
                */
                

                // if matching_geneids.len() == 1 && i > 3 {
                //     //eprintln!("we have one best gene identified! {}",matching_geneids[0]);
                //     return Some( matching_geneids[0] )
                //     //break;
                // }

                //eprintln!("Ill create a new tool here based on seq {}", entries.1);

                match &self.mapper[entries.0 as usize].get( &entries.1 ){

                    Some( (gene_ids, nw_value) ) => {
                        //tmp += &format!("And we identified some gene(s): {gene_ids:?}\n");
                        //println!("Got some gene ids (one?): {:?}", gene_ids);
                        for name_entry in gene_ids {
                            for gid in &name_entry.data{
                                let ext_gid = self.get_extern_id_from_intern( gid.0 );
                                genes.entry( (ext_gid, gid.1) ).or_insert(Vec::with_capacity(10)).push(*nw_value);
                            }
                        }
                        
                    },
                    None => {
                        
                        //println!("Got no gene id in the first run get() run -> find() instead:");
                        match self.mapper[entries.0 as usize].find(  &entries.1 ){
                            Some( (gene_ids, nw_value) ) => {
                                //println!("But in the second I got one: {gene_ids:?}");
                                //tmp += &format!("And we identified some gene(s) with less stingency: {gene_ids:?}\n");
                                for name_entry in &gene_ids{
                                    for gid in &name_entry.data{
                                        let ext_gid = self.get_extern_id_from_intern( gid.0 );
                                        genes.entry( (ext_gid, gid.1) ).or_insert(Vec::with_capacity(10)).push(nw_value);
                                    }
                                }
                            },
                            None => {  
                                na +=1;
                                //tmp +="no gene detceted!\n";
                                //eprintln!("And none in the find() either.");
                                //no_32bp_match +=1;
                            },
                        }
                        
                    },
                }
            }else { // my mapper has no data!?
                //tmp += &format!("the mapper had not indexed that u16 {}\n ",entries.0);
               // println!("The first id {} is not known to this index!", tool.u64_to_string( 8,&(entries.0 as u64) ) );
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

        //tmp += &format!("at the end we collected this gene info: {genes:?} ");
        match self.debug {
            true => {
                let names: Vec<String> = genes.keys()
                    .map(|&(id, _)| self.get_name(id) )
                    .collect();
                println!("get - testing {i} oligos I collected these genes: {genes:?} with these names: {names:?}")
            },
            false => (),
        };
        

        if genes.is_empty() {
            //println!("{tmp} No gene found");
            //println!("Nothing found - really nothing?!");
            return Err(MappingError::NoMatch)
        }
        
        /*let bad_gene = "ADA";
        let bad = self.get_id( bad_gene.to_string()) ;*/
        
        // check if there is only one gene //
        if self.get_best_gene( &genes, &mut matching_geneids ){
            // Cd3e <- keep that to fiond this place back
            /*if matching_geneids[0] == bad {
                println!("read mapping to {} - should not happen here?: {:?}\n{:?}", bad_gene, self.gene_names_for_ids( &matching_geneids ),String::from_utf8_lossy(seq) );
                //println!("This is our total matching set: {:?}", genes);
            }*/
            //println!("gene {matching_geneids:?} detected");
            if matching_geneids.len() == 1 {
                return Ok( matching_geneids )
            }else {
                return Err(MappingError::MultiMatch)
            }
            
        }
        //eprintln!("We got a mutlimatcher?! {matching_geneids:?}, {genes:?}");
        //println!("{tmp} multimapper?!");
        Err(MappingError::MultiMatch)
    }

    pub fn to_header( &self ) -> std::string::String {

        let mut ret= Vec::<std::string::String>::with_capacity( self.names.len() +4 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for obj in self.names.keys() {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push(  obj.to_string() ) ;
        }
        ret.push("Most likely name".to_string());
        ret.push("Faction total".to_string());
        ret.push("dist to nr.2 [%max]".to_string());
        "CellID\t".to_owned()+&ret.join("\t")
    }

    

    pub fn to_header_n( &self, names: &Vec<String> ) -> std::string::String {
        let mut ret= Vec::<std::string::String>::with_capacity( names.len() +5 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for name in names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( name.to_string() ) ;
        }
        ret.push("AssignedSampleName".to_string());
        ret.push("FractionTotal".to_string());
        ret.push("n".to_string());
        ret.push("dist to nr.2 [%max]".to_string());
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
                count = self.mapper[idx].with_data();
                match ofile.buff1.write( &count.to_le_bytes() ){
                    Ok(_) => (), 
                    //Ok(_) => println!("amount of 32 bp second kmers to the first 8bp: {} -> {:?} bytes",count, &count.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("count could not be written"),
                };
                for (_, tuple) in self.mapper[idx].map.iter(){
                    if ! tuple.keep {
                        continue;
                    }
                    // this will write the u64 and the u8 in one go.
                    match ofile.buff1.write( &tuple.key.to_le_bytes() ){
                        Ok(_) => (), 
                        //Ok(_) => println!("the u64 kmer: {} or {:b} binary -> {:?} bytes",tuple.0, tuple.0, &tuple.0.to_le_bytes()  ) ,
                        Err(_err) => return Err::<(), &str>("key could not be written"),
                    };
                    // now we need the NameEntry.len() to to_le_bytes()
                    match ofile.buff1.write( &tuple.data.len().to_le_bytes() ){
                        Ok(_) => (),
                        //Ok( len ) => println!("length of gene list attached to that 32bp kmer: {:?}",len   ) ,
                        Err(_err) => return Err::<(), &str>("value could not be written"),
                    };
                    for id in 0..tuple.data.len(){
                        match ofile.buff1.write( &tuple.data[id].0.to_le_bytes() ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                        match ofile.buff1.write( &tuple.data[id].1.to_le_bytes() ){
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
        let mut buff_seq_entry = [0_u8 ;9 ];

        let mut kmer:SecondSeq;
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
        
        eprintln!("Reading the gene names\n");
        for name in ifile.buff2.lines() {
            //println!("I read this name: {name:?}");
            //let mut buff = [0_u8 ;8 ];
            match name {
                Ok(n) => {
                    //eprintln!("the name of the gene: {n}\n");
                    if ! self.names.contains_key( &n ){
                        self.names.insert( n.clone(), self.last_count );
                        self.names_store.push( n.clone() );
                        self.names_count.push(0); // will be populated later.
                        self.max_id += 1;
                        self.last_count +=1;
                    }
                },
                Err(e) => panic!("I could not read from the genes table {e:?}"),
            };
        }
           
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
                ifile.buff1.read_exact(&mut buff_seq_entry).unwrap();
                kmer = SecondSeq::from_le_bytes( buff_seq_entry ).unwrap();
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
                    if self.mapper[idx].add( kmer, (gene_id, gene_level), EMPTY_VEC.clone() ){
                        self.names_count[gene_id] +=1;
                    }

                }
            }
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

        let mut ofile = Ofilesr::new( 1, "index", "Index.txt", "gene.txt.gz",  &path );

        
        let mut count:usize;
        let mut idx:usize;

        let version:u64 = self.version as u64;
        match writeln!( ofile.buff1, "version: {}", version ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("version could not be written"),
        };

        let seq_len:u64 = self.kmer_len as u64;
        match writeln!( ofile.buff1, "seq_len: {}", seq_len ){
            Ok(_) => (),
            Err(_err) => return Err::<(), &str>("seq_len could not be written"),
        };

        let with_data:u64 = self.with_data as u64;
        match writeln!( ofile.buff1, "with_data: {}", with_data ){
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
                match write!(ofile.buff1, "\ni: {} ({:016b}) {} ", i, i,nucl ){
                    Ok(_) => (), 
                    //Ok(_) => println!("8bp mapper: {:b} binary -> {:?} bytes",i, &i.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("i could not be written"),
                };
                //write the amount of downstream entries
                count = self.mapper[idx].with_data();
                match write!(ofile.buff1, "count {} ", count ){
                    Ok(_) => (), 
                    //Ok(_) => println!("amount of 32 bp second kmers to the first 8bp: {} -> {:?} bytes",count, &count.to_le_bytes()  ) ,
                    Err(_err) => return Err::<(), &str>("count could not be written"),
                };
                tuple_idx= 0;
                for (_, tuple)  in self.mapper[idx].map.iter(){
                    if ! tuple.keep {
                        continue;
                    }
                    //nucl.clear();
                    //self.tool.u64_to_str( 32, &tuple.key.0, &mut nucl ); 
                    match write!(ofile.buff1, "\n\t\t{tuple_idx}: {} ",&tuple.key){
                        Ok(_) => (), 
                        //Ok(_) => println!("the u64 kmer: {} or {:b} binary -> {:?} bytes",tuple.0, tuple.0, &tuple.0.to_le_bytes()  ) ,
                        Err(_err) => return Err::<(), &str>("key could not be written"),
                    };
                    // now we need the NameEntry.len() to to_le_bytes()
                    match write!(ofile.buff1, "n genes: {} ", tuple.data.len() ){
                        Ok(_) => (),
                        //Ok( len ) => println!("length of gene list attached to that 32bp kmer: {:?}",len   ) ,
                        Err(_err) => return Err::<(), &str>("value could not be written"),
                    };
                    for id in 0..tuple.data.len(){
                        match write!(ofile.buff1, "gene id: {} resp {} ", tuple.data[id].0, self.names_store[tuple.data[id].0] ){
                            Ok(_) => (),
                            //Ok(_) => println!("\tgene_id: {} -> {:?} bytes",id, &id.to_le_bytes()  ) ,
                            Err(_err) => return Err::<(), &str>("value could not be written"),
                        };
                        match write!(ofile.buff1, "gene level: {}  ", tuple.data[id].1 ){
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




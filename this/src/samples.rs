/// Samples Strings are a little short to be efficiently mapped using the 
// fast_mapper class. Inseated this should be handled by a different approach
// Likely a simple 8bp matching array.

use std::collections::BTreeMap;
//use std::collections::HashSet;
use std::collections::HashMap;
//use kmers::naive_impl::Kmer;
use crate::int_to_str::IntToStr;
use crate::traits::Index;


#[derive(Debug)]
pub struct SampleIds{
    data: Vec<usize>, // store gene names as vector
    
}

impl SampleIds{
	pub fn new( )-> Self {
		let data = Vec::new();
		Self{
			data,
			
		}
	}

	pub fn add(&mut self, value:usize ) ->usize {
		self.data.push(value);
		self.data.len()
	}

	pub fn entries( &self ) -> Vec<usize>{
		self.data.clone()
	}

}


#[derive(Debug)]
pub struct Samples{
    pub mapper: Vec<SampleIds>, // the position is the sequence!
    pub names : BTreeMap<std::string::String, usize>, // gene name and gene id
    pub names_store: Vec<String>, // store gene names as vector
    tool:IntToStr,
    max_id:usize,
    with_data:usize,
    // the sample_id entries that have been 
    neg:usize,
    min_counts: usize, // how many matches does a sample need to be reported?
}


impl Index for Samples{
	fn new(  _kmer_len:usize, _allocate:usize)-> Self {
		let mut mapper: Vec<SampleIds> = Vec::with_capacity(u16::MAX as usize);
  
        for _i in 0..u16::MAX{
            let b = SampleIds::new();
            mapper.push( b );
        }

        let names = BTreeMap::<std::string::String, usize>::new();
        let names_store = Vec::with_capacity( 4 );
        let max_id = 0;
        // init a tool with totally useless values (not needed here)
        let tool = IntToStr::new( b"AACCGGTT".to_vec() , 16 );
        Self {
        	mapper,
        	names,
        	names_store,
        	tool,
        	max_id,
        	with_data :0,
        	neg: 0,
        	min_counts:4,
        }
    }

    fn names(&self) -> Vec<String>{
    	self.names_store.clone()
    }

	fn change_start_id ( &mut self, new_start :usize ){
        self.max_id = new_start;
        for _i in 0..new_start{
            self.names_store.push("na".to_string());
        }
    }

    fn add(&mut self, seq: &Vec<u8>, name: std::string::String, _class_ids: Vec<String> ) -> usize{
    	if ! self.names.contains_key( &name ){
            self.names.insert( name.clone(), self.max_id );
            self.names_store.push( name.clone() );
            self.max_id += 1;
        }
        self.tool.from_vec_u8( seq.to_vec() );
        let sample_id = self.get_id( name.to_string() );
        let mut i = 0;
        while let Some(entry) = self.tool.iter_8bp(){
        	match self.mapper[entry as usize].add( sample_id ){
        		1 => self.with_data +=1,
        		2 => self.neg +=1,
        		_ => {},
        	};
        	i += 1;
        }
        return i
    }


    fn to_header_n(&self, names: &Vec<String> ) -> String {
    	let mut ret = names.to_vec();
    	ret.extend_from_slice(&["AsignedSampleName".to_string(), "FractionTotal".to_string(), "n".to_string()]);
    	format!("CellID\t{}", ret.join("\t"))
	}


    fn get( &self, seq: &[u8] ) -> Option<usize>{
    	
    	self.tool.from_vec_u8( seq.to_vec() );
    	let mut possible = HashMap::<usize, usize>::new();

    	while let Some(seq_id) = self.tool.iter_8bp_4bp_hops(){
    		// so now lets find them
    		for sample_id in self.mapper[seq_id as usize].entries(){
    			*possible.entry(sample_id).or_insert(0) += 1;
    		}
    	}

    	let max_val = possible.values().max().unwrap_or(&0);
    	if max_val < &self.min_counts{
    		return None
    	}
    	let mut good = Vec::<usize>::with_capacity( possible.len() );
    	for  (sample_id, matches ) in possible.iter(){
    		if matches == max_val{
    			good.push(*sample_id);
    		}
    	}
    	if good.len() == 1 {
    		return Some(good[0])
    	}
    	return None
    }

    fn get_id( &self, name: String ) -> usize{
        let id = match self.names.get( &name ) {
            Some( id ) => id,
            None => panic!("Gene {name} not defined in the GeneID object"),
        };
        id.clone()
    }

    fn print( &self ){
        if self.names_store.len() > 0 {
            println!("I have {} kmers for {} samples with {}% duplicate entries", self.with_data, self.names.len(), self.neg as f32 / (self.with_data + self.neg) as f32 );
            println!("gene names like '{}'", self.names_store[0]);
        }else {
            println!("This index is empty");
        }    
    }

    fn names_len( &self ) -> usize {
    	self.names_store.len()
    }

    fn names4sparse(&self) -> Vec<String> {
    	self.names_store.clone()
    }

    // not needed here
    fn reset_names4sparse( &mut self ){
    	panic!("Not implemented")
    }
    fn add_2_names4sparse( &mut self, name:&str ){
    	panic!("Not implemented")
    }
    fn max_id(&self) ->usize {
    	panic!("Not implemented")
    }

}

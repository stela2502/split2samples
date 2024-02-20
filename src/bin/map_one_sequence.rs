use clap::Parser;

use rustody::fast_mapper::FastMapper;
use rustody::int_to_str::IntToStr;

use std::path::Path;
use std::collections::HashMap;

use needletail::parse_fastx_file;

static EMPTY_VEC: Vec<String> = Vec::new();


#[derive(Parser)]
#[clap(version = "0.0.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the dna to seach for
    #[clap(short, long)]
    dna: String,
    /// the fasta database containing the genes
    #[clap(short, long)]
    expression: Option<String>,
    /// the fasta database containing the antibody tags
    #[clap(short, long)]
    antibody: Option<String>,
    /// a pre-defined index folder produced by the cerateIndex scipt
    #[clap(short, long)]
    index: Option<String>,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the path to write the text index to
    #[clap(short, long)]
    outpath: String,
    
}


// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    
    let opts: Opts = Opts::parse();

	// let mut cell_umi:HashSet<u128> = HashSet::new();
    //let mut genes :GeneIds = GeneIds::new(gene_kmers); // split them into 9 bp kmers
    let mut genes :FastMapper = FastMapper::new( 32, 100_000, 0 ); // split them into 9 bp kmers
    let mut samples :FastMapper = FastMapper::new( 32, 10_000 , 0 );
    let mut antibodies :FastMapper = FastMapper::new( 32, 10_000, 0  );

    let mut gene_count = 600;
    
    let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

    if let Some(i) = &opts.index {
    	println!("Loading index from path {i}");
    	match genes.load_index( i.to_string() ){
    		Ok(_r) => (),
    		Err(e) => panic!("Failed to load the index {e:?}")
    	}
    	genes.print();
    	gene_count = genes.names.len();
    	
    }

    let mut gene_names = Vec::new();
    for gname in &genes.names_store {
    	gene_names.push( gname.to_string());
    }

    let mut gene_names:Vec<String> = Vec::with_capacity(gene_count);

    for gene in genes.names.keys() {
    	gene_names.push(gene.to_string());
    }
    let mut ab_names:Vec<String> = Vec::with_capacity(30);

    let mut seq_temp:Vec::<u8>;

    if let Some(ex) = &opts.expression {
    	if Path::new(&ex).exists(){

	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
	                		if ! genes.names.contains_key(  id ){
	                			seq_temp = seqrec.seq().to_vec();
	                			//seq_temp.reverse();
		                    	genes.add( &seq_temp, id.to_string(), EMPTY_VEC.clone() );
	                    		gene_names.push( id.to_string() );
	                    	}
	                    	//genes2.add_unchecked( &seqrec.seq(), id.to_string() );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	eprintln!("Expression file could not be read - ignoring")
	    }
    }

    eprintln!("Changing the expression start gene id to {}", genes.last_count );
	antibodies.change_start_id( genes.last_count);

    if let Some(ab) = &opts.antibody {

	    if Path::new(&ab).exists(){


	   		let mut ab_file = parse_fastx_file(ab).expect("valid path/file");
	    	while let Some(ab_record) = ab_file.next() {
		        let seqrec = ab_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
	                		seq_temp = seqrec.seq().to_vec();
	                		//seq_temp.reverse();
		                    antibodies.add_small( &seq_temp, id.to_string(), EMPTY_VEC.clone() );
	                    	ab_names.push( id.to_string() );
	                    	//gene_names.push( id.to_string() );
	                    	//genes2.add_unchecked( &seqrec.seq(), id.to_string() );
	                	};
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	eprintln!("Antibody file could not be read - ignoring")
	    }

	}
	let mut id = 0;
	samples.change_start_id( antibodies.last_count );
	    if  opts.specie.eq("human") {
	        // get all the human sample IDs into this.
	        // GTTGTCAAGATGCTACCGTTCAGAGATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG
	        //                        b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG"
	        // Seams they have introduced new sample ids, too....
	        let sequences = [ b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG",
	        b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 
	        b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC",
	        b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT",
	        b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC",
	        b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG" ];

	        for seq in sequences{
	        	//seq.reverse();
	        	let mut seq_ext = b"GTTGTCAAGATGCTACCGTTCAGAG".to_vec();
	        	seq_ext.extend_from_slice( seq );
	        	samples.add( &seq_ext, format!("Sample{id}"),EMPTY_VEC.clone() );
	        	//sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }
	    }
	    else if opts.specie.eq("mouse") {
	        // and the mouse ones
	        let sequences = [b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 
	        b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
	        b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
	        b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
	        b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
	        b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC"];

	        for seq in sequences{
	        	//seq.reverse();
	        	//let mut seq_ext = b"GTTGTCAAGATGCTACCGTTCAGAG".to_vec();
	        	//seq_ext.extend_from_slice( seq );
	        	//samples.add_small( &seq_ext, format!("Sample{id}"),EMPTY_VEC.clone() );
	        	samples.add( &seq.to_vec(), format!("Sample{id}"),EMPTY_VEC.clone() );
	        	//sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }

	    } else {
	        println!("Sorry, but I have no primers for species {}", &opts.specie);
	        std::process::exit(1)
	    }

	// to understand the mapping a little better I now need the index written to a file
	let _ = genes.write_index_txt(  format!("{}/genes_index", &opts.outpath) );
	let _ = antibodies.write_index_txt(  format!("{}/antibodies_index", &opts.outpath) );
	let _ = samples.write_index_txt(  format!("{}/samples_index", &opts.outpath) );

	eprintln!("Ill check this dna: '{}'", &opts.dna );

	let mut collected:HashMap::<usize, usize>= HashMap::new();

	let mut ok = match antibodies.get( &opts.dna.as_bytes(), &mut tool , 0.2, 5){
		Ok(gene_ids) =>{
        	//eprintln!("gene id {gene_id:?} seq {:?}", String::from_utf8_lossy(&data[i].1) );
        	//eprintln!("I got an ab id {gene_id}");
        	
        	//println!("I found an antibody match: {:?}", gene_ids );
        	for gid in gene_ids {
        			match collected.get_mut( &gid) {
        				Some(gene_count) => {
        					*gene_count +=1;

        				},
        				None => {
		                    //eprintln!( "Adding a new gene {} with count 1 here!", gid.0);
		                    collected.insert( gid.clone(), 1);
		                },
		            };
		    }
		    true
		},
		Err(_) => {
			false
		}
	};
	if ! ok{
		ok = match samples.get( &opts.dna.as_bytes(),  &mut tool, 0.18, 15){
			Ok(gene_ids) =>{
	        	//eprintln!("gene id {gene_id:?} seq {:?}", String::from_utf8_lossy(&data[i].1) );
	        	//eprintln!("I got an ab id {gene_id}");
	        	
	        	//println!("I found an sample match: {:?}", gene_ids );
	        	for gid in gene_ids {
	        			match collected.get_mut( &gid) {
	        				Some(gene_count) => {
	        					*gene_count +=1;

	        				},
	        				None => {
			                    //eprintln!( "Adding a new gene {} with count 1 here!", gid.0);
			                    collected.insert( gid.clone(), 1);
			                },
			            };
			    }
			    true
			},
			Err(_) => {
				false
			}
		};
	}

	if ! ok{
		let _ = match genes.get( &opts.dna.as_bytes(),  &mut tool, 0.4, 5 ){
			Ok(gene_ids) =>{
	        	//eprintln!("gene id {gene_id:?} seq {:?}", String::from_utf8_lossy(&data[i].1) );
	        	//eprintln!("I got an ab id {gene_id}");
	        	
	        	//println!("I found an gene match: {:?}", gene_ids );
	        	for gid in gene_ids {
	        			match collected.get_mut( &gid) {
	        				Some(gene_count) => {
	        					*gene_count +=1;

	        				},
	        				None => {
			                    //eprintln!( "Adding a new gene {} with count 1 here!", gid.0);
			                    collected.insert( gid.clone(), 1);
			                },
			            };
			    }
			    true
			},
			Err(_) => {
				false
			}
		};
	}

	let mut names =Vec::<String>::with_capacity( collected.len());
	for key in collected.keys(){
		if key < &antibodies.last_count{
			names.push( antibodies.names_store[ key - antibodies.offset ].to_string());
		}
		else if key < &samples.last_count{
			names.push( samples.names_store[ key - samples.offset ].to_string());
		}
		else {
			names.push( genes.names_store[ key - genes.offset ].to_string());
		}
	}


    println!("I have collected these genes: {names:?}" );



    println!("Finished");
}

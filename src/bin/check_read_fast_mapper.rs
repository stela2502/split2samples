use clap::Parser;

use rustody::fast_mapper::FastMapper;
use rustody::int_to_str::IntToStr;

use std::path::Path;
use std::collections::HashMap;

use needletail::parse_fastx_file;

static EMPTY_VEC: Vec<String> = Vec::new();

/// Uses the logic from the quantify_ scripts to test exactly one read sequence (just the fasta part without the name line)

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
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
    /// Mapper how many times does a 40bp read combo need to match to any given gene to be reported (default=1)
    #[clap( long)]
    min_matches: Option<usize>,
    /// What is the highest acceptable needleman wush inspired cut off (default 0.5)
    #[clap( long)]
    highest_nw_val: Option<f32>,
    /// What is the highest acceptable humming distance to even run NW (default 0.6)
    #[clap( long)]
    highest_humming_val: Option<f32>,
    
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
    genes.debug( Some(true) );
    samples.debug( Some(true) );
    antibodies.debug( Some(true) );

    if let Some(min_matches) = opts.min_matches{
        genes.set_min_matches( min_matches );
        samples.set_min_matches( min_matches );
        antibodies.set_min_matches( min_matches );
        println!("Setting the mapper min_matches to {min_matches}")
    }
    if let Some(highest_nw_val) = opts.highest_nw_val{
       genes.set_highest_nw_val( highest_nw_val );
       samples.set_highest_nw_val( highest_nw_val );
       antibodies.set_highest_nw_val( highest_nw_val );
       println!("Setting the mapper highest_nw_val to {highest_nw_val}")
    }
    if let Some(highest_humming_val) = opts.highest_humming_val{
    	genes.set_highest_humming_val( highest_humming_val);
    	samples.set_highest_humming_val( highest_humming_val);
        antibodies.set_highest_humming_val( highest_humming_val);
        println!("Setting the mapper highest_humming_val to {highest_humming_val}")
    }
    
    let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

    if let Some(i) = &opts.index {
    	println!("Loading index from path {i}");
    	match genes.load_index( i ){
    		Ok(_r) => (),
    		Err(e) => panic!("Failed to load the index {e:?}")
    	}
    	genes.print();
    	
    }

    

    let mut seq_temp:Vec::<u8>;

    if let Some(ex) = &opts.expression {
    	if Path::new(&ex).exists(){

	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
                			seq_temp = seqrec.seq().to_vec();
                			//seq_temp.reverse();
	                    	genes.add( &seq_temp, id, id, EMPTY_VEC.clone() );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	eprintln!("Expression file could not be read - ignoring")
	    }
    }

	//let gene_names = genes.get_all_gene_names();

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
		                    antibodies.add( &seq_temp, id, id, EMPTY_VEC.clone() );
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
	//let ab_names = antibodies.get_all_gene_names();
	let mut id = 0;
	samples.change_start_id( genes.last_count+ antibodies.last_count );
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
	        	samples.add( &seq_ext, &format!("SampleTag{id:02}_hs"), &format!("SampleTag{id:02}_hs"), EMPTY_VEC.clone() );
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
	        	samples.add( &seq.to_vec(), &format!("SampleTag{id:02}_mm"), &format!("SampleTag{id:02}_mm"), EMPTY_VEC.clone() );
	        	//sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }

	    } else {
	        println!("Sorry, but I have no primers for species {}", &opts.specie);
	        std::process::exit(1)
	    }

		genes.make_index_te_ready();
		antibodies.make_index_te_ready();
		samples.make_index_te_ready();

	// to understand the mapping a little better I now need the index written to a file
	let _ = genes.write_index_txt(  &format!("{}/genes_index", &opts.outpath) );
	let _ = antibodies.write_index_txt(  &format!("{}/antibodies_index", &opts.outpath) );
	let _ = samples.write_index_txt(  &format!("{}/samples_index", &opts.outpath) );


	println!("Ill check this dna: '{}'", &opts.dna );

	let mut collected:HashMap::<usize, usize>= HashMap::new();

	println!("\nChecking antibodies:");
	let mut ok = match antibodies.get( opts.dna.as_bytes(), &mut tool ){
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
		                    collected.insert( gid, 1);
		                },
		            };
		    }
		    true
		},
		Err(_) => {
			false
		}
	};
	println!("I have collected these genes: {collected:?}" );
	if ! ok{
		println!("\nChecking samples:");
		ok = match samples.get( opts.dna.as_bytes(),  &mut tool ){
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
			                    collected.insert( gid, 1);
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
	println!("I have collected these genes: {collected:?}" );

	if ! ok{
		println!("\nChecking genes:");
		let _ = match genes.get( opts.dna.as_bytes(),  &mut tool ){
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
			                    collected.insert( gid, 1);
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


    println!("I have collected these genes: {collected:?}" );



    println!("Finished");
}

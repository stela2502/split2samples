use clap::Parser;

use rustody::fast_mapper::FastMapper;
use rustody::int_to_str::IntToStr;

use std::path::Path;
use std::time::SystemTime;

use needletail::parse_fastx_file;

static EMPTY_VEC: Vec<String> = Vec::new();


#[derive(Parser)]
#[clap(version = "0.0.3", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    sequence: String,
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

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

	// let mut cell_umi:HashSet<u128> = HashSet::new();
    //let mut genes :GeneIds = GeneIds::new(gene_kmers); // split them into 9 bp kmers
    let mut genes :FastMapper = FastMapper::new( 32, 100_000 ); // split them into 9 bp kmers
    let mut samples :FastMapper = FastMapper::new( 32, 10_000  );
    let mut antibodies :FastMapper = FastMapper::new( 32, 10_000  );

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


	/// to understand the mapping a little better I now need the index written to a file
	genes.write_index_txt(  format!("{}/genes.index.txt", &opts.outpath) );
	antibodies.write_index_txt(  format!("{}/antibodies.index.txt", &opts.outpath) );
	samples.write_index_txt(  format!("{}/samples.index.txt", &opts.outpath) );

	eprintln!("Ill check this sequence: '{}'", &opts.sequence );
	let mut ok = match antibodies.get_strict( &opts.sequence.as_bytes(), &mut tool ){
		Some(gene_id) =>{
        	//eprintln!("gene id {gene_id:?} seq {:?}", String::from_utf8_lossy(&data[i].1) );
        	//eprintln!("I got an ab id {gene_id}");
        	if gene_id.len() == 1 {
        		println!("I found an antibody match: {:?}", gene_id );
        		true
        	}else {
        		false
        	}

        },
        None => {
        	false
        }
    };
    if ! ok{
    	ok = match samples.get_strict( &opts.sequence.as_bytes(),  &mut tool ){
            Some(gene_id) =>{
            	//println!("sample ({gene_id:?}) with {:?}",String::from_utf8_lossy(&data[i].1) );
            	//eprintln!("I got a sample umi id {umi}");
            	if gene_id.len() == 1 {
            		println!("I found a sample id match: {:?}", gene_id );
                    true
                }else {
                	false
                }
                
            },
            None => {
				false
            }
        };
    }

    if ! ok{
    	
        match genes.get( &opts.sequence.as_bytes(),  &mut tool ){
        	Some(gene_id) =>{
                if gene_id.len() == 1 {
                    println!("I found a gene id match: {:?}", gene_id );
                }else {
                	println!("I found multiple genes ids match: {:?}", gene_id );
                }
                
            },
            None => {
            }
        };
    }



    println!("Finished");
}

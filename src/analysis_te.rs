use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
//use std::collections::HashSet;

use crate::singlecelldata::SingleCellData;
use crate::singlecelldata::cell_data::GeneUmiHash;
//use crate::geneids::GeneIds;
use crate::fast_mapper::FastMapper;
use crate::int_to_str::IntToStr;
//use crate::traits::te_index;
use crate::errors::{ FilterError, MappingError }; // not really errors but enums!
use crate::cellids::CellIds;
use crate::cellids10x::CellIds10x;
use crate::traits::CellIndex;

use crate::mapping_info::MappingInfo;
//use std::io::BufReader;
//use flate2::write::GzDecoder;

//use this::last5::Last5;

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
//use std::time::Duration;

use needletail::errors::ParseError;

use std::path::PathBuf;
use std::fs::File;
use std::path::Path;
use std::io::Write;

use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;

use crate::ofiles::{Ofiles, Fspot};
//use std::fs::write;


static EMPTY_VEC: Vec<String> = Vec::new();


fn mean_u8( data:&[u8] ) -> f32 {
    let this = b"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

    let mut total = 0;
    let mut sum = 0;
    for ent in data {
        for (i, entry) in this.iter().enumerate() {
            if ent == entry {
                sum += i;
                total += 1;
                break;
            }
        }
        //total += 1;
    }
    sum as f32 / total as f32
}

/// the analysis calss is a wrapper around my 'old' quantify_rhapsody main funtion.
/// I started it to easily process multiple fastq files in a row.
pub struct AnalysisTE{
	expr_index_obj: FastMapper,
	samples: FastMapper,
	te_index_obj: FastMapper,
	cells: Box<dyn CellIndex + Sync>,
	gex: SingleCellData,
	sample_names:Vec<String>,
	gene_names:Vec<String>,
	te_names:Vec<String>,
	num_threads:usize,
}


impl AnalysisTE{


	pub fn new(gene_kmers:usize, version:String, te_index:Option<String>, 
		expression_index:Option<String>, specie:String, num_threads:usize, exp:&str  ) -> Self{
		//let sub_len = 9;
	    //let mut cells = SampleIds::new( sub_len );// = Vec::with_capacity(12);
	    //cells.init_rhapsody( &opts.specie );

	    // let mut cell_umi:HashSet<u128> = HashSet::new();
	    //let mut expr_index_obj :GeneIds = GeneIds::new(gene_kmers); // split them into 9 bp kmers
	    
	    let mut te_index_obj :FastMapper = FastMapper::new( gene_kmers, 10_000, 0  );

	    let mut gene_count = 600;
	    
	    if let Some(i) = te_index {
	    	println!("Loading te_index from path {i}");
	    	match te_index_obj.load_index( i ){
	    		Ok(_r) => (),
	    		Err(e) => panic!("Failed to load the te_index {e:?}")
	    	}
	    	te_index_obj.print();
	    	gene_count = te_index_obj.names.len();
	    	
	    }

	    let mut expr_index_obj :FastMapper = FastMapper::new( gene_kmers, 100_000, gene_count ); // split them into 9 bp kmers

	    if let Some(expr_path) = expression_index {

		    if Path::new(&expr_path).exists(){
		    	println!("Loading expression index from path {expr_path}");
		    	match expr_index_obj.load_index( expr_path ){
		    		Ok(_r) => (),
		    		Err(e) => panic!("Failed to load the te_index {e:?}")
		    	}
		    	expr_index_obj.print();
		    	//gene_count = expr_index_obj.names.len();
		    }else {
		    	eprintln!("expression_index file could not be read - ignoring")
		    }

		} 

		let te_names = te_index_obj.names.keys().cloned().collect::<Vec<String>>();
		let gene_names = expr_index_obj.names.keys().cloned().collect::<Vec<String>>();

	    //  now we need to get a CellIDs object, too
	    let cells:  Box::<dyn CellIndex + Sync> = match exp {
	    	"bd" => Box::new(CellIds::new( &version)),
	    	"10x" => Box::new(CellIds10x::new( &version )),
	    	_ => panic!("Only 'bd' or '10x' systems are supported in the exp option")
	    };

	    // that is a class to strore gene expression data.
	    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
	    // and umi's are u64 again
	    // here I need the cell kmer site.
	    let gex = SingleCellData::new( num_threads );

	    let mut sample_names:Vec<String> = Vec::with_capacity(12);

	    let mut id = 1;

	   	let mut samples :FastMapper = FastMapper::new( gene_kmers, 10_000, te_index_obj.names.len() +  expr_index_obj.names.len()  );

	    if  specie.eq("human") {
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
	        	//let mut seq_ext = b"GTTGTCAAGATGCTACCGTTCAGAG".to_vec();
	        	//seq_ext.extend_from_slice( seq );
	        	samples.add_small( &seq.to_vec(), format!("Sample{id}"),EMPTY_VEC.clone() );
	        	sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }
	    }
	    else if specie.eq("mouse") {
	        // and the mouse ones
	        let sequences = [
	        b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 
	        b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
	        b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
	        b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
	        b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
	        b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC"];

	        /*
	        " GTTGTCAAGATGCTACCGTTCAGAG ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC"
			"GTTGTCAAGATGCTACCGTTCAGGG ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCCAAAA"
			"                          AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC"
			"GTTGTCAAGATGCTATCGTTCAGAG AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTCAAAA"
			"                          TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT"
			"GTTGTCAAGATG-TACCGTTCAGAG TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCTAAAAA"
			"                          ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC"
			"GA GTCAAGATGCTACCGTTCAGAG ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCCAAAAA"
	        */
	        for seq in sequences{
	        	//seq.reverse();
	        	//let mut seq_ext = b"GTTGTCAAGATGCTACCGTTCAGAG".to_vec();
	        	//seq_ext.extend_from_slice( seq );
	        	//samples.add_small( &seq_ext, format!("Sample{id}"),EMPTY_VEC.clone() );
	        	samples.add_small( &seq.to_vec(), format!("Sample{id}"),EMPTY_VEC.clone() );
	        	sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }

	    } else {
	        println!("Sorry, but I have no primers for species {}", specie);
	        std::process::exit(1)
	    }

	    println!("After indexing all fastq files we have the following indices:");
	    println!("the TE index:");
		te_index_obj.print();
		println!("the expression index:");
		expr_index_obj.print();
		println!("and the sample id index:");
		samples.print();
		
	    
		Self{
			expr_index_obj,
			samples,
			te_index_obj,
	//		expr_index_obj2,
			cells,
			gex,
			sample_names,
			gene_names: gene_names.to_vec(),
			te_names: te_names.to_vec(),
			num_threads,
		}
	}

	pub fn write_index(&mut self, path:&String ){
		self.expr_index_obj.write_index( path.to_string() ).unwrap();
		self.expr_index_obj.write_index_txt( path.to_string() ).unwrap();
	}


    fn quality_control(   read1:&SequenceRecord, read2:&SequenceRecord, min_quality:f32, pos: &[usize;8], min_sizes: &[usize;2]  ) -> Result<(), FilterError> {
        
        // the quality measurements I get here do not make sense. I assume I would need to add lots of more work here.
        //println!("The read1 mean quality == {}", mean_u8( seqrec.qual().unwrap() ));

        if mean_u8( read1.qual().unwrap() ) < min_quality {
            //println!("filtered R1 by quality" );
            return Err(FilterError::Quality);
        }
        if mean_u8( read2.qual().unwrap() ) < min_quality {
            //println!("filtered R2 by quality" );
            return Err(FilterError::Quality);
        }

        // totally unusable sequence
        // as described in the BD rhapsody Bioinformatic Handbook
        // GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf (google should find this)
        if read1.seq().len() < min_sizes[0] {
        	//println!("R1 too short (< {}bp)\n", min_sizes[0] );
            return Err(FilterError::Length);
        }
        if read2.seq().len() < min_sizes[1] {
        	//println!("R2 too short (< {}bp)\n",  min_sizes[1]);
            return Err(FilterError::Length);
        }
        //let seq = seqrec.seq().into_owned();
        for nuc in &read2.seq()[pos[6]..pos[7]] {  
            if *nuc ==b'N'{
            	//println!("N nucs\n");
                return Err(FilterError::Ns);
            }
        }

        /*
        // I have huge problems with reds that contain a lot of ploy A in the end
        // None of them match using blast - so I think we can savely just dump them.
        let mut bad = 0;

        let mut last_a = false;

        for nuc in &read2.seq()[ (read2.seq().len()-30).. ] {
        	if *nuc ==b'A'{
        		if last_a{
        			bad += 1;
        		}
        		last_a = true;
        	}else {
        		last_a = false;
        	}
        }
        if bad >15 {
        	//eprintln!("Sequence rejected due to high polyA content \n{:?}", String::from_utf8_lossy(&read2.seq()) );
        	return Err(FilterError::PolyA)
        }
        */

        Ok(())
    }

    fn analyze_paralel( &self, data:&[(Vec<u8>, Vec<u8>)], report:&mut MappingInfo, pos: &[usize;8] ) -> SingleCellData{
    	


        // first match the cell id - if that does not work the read is unusable
        //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
        let mut gex = SingleCellData::new( self.num_threads );
        let mut ok : bool;

        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

        // lets tag this with the first gene I was interested in: Cd3e
        //let goi_id = self.expr_index_obj.get_id("ADA".to_string());

        for i in 0..data.len() {

        	match &self.cells.to_cellid( &data[i].0 ){
	            Ok( (cell_id, umi) ) => {
	            	//let tool = IntToStr::new( data[i].0[(pos[6]+add)..(pos[7]+add)].to_vec(), 32 );
	            	report.cellular_reads +=1;

	            	// now I have three possibilites here:
	            	// an expression_index tag match
	            	// a sample id match
	            	// or a mRNA match
	            	// And of casue not a match at all


	            	ok = match &self.expr_index_obj.get( &data[i].1, &mut tool ){
	                    Ok(gene_id) =>{
	                    	//eprintln!("gene id {gene_id:?} seq {:?}", String::from_utf8_lossy(&data[i].1) );
	                    	//eprintln!("I got an ab id {gene_id}");
	                    	report.iter_read_type( "expression reads" );
	                    	
                    		let data = GeneUmiHash( gene_id[0], *umi);
	                    	if ! gex.try_insert( 
	                        	&(*cell_id as u64),
	                        	data,
	                        	report
	                        ){
	                        	report.pcr_duplicates += 1 
	                        }
	                        true
	                    },
	                    Err(MappingError::NoMatch) => {
	                    	false
	                    },
	                    Err(MappingError::MultiMatch) => {
	                    	// this is likely not mapping to anyting else if we alredy have mult matches here!
	                    	report.multimapper +=1;
	                    	report.no_data +=1;
	                    	continue
	                    }
	                };

	                if ! ok{
	                	ok = match &self.samples.get_strict( &data[i].1,  &mut tool ){
		                    Ok(gene_id) =>{
		                    	//println!("sample ({gene_id:?}) with {:?}",String::from_utf8_lossy(&data[i].1) );
		                    	//eprintln!("I got a sample umi id {umi}");
		                    	report.iter_read_type( "sample reads" );
		                    	
	                    		let data = GeneUmiHash( gene_id[0], *umi);
		                        if ! gex.try_insert( 
		                        	&(*cell_id as u64),
		                        	data,
		                        	report
		                        ) { 
		                        	report.pcr_duplicates += 1 
		                        }
		                        true 
		                    },
		                    Err(MappingError::NoMatch) => {
		                    	false
		                    },
		                    Err(MappingError::MultiMatch) => {
		                    	// this is likely not mapping to anyting else if we alredy have mult matches here!
		                    	report.multimapper +=1;
		                    	report.no_data +=1;
		                    	continue
		                    }
		                };
	                }

	                if ! ok{
	                	
		                match &self.te_index_obj.get_strict( &data[i].1,  &mut tool ){
		                	Ok(gene_id) =>{
		                		// will now always be exactly ONE!
		                		report.iter_read_type( "antibody reads" );
		                		//println!("gene found strict, but not lax?! ({gene_id:?}) with {:?}",String::from_utf8_lossy(&data[i].1) );
		                    	//eprintln!("I got a sample umi id {umi}");
			                    
		                    	let data = GeneUmiHash( gene_id[0], *umi);
		                        if ! gex.try_insert( 
		                        	&(*cell_id as u64),
		                        	data,
		                        	report
		                        ){
		                        	report.pcr_duplicates += 1 
		                        }
		                        /*if gene_id[0] == goi_id{
		                        	println!( "UMI for gene id {goi_id} and cell {cell_id} is {umi}" );
		                        }*/
		                    },
		                    Err(MappingError::NoMatch) => {
		                    	// I want to be able to check why this did not work
		                    	report.write_to_ofile( Fspot::Buff1, 
		                    		format!(">Cell{cell_id} no gene detected\n{:?}\n", &data[i].1) 
		                    	);
								
								report.write_to_ofile( Fspot::Buff2, 
									format!(">Cell{cell_id} no gene detected\n{:?}\n", &data[i].0 ) 
								);
		                    	report.no_data +=1;
		                    },
		                    Err(MappingError::MultiMatch) => {
		                    	// this is likely not mapping to anyting else if we alredy have mult matches here!
		                    	report.multimapper +=1;
		                    	report.no_data +=1;
		                    }
		                };

		            }
	            },
	            Err(_err) => {
	            	// this is fucked up - the ids are changed!
	            	let _ =report.write_to_ofile( Fspot::Buff1, 
	            		format!(">No Cell detected\n{:?}\n{:?}\n{:?}\n{:?}\n", 
	            			&data[i].0, 
	            			&data[i].0[pos[0]..pos[1]],
	            			&data[i].0[pos[2]..pos[3]],
	            			&data[i].0[pos[4]..pos[5]],
	            		)
	            	);

                	report.write_to_ofile( Fspot::Buff2, 
                		format!(">No Cell detected\n{:?}\n", &data[i].1 )
                	);

	                report.no_sample +=1;
	                continue
	            }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
       		};
    	}
        gex
    }


    /// Analze BPO Rhapsody data in a paralel way.
    pub fn parse_parallel(&mut self,  f1:&str, f2:&str,  
    	report:&mut MappingInfo,pos: &[usize;8], min_sizes: &[usize;2], outpath: &str, max_reads:usize, chunk_size:usize ){

    	println!("I am using {} cpus", self.num_threads);

    	report.start_counter();

    	let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = match parse_fastx_file(&f1) {
        	Ok(reader) => reader,
        	Err(err) => {
            	panic!("File 'reads' {f1} Error: {}", err);
        	}
   		};
        let mut readefile = match parse_fastx_file(&f2) {
        	Ok(reader) => reader,
        	Err(err) => {
            	panic!("File 'file' {f2} Error: {}", err);
        	}
   		};
        let m = MultiProgress::new();
        let pb = m.add(ProgressBar::new(5000));
        pb.set_style(spinner_style);

        //let reads_perl_chunk = 1_000_000;
        //eprintln!("Starting with data collection");
        let mut good_reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity( chunk_size * self.num_threads );
        let mut good_read_count = 0;

        'main: while let (Some(record1), Some(record2)) = (&readereads.next(), &readefile.next())  {
        	if report.total > max_reads{
        		break 'main
        	}
        	if good_read_count < chunk_size*self.num_threads {
        		report.total += 1;
	            let read2 = match record2{
	                Ok( res ) => res,
	                Err(err) => {
	                    eprintln!("could not read from R2:\n{err}");
	                    continue 'main;
	                }
	            };
	            let read1 = match record1{
	                Ok(res) => res,
	                Err(err) => {
	                    eprintln!("could not read from R1:\n{err}");
	                    continue 'main;
	                }
	            };

        		match Self::quality_control( &read1, &read2, report.min_quality, pos, min_sizes){
        			Ok(()) => {
        				good_read_count +=1;
        				good_reads.push( (read1.seq().to_vec(), read2.seq().to_vec() ) );
        			},
        			/*Err(FilterError::PolyA)=> {
        				report.poly_a +=1;
        				continue 'main;
        			},*/
        			Err(FilterError::Quality) => {
        				report.quality +=1;
        				continue 'main;
        			},
        			Err(FilterError::Length) => {
        				report.length +=1;
        				continue 'main;
        			},
        			Err(FilterError::Ns) => {
        				report.n_s +=1;
        				continue 'main;
        			}
        		};
    		}
    		else {
    			report.stop_file_io_time();
    			//eprintln!("Mapping one batch");
    			pb.set_message( format!("mapping reads - {}", report.log_str() ) );
            	//eprintln!("I have {} lines of data and {} threads", good_reads.len(), self.num_threads);

    			good_read_count = 0;
		    	let total_results: Vec<(SingleCellData, MappingInfo)> = good_reads
			        .par_chunks(good_reads.len() / self.num_threads + 1) // Split the data into chunks for parallel processing
		    	    .map(|data_split| {
		    	    	//eprintln!("I am processing {} lines of data", data_split.len() );
		    	    	// Get the unique identifier for the current thread
		            	let thread_id = thread::current().id();
		            
		            	// Convert the thread ID to a string for use in filenames or identifiers
		            	let thread_id_str = format!("{:?}",thread_id );
		            	let ofile = Ofiles::new( 1, &("Umapped_with_cellID".to_owned()+&thread_id_str), "R2.fastq.gz", "R1.fastq.gz",  outpath );
		            	let log_file_str = PathBuf::from(outpath).join(
		        			format!("Mapping_log_{}.txt",thread_id_str )
		        		);
			    
			    		let log_file = match File::create( log_file_str ){
			        			Ok(file) => file,
			        			Err(err) => {
			           			panic!("thread {thread_id_str} Error: {err:#?}" );
			        		}
			    		};

			    		let mut rep = MappingInfo::new( Some(log_file), report.min_quality, report.max_reads, Some(ofile) );
			    		rep.write_to_log( format!("I am processing {} lines of data", data_split.len() ));
			            // Clone or create a new thread-specific report for each task
			    	    let res = self.analyze_paralel(&data_split, &mut rep, pos );
			    	    (res, rep)

		    	    }) // Analyze each chunk in parallel
		        .collect(); // Collect the results into a Vec
		        //eprintln!("sum up temp results");

		        report.stop_multi_processor_time();

			    for gex in total_results{
			    	self.gex.merge(&gex.0);
			       	report.merge( &gex.1 );
			    }
			    //eprintln!("Collecting more reads");
			    good_reads.clear();
			    //println!("{}", report.log_str());
			    report.stop_single_processor_time();
			}
			report.log(&pb);
		}
		

	    if good_read_count > 0{
	    	report.stop_file_io_time();
	    	pb.set_message( format!("mapping reads - {}", report.log_str() ) );
	    	// there is data in the good_reads
	    	let total_results: Vec<(SingleCellData, MappingInfo)> = good_reads
		        .par_chunks(good_reads.len() / self.num_threads + 1) // Split the data into chunks for parallel processing
	    	    .map(|data_split| {
	    	    	// Get the unique identifier for the current thread
	            	let thread_id = thread::current().id();
	            
	            	// Convert the thread ID to a string for use in filenames or identifiers
	            	let thread_id_str = format!("{:?}",thread_id );
	            	let ofile = Ofiles::new( 1, &("Umapped_with_cellID".to_owned()+&thread_id_str), "R2.fastq.gz", "R1.fastq.gz",  outpath );
	            	let log_file_str = PathBuf::from(outpath).join(
	        			format!("Mapping_log_{}.txt",thread_id_str )
	        		);
		    
		    		let log_file = match File::create( log_file_str ){
		        			Ok(file) => file,
		        			Err(err) => {
		           			panic!("thread {thread_id_str} Error: {err:#?}" );
		        		}
		    		};
		    		let mut rep = MappingInfo::new( Some(log_file), report.min_quality, report.max_reads, Some(ofile) );
		            // Clone or create a new thread-specific report for each task
		    	    let res = self.analyze_paralel(data_split, &mut rep, pos );
		    	    (res, rep)

	    	    }
	    	    ) // Analyze each chunk in parallel
	        .collect(); // Collect the results into a Vec

	        report.stop_multi_processor_time();

	        for gex in total_results{
	        	self.gex.merge(&gex.0);
	        	report.merge( &gex.1 );
	        }
	        pb.set_message( report.log_str().clone() );
	    }

	    report.stop_single_processor_time();

	    let log_str = report.log_str();

        pb.finish_with_message( log_str );

    }

/*
    pub fn parse(&mut self,  f1:String, f2:String,  report:&mut MappingInfo,pos: &[usize;8], min_sizes: &[usize;2]  ){

    	let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
    	.unwrap()
    	.tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&f1).expect("valid path/file");
        let mut readefile = parse_fastx_file(&f2).expect("valid path/file");
        let mut tool = IntToStr::new(b"AAGGCCTT".to_vec(), 32);
        let m = MultiProgress::new();
        let pb = m.add(ProgressBar::new(5000));
        pb.set_style(spinner_style);
        //pb.set_prefix(format!("[{}/?]", i + 1));

        //let report_gid = self.expr_index_obj.get_id( "Sample1".to_string() );
        'main: while let Some(record2) = readefile.next() {
        	if let Some(record1) = readereads.next() {
        		report.total += 1;

        		let seqrec2 = match record2{
        			Ok( res ) => res,
        			Err(err) => {
        				eprintln!("could not read from R2:\n{err}");
        				continue 'main;
        			}
        		};
        		let seqrec1 = match record1{
        			Ok(res) => res,
        			Err(err) => {
        				eprintln!("could not read from R1:\n{err}");
        				continue 'main;
        			}
        		};

                //eprintln!("Processing line {}", report.total);
            	//let id_fq = String::from_utf8_lossy(&seqrec.seq()).to_string();
            	//let id_fq1 = String::from_utf8_lossy(&seqrec1.seq()).to_string();
            	//eprintln!("These are my seqs:\n{:?},\n{:?}", id_fq,id_fq1);

                // the quality measurements I get here do not make sense. I assume I would need to add lots of more work here.
                //println!("The read1 mean quality == {}", mean_u8( seqrec.qual().unwrap() ));
                match Self::quality_control( &seqrec1, &seqrec2, report.min_quality, pos, min_sizes){
        			Ok(()) => {
        				//good_read_count +=1;
        			},
        			/*Err(FilterError::PolyA)=> {
        				report.poly_a +=1;
        				continue 'main;
        			},*/
        			Err(FilterError::Quality) => {
        				report.quality +=1;
        				continue 'main;
        			},
        			Err(FilterError::Length) => {
        				report.length +=1;
        				continue 'main;
        			},
        			Err(FilterError::Ns) => {
        				report.n_s +=1;
        				continue 'main;
        			}
        		};
                
                let mut ok: bool;
                // first match the cell id - if that does not work the read is unusable
                //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                match &self.cells.to_cellid( &seqrec1.seq() ){
                	Ok( (cell_id, umi) ) => {
		            	report.cellular_reads +=1;

		            	// now I have three possibilites here:
		            	// an expression_index tag match
		            	// a sample id match
		            	// or a mRNA match
		            	// And of casue not a match at all

		            	ok = match &self.samples.get_strict( &seqrec2.seq(),  &mut tool ){
	                		Some(gene_id) =>{
	                			report.iter_read_type( "sample reads" );
		                    	//println!("sample ({gene_id:?}) with {:?}",String::from_utf8_lossy(&seqrec.seq()) );

		                    	//eprintln!("Got a samples match! {gene_id}");
		                    	//eprintln!( "{:?} got {gene_id} resp {:?} ", String::from_utf8_lossy( &seqrec.seq() ), &self.samples.names_store[*gene_id] );
		                    	if gene_id.len() == 1 {
		                    		let data = GeneUmiHash( gene_id[0], *umi);
			                    	self.gex.try_insert( 
			                    		&(*cell_id as u64),
			                    		data,
			                    		report
			                    		);
			                    	true
			                    }
			                    else {
			                    	false
			                    }
		                    },
		                    None => {
		                    	false
		                    }
		                };
		            	

		                if ! ok{
		                	ok = match &self.te_index_obj.get_strict( &seqrec2.seq(),  &mut tool ){
			            		Some(gene_id) =>{
			            			report.iter_read_type( "expression_index reads" );
			                    	if gene_id.len() == 1 {
			                    		let data = GeneUmiHash( gene_id[0], *umi);
				                    	self.gex.try_insert( 
				                    		&(*cell_id as u64),
				                    		data,
				                    		report
				                    		);
				                    	true
				                    }
				                    else {
				                    	false

				                    }
			                    },
			                    None => {
			                    	false
			                    }
			                }; 	
			            }

			            if ! ok{
			            	// multimapper are only allowed for expression reads - 
			            	// if either sample ids or AB reads would have this it would be a library design fault!
			            	match &self.expr_index_obj.get_strict( &seqrec2.seq(),  &mut tool ){
			            		Some(gene_id) =>{
			            			report.iter_read_type( "expression reads" );
			            			if gene_id.len() == 1 {
			            				let data = GeneUmiHash( gene_id[0], *umi);
				            			self.gex.try_insert( 
				            				&(*cell_id as u64),
				            				data,
				            				report
				            			);
				            		}
				            		else {
				            			panic!("Multimapper? - I have removed this option!");
			            				/*self.gex.try_insert_multimapper(
				                    		&(*cell_id as u64),
				                    		gene_id,
				                    		umi,
				                    		report
				                    	);*/
				            		}
			                        //println!("R2 {}",String::from_utf8_lossy(&seqrec.id()).to_owned() );
			                        //println!("R1 {}",String::from_utf8_lossy(&seqrec1.id()).to_owned() );
			                    },
			                    None => {
			                    	report.no_data +=1;
			                    }
			                };
			            }
			        },
			        Err(_err) => {
			        	/*println!("R1 di not match to any cell: {} {} {}",
			        		String::from_utf8_lossy(&seqrec1.seq()[0..9]).to_owned(),
			        		String::from_utf8_lossy(&seqrec1.seq()[21..31]).to_owned(),
			        		String::from_utf8_lossy(&seqrec1.seq()[43..51]).to_owned() );*/
			        	report.no_sample +=1;
                    }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                };
                if report.ok_reads == report.max_reads{
                	break 'main;
                }
                report.log(&pb);

            } else {
            	println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        let log_str = report.log_str();

        pb.finish_with_message( log_str );

    }
*/
    /*
    c1 = (0, 9);
    c2 = (21,30);
    c3 = (43,52);
    */

    pub fn write_data( &mut self, outpath:String, results:&mut MappingInfo, min_umi : usize ) {
		    // calculating a little bit wrong - why? no idea...
	    //println!( "collected sample info:i {}", gex.mtx_counts( &mut expr_index_obj, &gene_names, opts.min_umi ) );

	    //let fp1 = PathBuf::from(opts.reads.clone());
	    //println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
	    let file_path = PathBuf::from(&outpath).join(
	    	"SampleCounts.tsv"
	    	);

	    let file_path_sp = PathBuf::from(&outpath).join(
	    	"BD_Rhapsody_expression"
	    	);

	    // this always first as this will decide which cells are OK ones!
	    results.stop_file_io_time();
	    
	    println!("filtering cells");
	    self.gex.mtx_counts( &mut self.te_index_obj, &self.te_names, min_umi, self.gex.num_threads ) ;
	    
	    results.stop_multi_processor_time();
	    println!("writing expression sets:");

	    if self.gene_names.len() > 0{ 
		    println!("writing gene expression");

		    match self.gex.write_sparse_sub ( file_path_sp, &mut self.expr_index_obj , &self.gene_names, min_umi ) {
		    	Ok(_) => (),
		    	Err(err) => panic!("Error in the data write: {err}")
		    };
		}

	    let file_path_sp = PathBuf::from(&outpath).join(
	    	"BD_Rhapsody_te_expression"
	    	);

	    println!("Writing TE counts");
	    match self.gex.write_sparse_sub ( file_path_sp, &mut self.te_index_obj, &self.te_names, 0 ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}")
	    };

	    println!("Writing samples table");

	    match self.gex.write_sub ( file_path, &mut self.samples, &self.sample_names, 0 ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}" )
	    };

	    
	    let reads_expr_index_obj = self.gex.n_reads( &mut self.expr_index_obj , &self.gene_names );
	    let reads_ab = self.gex.n_reads( &mut self.te_index_obj , &self.te_names );
	    let reads_samples = self.gex.n_reads( &mut self.samples , &self.sample_names );

	    println!( "{}",results.summary( reads_expr_index_obj, reads_ab, reads_samples) );

	    let file_path2 = format!("{}/SampleCounts.tsv", outpath );
	    println!( "\nCell->Sample table written to {file_path2:?}\n" );
	    results.stop_file_io_time();
	}
	
}


/// taken from https://github.com/onecodex/needletail/blob/eea49d01e5895a28ea85ba23a5e1a20a2fd7b378/src/parser/record.rs
/// as the library defines this part as private? Or at least I was not able to load it from there...  :-(
pub fn write_fastq(
	id: &[u8],
	seq: &[u8],
	qual: Option<&[u8]>,
	writer: &mut dyn Write,
	) -> Result<(), ParseError> {
	let ending = b"\n".to_owned();
	writer.write_all(b"@")?;
	writer.write_all(id)?;
	writer.write_all(&ending)?;
	writer.write_all(seq)?;
	writer.write_all(&ending)?;
	writer.write_all(b"+")?;
	writer.write_all(&ending)?;
    // this is kind of a hack, but we want to allow writing out sequences
    // that don't have qualitys so this will mask to "good" if the quality
    // slice is empty
    if let Some(qual) = qual {
    	writer.write_all(qual)?;
    } else {
    	writer.write_all(&vec![b'I'; seq.len()])?;
    }
    writer.write_all(&ending)?;
    Ok(())
}

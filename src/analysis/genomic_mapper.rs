use crate::analysis::MinimalSam;

use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
//use std::collections::HashSet;
use crate::singlecelldata::SingleCellData;
use crate::singlecelldata::cell_data::GeneUmiHash;
//use crate::geneids::GeneIds;
use crate::genes_mapper::GenesMapper;
//#[cfg(debug_assertions)]
//use crate::genes_mapper::CigarEndFix;
use crate::genes_mapper::{ SeqRec, NeedlemanWunschAffine};

//use crate::traits::BinaryMatcher;
// to access the command that was used to run this!
use std::env;
//use cargo_metadata::MetadataCommand;

use crate::errors::MappingError;

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
//use std::path::Path;
use std::io::{BufWriter, Write};
//use std::fs;

use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;

use crate::ofiles::{Ofiles, Fspot};
//use std::fs::write;
use glob::glob;

const VERSION: &str = env!("CARGO_PKG_VERSION");

//static EMPTY_VEC: Vec<String> = Vec::new();

#[derive(Debug)]
enum FilterError {
    Length,
    //PolyA,
    Ns,
    Quality,
}

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

/// the AnalysisGeneMapper cals is a wrapper around my 'old' quantify_rhapsody main funtion.
/// I started it to easily process multiple fastq files in a row.
pub struct AnalysisGenomicMapper{
	genes: GenesMapper,
	samples: GenesMapper,
	cells: Box<dyn CellIndex + Sync>,
	gex: SingleCellData,
	sample_names:Vec<String>,
	gene_names:Vec<String>,
	num_threads:usize,
}


impl AnalysisGenomicMapper{


	pub fn new(_gene_kmers:usize, version:String, specie: String, index:Option<String>, num_threads:usize, exp:&str, _debug: bool  ) -> Self{
		//let sub_len = 9;
	    //let mut cells = SampleIds::new( sub_len );// = Vec::with_capacity(12);
	    //cells.init_rhapsody( &opts.specie );

	    // let mut cell_umi:HashSet<u128> = HashSet::new();
	    //let mut genes :GeneIds = GeneIds::new(gene_kmers); // split them into 9 bp kmers
	    let mut genes :GenesMapper = GenesMapper::new( 0 ); // split them into 9 bp kmers
	    /*if debug {
	    	genes.debug(Some(debug));
	    }*/
	    //let mut gene_count = 0;
	    
	    if let Some(i) = index {
	    	println!("Loading index from path {i}");
	    	genes = match  GenesMapper::load_index( &i ){
	    		Ok(r) => r,
	    		Err(e) => panic!("Failed to load the index {e:?}")
	    	};
	    	genes.print();
	    	//gene_count = genes.names.len();
	    	
	    }

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

	    //panic!("Antibody count {} and expression count {}", antibodies.get_gene_count(), genes.get_gene_count());
	    let mut samples= GenesMapper::new( genes.get_gene_count() );
	    //let mut sample_names:Vec<String> = Vec::with_capacity(12);
	    /*if debug {
	    	samples.debug(Some(debug));
	    }*/
	    let mut id = 1;
	    if  specie.eq("human") {
	        // get all the human sample IDs into this.
	        // GTTGTCAAGATGCTACCGTTCAGAGATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG
	        //                        b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG"
	        // Seams they have introduced new sample ids, too....
	        let sequences = [ 
	        b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG",
	        b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", 
	        b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC",
	        b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT",
	        b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC",
	        b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG" ];

	        for seq in sequences{
	        	//seq.reverse();
	        	//let mut seq_ext = b"GTTGTCAAGATGCTACCGTTCAGAG".to_vec();
	        	//seq_ext.extend_from_slice( seq );
	        	samples.add( &seq.to_vec(), &format!("SampleTag{id:02}_hs"), &format!("SampleTag{id:02}_hs"), "human_sample", 0 );
	        	//sample_names.push( format!("Sample{id}") );
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

	        /* real world examples
	        "GTTGTCAAGATGCTACCG TTCAGAGTTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCTAAAA"
	        "                         ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC"
	        "GTTGTCAAGATGCTACCG TTCAGGACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCCAAAAA
	        "                          TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT"
			"GTTGTCAAGATGCTACCG TTCAGAGTTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT"
			"                          AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC"
			"GTTGTCAAGATGCTACG  TTCAGAGAGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTCAAGAA"
			*/
	        for seq in sequences{
	        	//seq.reverse();
	        	//let mut seq_ext = b"GTTGTCAAGATGCTACCGTTCAGAG".to_vec();
	        	//seq_ext.extend_from_slice( seq );
	        	//samples.add_small( &seq_ext, format!("Sample{id}"),EMPTY_VEC.clone() );
	        	samples.add( &seq.to_vec(), &format!("SampleTag{id:02}_mm"), &format!("SampleTag{id:02}_mm"), "mouse_sample", 0 );
	        	//sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }

	    } else {
	        println!("Sorry, but I have no primers for species {}", specie);
	        std::process::exit(1)
	    }

	    //samples.set_min_matches(2);
	    samples.set_highest_nw_val(0.2);
	    samples.set_min_matches( 5 );
	    genes.set_min_matches( 8 );


	    //samples.set_highest_humming_val(0.4);

	    println!("After indexing all fastq files we have the following indices:");
		println!("the mRNA index:");
		genes.print();
		println!("the sample id index:");
		samples.print();

		//this is problematic as it does not work with a genomic index! 
		//samples.make_index_te_ready();
		//genes.make_index_te_ready();
		//antibodies.make_index_te_ready();
		let sample_names: Vec<String>  = samples.get_all_gene_names();
		let gene_names: Vec<String>  = genes.get_all_gene_names();
	    genes.report4( &gene_names.iter().map(|s| s.as_str()).collect::<Vec<&str>>() );
	    samples.set_small_entries();

		Self{
			genes,
			samples,
	//		genes2,
			cells,
			gex,
			sample_names,
			gene_names,
			num_threads,
		}
	}

	pub fn report4gname( &mut self, gname: &[&str] ){
		//eprintln!("report4gname - Not supported at the moment");
		self.genes.report4(gname);
		self.samples.report4(gname);

	}

	// Function to set min_matches
    pub fn set_min_matches(&mut self, _value: usize) {
    	panic!("set_min_matches - Not supported at the moment");
		//self.genes.set_min_matches(value);
    }

    // Function to set highest_nw_val
    pub fn set_highest_nw_val(&mut self, value: f32) {
		self.genes.set_highest_nw_val(value);
    }

    // Function to set highest_humming_val
    pub fn set_highest_humming_val(&mut self, value: f32) {
    	//eprintln!("set_highest_humming_val - Not supported at the moment");
		self.genes.set_highest_humming_val(value);
    }

    pub fn debug( &mut self, value:Option<bool> ) {
		self.genes.debug(value);
		self.samples.debug(value);
    }

	pub fn write_index(&mut self, path:&str ){
		self.genes.write_index( path ).unwrap();
		self.genes.write_index_txt( path ).unwrap();
	}


    fn quality_control( read1:&SequenceRecord, read2:&SequenceRecord, min_quality:f32, pos: &[usize;8], min_sizes: &[usize;2]  ) -> Result<(), FilterError> {
        
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

        Ok(())
    }


    pub fn analyze_paralel( &self, data:&[(SeqRec, SeqRec)], report:&mut MappingInfo, _pos: &[usize;8] ) -> (SingleCellData, Vec<String>){
    	

        // first match the cell id - if that does not work the read is unusable
        //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
        let mut gex = SingleCellData::new( self.num_threads );
        let mut ok : bool;
        let minimal_sam = MinimalSam::new();

        let mut nwa = NeedlemanWunschAffine::new();

        //let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

        // lets tag this with the first gene I was interested in: Cd3e
        //let goi_id = self.genes.get_id("ADA".to_string());

        let mut bam = Vec::<String>::with_capacity(10_000);
        for i in 0..data.len() {
        	
        	match &self.cells.to_cellid( &data[i].0 ){
        		//   u32      u64
	            Ok( (cell_id, umi, cell_seq, umi_seq) ) => {
	            	//println!("cell {cell_id}: {}", String::from_utf8_lossy(&data[i].0).into_owned() );
	            	//let tool = IntToStr::new( data[i].0[(pos[6]+add)..(pos[7]+add)].to_vec(), 32 );
	            	report.cellular_reads +=1;

	            	// now I have three possibilites here:
	            	// an antibody tag match
	            	// a sample id match
	            	// or a mRNA match
	            	// And of casue not a match at all

	                ok = match &self.samples.get_strict( &data[i].1.seq().to_vec(), *cell_id as u32, &mut nwa){
	                    Ok(gene_id) =>{
	                    	report.iter_read_type( "sample reads" );
	                    	
                    		let guh = GeneUmiHash( gene_id[0].gene_id(), *umi as u64);
	                        if ! gex.try_insert( 
	                        	&(*cell_id as u64),
	                        	guh,
	                        	report
	                        ) { 
	                        	report.pcr_duplicates += 1 
	                        }
	                        // not interested in that!!!
	                        //bam.push(self.build_sam_record( &data[i].1, gene_id, cell_seq, umi_seq ));
	                        true
	                    },
	                    Err(MappingError::NoMatch) => {
	                    	false
	                    },
	                    Err(MappingError::MultiMatch) => {
	                    	// this is likely not mapping to anyting else if we alredy have mult matches here!
	                    	report.multimapper +=1;
	                    	report.no_data +=1;
	                    	false
	                    },
	                    Err(MappingError::OnlyCrap) => {
	                    	report.no_data +=1;
	                    	continue
	                    },
	                };
	                

	                if ! ok{
	                	
		                match &self.genes.get_strict( &data[i].1.seq().to_vec(), *cell_id as u32, &mut nwa){
		                	Ok(gene_id) =>{
		                		report.iter_read_type( "expression reads" );

		                    	let guh = GeneUmiHash( gene_id[0].gene_id(), *umi as u64);
		                        if ! gex.try_insert( 
		                        	&(*cell_id as u64),
		                        	guh,
		                        	report
		                        ){
		                        	report.pcr_duplicates += 1 
		                        }
		                        match minimal_sam.to_sam_line( &data[i].1, gene_id, cell_seq, umi_seq, &self.genes ) {
		                        	Some(sam_line) => bam.push( sam_line ),
		                        	None => {
		                        		eprintln!("There has been an error in the build_sam_record() function - please check what went wrong with this sequence data:\n{}\nand this the cell sequence:\n{}\n",&data[i].1, &data[i].0 );
		                        	}
		                        }
		                    },
		                    Err(MappingError::NoMatch) => {
		                    	// I want to be able to check why this did not work
		                    	/*report.write_to_ofile( Fspot::Buff1, 
		                    		format!(">Cell{cell_id} no gene detected\n{}\n", &data[i].1) 
		                    	);
								
								report.write_to_ofile( Fspot::Buff2, 
									format!(">Cell{cell_id} no gene detected\n{}\n", &data[i].0 ) 
								);*/
		                    	report.no_data +=1;
		                    },
		                    Err(MappingError::MultiMatch) => {
		                    	// this is likely not mapping to anyting else if we alredy have mult matches here!
		                    	report.multimapper +=1;
		                    	report.no_data +=1;
		                    },
		                    Err(MappingError::OnlyCrap) => {
		                    	report.no_data +=1;
		                    	continue
		                    },
		                };
		            }
	            },
	            Err(_err) => {
	            	// this is fucked up - the ids are changed!
	            	/*
	            	report.write_to_ofile( Fspot::Buff1, 
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
					*/
	                report.no_sample +=1;
	                continue
	            }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
       		};
    	}
        (gex, bam)
    }


    /// Analze BPO Rhapsody data in a paralel way.
    pub fn parse_parallel(&mut self,  f1:&str, f2:&str,  
    	report:&mut MappingInfo,pos: &[usize;8], min_sizes: &[usize;2], outpath: &str, max_reads:usize, chunk_size:usize ){

    	println!("I am using {} cpus", self.num_threads);

    	report.start_counter();
		let pattern = format!("{}/MappedReads*.sam", outpath);
    	let num_existing_files = match glob(&pattern) {
	        Ok(files) => files.count(),
	        Err(_err) => {
	            0
	        }
	    };

		let filename = format!("{}/MappedReads_{}.sam", outpath, num_existing_files +1 );

    	let file = match File::create(filename.clone()){
            Ok(f1) => f1,
            Err(err) => panic!("The file {filename} cound not be created: {err}" )
        };

    	let mut writer = BufWriter::new(file);
    	let mut header = "@HD\tVN:1.4\tSO:coordinate\n".to_string();
    	let mut fasta= "".to_string();
    	if let Some((h,f)) = self.genes.sam_header(){
    		header += &h;
    		fasta += &f;
    	}
    	if ! fasta.is_empty(){

			let fasta_f = match File::create(outpath.to_string() +"/GenesOfInterest.fasta"){
	            Ok(f1) => f1,
	            Err(err) => panic!("The file {} cound not be created: {err}", outpath.to_string() +"/GenesOfInterest.fasta" )
	        };
	        let mut fasta_w = BufWriter::new(fasta_f);
	        match write!(fasta_w, "{}", fasta) {
	    		Ok(_) => (),
	    		Err(err) => panic!("Could not write the fasta entry: {err:?}"),
	    	};	
	        
	    }
    	// this is just a copy from a real Illumina bam file - I also use the Sample4 as sample in my bam export.
    	// So if that is changed this header part needs to also change!
    	header += "@RG\tID:Sample4:0:1:HN2CKBGX9:1\tSM:Sample4\tLB:0.1\tPU:Sample4:0:1:HN2CKBGX9:1\tPL:ILLUMINA\n";
    	let args: Vec<String> = env::args().collect();
    	let program = args[0].to_string();
		let command_line: String = args.join(" ");

		header += &format!("@PG\tPN:{}\tID:{}\tVN:{}\tCL:{}",&program, &program, &VERSION, command_line);

		match writeln!(writer, "{}", header){
    		Ok(_) => (),
    		Err(err) => panic!("Could not write the header line: {err:?}"),
    	};

    	let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = match parse_fastx_file(f1) {
        	Ok(reader) => reader,
        	Err(err) => {
            	panic!("File 'reads' {f1} Error: {}", err);
        	}
   		};
        let mut readefile = match parse_fastx_file(f2) {
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
        let mut good_reads: Vec<(SeqRec, SeqRec)> = Vec::with_capacity( chunk_size * self.num_threads );
        let mut good_read_count = 0;

        'main: while let (Some(record1), Some(record2)) = (&readereads.next(), &readefile.next())  {
        	if report.total > max_reads{
        		break 'main
        	}
        	if good_read_count < chunk_size *self.num_threads {
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

        		match Self::quality_control( read1, read2, report.min_quality, pos, min_sizes){
        			Ok(()) => {
        				good_read_count +=1;
        				let r1 = SeqRec::new( read1.id(), &read1.seq(), read1.qual().unwrap() );
        				let r2 = SeqRec::new( read2.id(), &read2.seq(), read2.qual().unwrap() );
        				good_reads.push( (r1, r2 ) );
        				//good_reads.push( (read1.seq().to_vec(), read2.seq().to_vec() ) );
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
		    	let total_results: Vec<((SingleCellData,Vec<String>), MappingInfo)> = good_reads
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
			    	    let res = self.analyze_paralel(data_split, &mut rep, pos );
			    	    (res, rep)

		    	    }) // Analyze each chunk in parallel
		        .collect(); // Collect the results into a Vec
		        //eprintln!("sum up temp results");

		        report.stop_multi_processor_time();

			    for gex in total_results{
			    	self.gex.merge(gex.0.0);
			    	for line in gex.0.1{
			    		match writeln!(writer, "{}", line){
			        		Ok(_) => (),
			        		Err(err) => panic!("parse_parallel could not write the bam line: {err:?}"),
			        	}
			    		;
			    	}
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
	    	let total_results: Vec<((SingleCellData,Vec<String>), MappingInfo)> = good_reads
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
	        	self.gex.merge(gex.0.0);
	        	for line in gex.0.1{
		    		match writeln!(writer, "{}", line){
		        		Ok(_) => (),
		        		Err(err) => panic!("could not write to sam file? {err:?}"),
		        	}
		    	}
	        	report.merge( &gex.1 );
	        }
	        pb.set_message( report.log_str().clone() );
	    }

	    report.stop_single_processor_time();

	    let log_str = report.log_str();

        pb.finish_with_message( log_str );

    }


    /*
    c1 = (0, 9);
    c2 = (21,30);
    c3 = (43,52);
    */

    pub fn write_data( &mut self, outpath:String, results:&mut MappingInfo, min_umi : usize ) {
		    // calculating a little bit wrong - why? no idea...
	    //println!( "collected sample info:i {}", gex.mtx_counts( &mut genes, &gene_names, opts.min_umi ) );

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

	    let genes_idx = &self.genes.as_indexed_genes();
	    let samples_idx = &self.samples.as_indexed_genes();
	    
	    println!("filtering cells");
	    self.gex.update_genes_to_print( genes_idx, &self.gene_names );
	    self.gex.mtx_counts( genes_idx, min_umi, self.gex.num_threads ) ;
	    
	    results.stop_multi_processor_time();
	    println!("writing gene expression");

	    match self.gex.write_sparse_sub ( file_path_sp, genes_idx , &self.gene_names, min_umi ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}")
	    };

	    println!("Writing samples table");

	    match self.gex.write_sub ( file_path, samples_idx, &self.sample_names, 0 ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}" )
	    };

	    
	    let reads_genes = self.gex.n_reads( genes_idx , &self.gene_names );
	    let reads_samples = self.gex.n_reads( samples_idx , &self.sample_names );

	    println!( "{}",results.summary( reads_genes, 0,  reads_samples) );

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

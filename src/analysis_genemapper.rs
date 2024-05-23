use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
//use std::collections::HashSet;
use crate::singlecelldata::SingleCellData;
use crate::singlecelldata::cell_data::GeneUmiHash;
//use crate::geneids::GeneIds;
use crate::genes_mapper::GenesMapper;
use crate::genes_mapper::{ MapperResult, SeqRec, CigarEndFix};
//use crate::traits::BinaryMatcher;
// to access the command that was used to run this!
use std::env;

//use crate::int_to_str::IntToStr;
use crate::errors::MappingError;

use crate::cellids::CellIds;
use crate::cellids10x::CellIds10x;
use crate::traits::CellIndex;
 use crate::genes_mapper::NeedlemanWunschAffine;

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
use std::io::{ BufWriter, Write};

use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;

use crate::ofiles::{Ofiles, Fspot};
//use std::fs::write;

//static EMPTY_VEC: Vec<String> = Vec::new();
const VERSION: &str = env!("CARGO_PKG_VERSION");


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
pub struct AnalysisGeneMapper{
	genes: GenesMapper,
	samples: GenesMapper,
	antibodies: GenesMapper,
	cells: Box<dyn CellIndex + Sync>,
	gex: SingleCellData,
	sample_names:Vec<String>,
	gene_names:Vec<String>,
	ab_names:Vec<String>,
	num_threads:usize,
}


impl AnalysisGeneMapper{


	pub fn new(_gene_kmers:usize, version:String, expression:Option<String>, 
		antibody:Option<String>, specie:String, index:Option<String>, num_threads:usize, exp:&str, _debug:bool  ) -> Self{
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
	    	genes = match  GenesMapper::load_index( i ){
	    		Ok(r) => r,
	    		Err(e) => panic!("Failed to load the index {e:?}")
	    	};
	    	genes.print();
	    	//gene_count = genes.names.len();
	    	
	    }

	    if let Some(ex) = expression {
	    	if Path::new(&ex).exists(){

		    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

		    	while let Some(e_record) = expr_file.next() {
			        let seqrec = e_record.expect("invalid record");
		        	match std::str::from_utf8(seqrec.id()){
			            Ok(st) => {
		                	if let Some(id) = st.to_string().split('|').next(){
			                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
		                	}
		            	},
		            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
		        	}
		        }
		    }else {
		    	eprintln!("Expression file could not be read - ignoring")
		    }
	    }


	    eprintln!("Changing the expression start gene id to {}", genes.get_gene_count() );
	    let mut antibodies :GenesMapper = GenesMapper::new( genes.get_gene_count()  );

	    /*if debug {
	    	antibodies.debug(Some(debug));
	    }*/

	    if let Some(ab) = antibody {

		    if Path::new(&ab).exists(){


		   		let mut ab_file = parse_fastx_file(ab).expect("valid path/file");
		    	while let Some(ab_record) = ab_file.next() {
			        let seqrec = ab_record.expect("invalid record");
		        	match std::str::from_utf8(seqrec.id()){
			            Ok(st) => {
		                	if let Some(id) = st.to_string().split('|').next(){
		                		//seq_temp = seqrec.seq().to_vec();
		                		//seq_temp.reverse();
			                    antibodies.add( &seqrec.seq().to_vec(), id.to_string(), "antibody_tag".to_string(), 0 );
		                    	//ab_names.push( id.to_string() );
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
	    let mut samples= GenesMapper::new( antibodies.get_gene_count() + genes.get_gene_count() );
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
	        	samples.add( &seq.to_vec(), format!("SampleTag{id:02}_hs"), "human_sample".to_string(), 0 );
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
	        	samples.add( &seq.to_vec(), format!("SampleTag{id:02}_mm"), "mouse_sample".to_string(), 0 );
	        	//sample_names.push( format!("Sample{id}") );
	        	id +=1;
	        }

	    } else {
	        println!("Sorry, but I have no primers for species {}", specie);
	        std::process::exit(1)
	    }

	    //samples.set_min_matches(2);
	    samples.set_highest_nw_val(0.1);
	    genes.set_min_matches( 8 );
	    samples.set_min_matches( 5 );


	    //samples.set_highest_humming_val(0.4);

	    println!("After indexing all fastq files we have the following indices:");
		println!("the mRNA index:");
		genes.print();
		println!("the sample id index:");
		samples.print();
		println!("and the antibodies index:");
		antibodies.print();
		//this is problematic as it does not work with a genomic index! 
		//samples.make_index_te_ready();
		//genes.make_index_te_ready();
		//antibodies.make_index_te_ready();
		let sample_names: Vec<String>  = samples.get_all_gene_names();
		let gene_names: Vec<String>  = genes.get_all_gene_names();
		let ab_names: Vec<String>  = antibodies.get_all_gene_names();
		genes.report4( &gene_names.iter().map(|s| s.as_str()).collect::<Vec<&str>>() );
	    
		Self{
			genes,
			samples,
			antibodies,
	//		genes2,
			cells,
			gex,
			sample_names,
			gene_names,
			ab_names,
			num_threads,
		}
	}

	pub fn report4gname( &mut self, gname: &[&str] ){
		//eprintln!("report4gname - Not supported at the moment");
		self.antibodies.report4(gname);
		self.genes.report4(gname);
		self.samples.report4(gname);

	}

	// Function to set min_matches
    pub fn set_min_matches(&mut self, value: usize) {
    	//eprintln!("set_min_matches - Not supported at the moment");
    	self.antibodies.set_min_matches(value);
		//self.genes.set_min_matches(value);
    }

    // Function to set highest_nw_val
    pub fn set_highest_nw_val(&mut self, value: f32) {
        self.antibodies.set_highest_nw_val(value);
		self.genes.set_highest_nw_val(value);
    }

    // Function to set highest_humming_val
    pub fn set_highest_humming_val(&mut self, value: f32) {
    	//eprintln!("set_highest_humming_val - Not supported at the moment");
        self.antibodies.set_highest_humming_val(value);
		self.genes.set_highest_humming_val(value);
    }

    pub fn debug( &mut self, value:Option<bool> ) {
    	self.antibodies.debug(value);
		self.genes.debug(value);
		self.samples.debug(value);
    }

	pub fn write_index(&mut self, path:&String ){
		self.genes.write_index( path.to_string() ).unwrap();
		self.genes.write_index_txt( path.to_string() ).unwrap();
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

    fn build_sam_record ( &self, read:&SeqRec, gene_id:&Vec<MapperResult>, cell_id:&SeqRec, umi:&SeqRec ) -> Option<String>{

    	let mut read2 = read.clone();

    	match gene_id[0].cigar() {
    		Some(cigar) => 
    		{
    			let (mine, _other) = cigar.calculate_covered_nucleotides( &cigar.to_string() );
    			if mine < read2.len(){
    				read2=read2.slice(0, mine ).unwrap();

    			}
    			if mine + gene_id[0].start() > gene_id[0].db_length() {
    				panic!("Cigar suggest longer match than db_length allows!");
    			}
    		},
    		None=> {
    		},
    	};

    	let mut record = "".to_string();
    	// starting
    	record += &String::from_utf8_lossy(read2.id());
    	record += "\t";
    	// the ID
    	record += &format!("{}\t", 0b00000000000);
    	/*
            0x1 (1): PAIRED - Read is paired in sequencing.
		    0x2 (2): PROPER_PAIR - Read is mapped in a proper pair.
		    0x4 (4): UNMAPPED - Read is unmapped.
		    0x8 (8): MATE_UNMAPPED - Mate is unmapped.
		    0x10 (16): REVERSE_STRAND - Read is mapped to the reverse strand.
		    0x20 (32): MATE_REVERSE_STRAND - Mate is mapped to the reverse strand.
		    0x40 (64): READ1 - Read is the first read in a pair.
		    0x80 (128): READ2 - Read is the second read in a pair.
		    0x100 (256): SECONDARY - Secondary alignment.
		    0x200 (512): QC_FAIL - Read fails quality checks.
		    0x400 (1024): DUPLICATE - PCR or optical duplicate.
		 */
		record += &format!("{}\t", gene_id[0].get_name() );
		// the name of the database sequence
    	record += &format!("{}\t", gene_id[0].start() +1);
    	// the start position of the read
    	record += &format!("{}\t", gene_id[0].mapq() );
    	// the map quality
    	let cig = gene_id[0].get_cigar();
    	let (mine, _other) = cig.calculate_covered_nucleotides( &cig.to_string() );
    	if mine != &read_mod.len(){
    		panic!("I am trying to create a bam line and found a discrepancy between cigar length and sequence length: {cig}\n{read2}");
    	}
		let read_mod = match cig.fixed {
			Some(CigarEndFix::StartInsert) => {
				//shit - the read needs to be adjusted!
				let (cig_mapping_length, _) = cig.calculate_covered_nucleotides( &cig.cigar );
				let bp_to_clip = read2.len() - cig_mapping_length;
				read2.slice( bp_to_clip, cig_mapping_length ).unwrap()
			},
			_ => {
				read2.clone()
			}
		};
    	let cigar_string = match gene_id[0].cigar() {
		    Some(cigar) => cigar.to_string(),
		    None => {
		    	format!("{}S", read_mod.seq().len())
		    }
		};
		record += &cigar_string ;
    	record +="\t";
		// the Cigar string
		record += "*\t";
		// the sequence name of the mate/next read
		record += "0\t";
		//1-based position of the mate/next read
		record += "0\t";
		//Template length
		record += &String::from_utf8_lossy( &read_mod.seq() ) ;
    	record +="\t";
		//Read sequence
		record += &String::from_utf8_lossy( &read_mod.qual() ) ;
    	record +="\t";
		//Phred quality scores

		// and here come the 10x specific tags - might be needed later on
		record += &format!("NH:i:{}\t", gene_id.len());
	    /*NH:i:9: Number of reported alignments.
	        NH: Tag identifier, stands for "Number of reported alignments."
	        i: Tag data type, denotes a 32-bit signed integer.
	        9: Value associated with the tag, indicating the number of reported alignments for the read.*/
	    record += "HI:i:1\t";
		/*HI:i:2: Alignment hit index.
	        HI: Tag identifier, stands for "Alignment hit index."
	        i: Tag data type, denotes a 32-bit signed integer.
	        2: Value associated with the tag, representing the alignment hit index for the read. This index is typically used to distinguish between multiple hits for the same query in the SAM/BAM file.
		*/
	    record += &format!("AS:i:{}\t", gene_id[0].mapq() ); // AS:i:89
	    /*AS:i:89: Alignment score.
	        AS: Tag identifier, stands for "Alignment score."
	        i: Tag data type, denotes a 32-bit signed integer.
	        89: Value associated with the tag, indicating the alignment score for the read. The alignment score represents the quality or confidence of the alignment.
		*/
	    record += &format!("nM:i:{}\t", gene_id[0].edit_dist() ); // nM:i:0
	    /*nM:i:0: Edit distance.
	        nM: Tag identifier, stands for "Edit distance."
	        i: Tag data type, denotes a 32-bit signed integer.
	        0: Value associated with the tag, representing the edit distance between the read and the reference sequence. The edit distance is the number of mismatches and gaps (insertions and deletions) between the read and the reference sequence.
		*/
	    record += "RE:A:I\t"; // RE:A:I - I do not report anything else
	    /*
	    RE:A:I: Read extension.
	        RE: Tag identifier, stands for "Read extension."
	        A: Tag data type, denotes a single character (ASCII).
	        I: Value associated with the tag, indicating the type of read extension. In this case, I represents an insert read extension.
	    Possible values:
	    1. Insert Read Extension (RE:A:I): Indicates that the read is an insert read, meaning it is derived from the original DNA fragment without any modifications or rearrangements.
	    2. Long Fragment Read Extension (RE:A:L): Indicates that the read is a long fragment read, suggesting that it spans a longer portion of the original DNA fragment compared to other reads.
	    3. Linked Read Extension (RE:A:Q): Indicates that the read is a linked read, meaning it is part of a set of reads that are linked together, typically by sharing a common barcode or unique molecular identifier (UMI).
	    4. Chimeric Read Extension (RE:A:C): Indicates that the read is a chimeric read, meaning it contains sequences from two or more distinct DNA molecules or genomic regions.
	    */
	    //record += &format!("xf:i:{}\t",gene_id[0].edit_dist()); // li:i:0 // not linked

	    record += "li:i:0\t"; // li:i:0 // not linked
	    /*li:i:0: Linked-read identifier.
	        li: Tag identifier, stands for "Linked-read identifier."
	        i: Tag data type, denotes a 32-bit signed integer.
	        0: Value associated with the tag, representing the linked-read identifier for the read. Linked-reads are reads that are part of the same molecule or DNA fragment.
		*/
	    record += &format!("BC:Z:{}\t", umi.as_dna_string() ); // BC:Z:GCCATTCC
	    /*BC:Z:GCCATTCC: Barcode sequence.
	        BC: Tag identifier, stands for "Barcode sequence."
	        Z: Tag data type, denotes a string.
	        GCCATTCC: Value associated with the tag, representing the barcode sequence assigned to the read.
		*/
	    record += &format!("QT:Z:{}\t", &String::from_utf8_lossy(umi.qual()) ); // QT:Z:AAAAAEEE
	    /*QT:Z:AAAAAEEE: Quality tags.
	        QT: Tag identifier, stands for "Quality tags."
	        Z: Tag data type, denotes a string.
	        AAAAAEEE: Value associated with the tag, representing the quality scores of the read bases.
		*/
		record += &format!("CR:Z:{}\t",cell_id.as_dna_string() ); // CR:Z:CTCCTTTCATACTGTG
	    /*CR:Z:CTCCTTTCATACTGTG: Chromium cellular barcode sequence.
	        CR: Tag identifier, stands for "Chromium cellular barcode sequence."
	        Z: Tag data type, denotes a string.
	        CTCCTTTCATACTGTG: Value associated with the tag, representing the cellular barcode sequence assigned by the Chromium platform.
		*/
	    record += &format!("CY:Z:{}\t",  String::from_utf8_lossy(cell_id.qual()) ); // CY:Z:AAAAAEEEEEEEEEEE
	    /*CY:Z:AAAAAEEEEEEEEEEE: Chromium cellular quality scores.
	        CY: Tag identifier, stands for "Chromium cellular quality scores."
	        Z: Tag data type, denotes a string.
	        AAAAAEEEEEEEEEEE: Value associated with the tag, representing the quality scores of the cellular barcode bases.
		*/
	    record += &format!("CB:Z:{}-1\t", cell_id.as_dna_string() ); // CB:Z:CTCCTTTCATACTGTG-1
	    /*CB:Z:CTCCTTTCATACTGTG-1: Chromium gem group barcode sequence.
	        CB: Tag identifier, stands for "Chromium gem group barcode sequence."
	        Z: Tag data type, denotes a string.
	        CTCCTTTCATACTGTG-1: Value associated with the tag, representing the gem group barcode sequence assigned by the Chromium platform.
		*/
	    record += &format!("UR:Z:{}\t", umi.as_dna_string() ); // UR:Z:TTTAGGTCTTGG
	    /*UR:Z:TTTAGGTCTTGG: Chromium UMI sequence.
	        UR: Tag identifier, stands for "Chromium UMI sequence."
	        Z: Tag data type, denotes a string.
	        TTTAGGTCTTGG: Value associated with the tag, representing the unique molecular identifier (UMI) sequence assigned by the Chromium platform.
		*/
	    record += &format!("UZ:Z:{}\t", String::from_utf8_lossy(umi.qual()) ); // UY:Z:EEEEEEEEEEEE
	    /*UY:Z:EEEEEEEEEEEE: Chromium UMI quality scores.
	        UY: Tag identifier, stands for "Chromium UMI quality scores."
	        Z: Tag data type, denotes a string.
	        EEEEEEEEEEEE: Value associated with the tag, representing the quality scores of the UMI bases.
		*/
	    record += &format!("UB:Z:{}\t", umi.as_dna_string()); // UB:Z:TTTAGGTCTTGG
	    /*UB:Z:TTTAGGTCTTGG: Chromium corrected UMI sequence.
	        UB: Tag identifier, stands for "Chromium corrected UMI sequence."
	        Z: Tag data type, denotes a string.

		*/	
	    record += &format!("RG:Z:{}", "Sample4:0:1:HN2CKBGX9:1"); // RG:Z:Sample4:0:1:HN2CKBGX9:1

	    Some(record)
    }

    pub fn analyze_paralel( &self, data:&[(SeqRec, SeqRec)], report:&mut MappingInfo, _pos: &[usize;8] ) -> (SingleCellData, Vec<String>){
    	

        // first match the cell id - if that does not work the read is unusable
        //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
        let mut gex = SingleCellData::new( self.num_threads );
        let mut ok : bool;
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


	            	ok = match &self.antibodies.get_strict( &data[i].1.seq().to_vec(), *cell_id as u32, &mut nwa ){
	                    Ok(gene_id) =>{

	                    	report.iter_read_type( "antibody reads" );
                    		
                    		let guh = GeneUmiHash( gene_id[0].gene_id(), *umi as u64);
	                    	if ! gex.try_insert( 
	                        	&(*cell_id as u64),
	                        	guh,
	                        	report
	                        ){
	                        	report.pcr_duplicates += 1 
	                        }


	                        if gene_id[0].save(){
	                        	//build_sam_record ( &self, read2:&SeqReq, gene_id:&Vec<MapperResult>, cell_id:&SeqReq, umi:&SeqReq ) 
							    //bam.push(self.build_sam_record( &data[i].1, gene_id, cell_seq, umi_seq ));
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
	                    	false
	                    }
	                };

	                if ! ok{
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

		                        if gene_id[0].save(){
		                        	//build_sam_record ( &self, read2:&SeqReq, gene_id:&Vec<MapperResult>, cell_id:&SeqReq, umi:&SeqReq ) 
								    //bam.push(self.build_sam_record( &data[i].1, gene_id, cell_seq, umi_seq ));
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

		                        if gene_id[0].save(){
		                        	//build_sam_record ( &self, read2:&SeqReq, gene_id:&Vec<MapperResult>, cell_id:&SeqReq, umi:&SeqReq )
		                        	match self.build_sam_record( &data[i].1, gene_id, cell_seq, umi_seq ) {
			                        	Some(sam_line) => bam.push( sam_line ),
			                        	None => {
			                        		eprintln!("There has been an error in the build_sam_record() function - please check what went wrong with this sequence:\n{}.",&data[i].1 );
			                        	}
			                        }
		                        }
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

    	let file = match File::create(outpath.to_string() +"/AlignedReadsOfInterest.sam"){
            Ok(f1) => f1,
            Err(err) => panic!("The file {} cound not be created: {err}", outpath.to_string() +"/AlignedReadsOfInterest.sam" )
        };

    	let mut sam_file_writer = BufWriter::new(file);
    	let mut header = "@HD\tVN:1.4\tSO:coordinate\n".to_string();
    	let mut fasta= "".to_string();
    	if let Some((h,f)) = self.samples.sam_header(){
    		header += &h;
    		fasta += &f;
    	}
    	if let Some((h,f)) = self.antibodies.sam_header(){
    		header += &h;
    		fasta += &f;
    	}
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
	        
	        match write!(fasta_w, "{}", fasta){
	    		Ok(_) => (),
	    		Err(err) => panic!("could not write the fasta entry? {err:?}"),
	    	};
	    }
    	// this is just a copy from a real Illumina bam file - I also use the Sample4 as sample in my bam export.
    	// So if that is changed this header part needs to also change!
    	header += "@RG\tID:Sample4:0:1:HN2CKBGX9:1\tSM:Sample4\tLB:0.1\tPU:Sample4:0:1:HN2CKBGX9:1\tPL:ILLUMINA\n";
    	let args: Vec<String> = env::args().collect();
    	let program = args[0].to_string();
		let command_line: String = args.join(" ");
		// Get metadata about the package
		header += &format!("@PG\tPN:{}\tID:{}\tVN:{}\tCL:{}",&program, &program, &VERSION, command_line);

    	match writeln!(sam_file_writer, "{}", header){
    		Ok(_) => (),
    		Err(err) => panic!("could not write the header line? {err:?}"),
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
			    	self.gex.merge(&gex.0.0);
			    	for line in gex.0.1{
			    		match writeln!(sam_file_writer, "{}", line){
			        		Ok(_) => (),
			        		Err(err) => panic!("could not write to sam file? {err:?}"),
			        	}
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
	        	self.gex.merge(&gex.0.0);
	        	for line in gex.0.1{
		    		match writeln!(sam_file_writer, "{}", line){
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
	    println!("Indexed the genes names");
	    let ab_idx = &self.antibodies.as_indexed_genes();
	    println!("Indexed the antibody names");
	    let samples_idx = &self.samples.as_indexed_genes();
	    println!("Indexed the samples names");
	    
	    println!("filtering cells");
	    self.gex.update_genes_to_print( genes_idx, &self.gene_names );
	    self.gex.mtx_counts( genes_idx, min_umi, self.gex.num_threads ) ;
	    
	    results.stop_multi_processor_time();
	    println!("writing gene expression");

	    match self.gex.write_sparse_sub ( file_path_sp, genes_idx , &self.gene_names, min_umi ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}")
	    };

	    let file_path_sp = PathBuf::from(&outpath).join(
	    	"BD_Rhapsody_antibodies"
	    	);

	    println!("Writing Antibody counts");
	    match self.gex.write_sparse_sub ( file_path_sp, ab_idx, &self.ab_names, 0 ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}")
	    };

	    println!("Writing samples table");

	    match self.gex.write_sub ( file_path, samples_idx, &self.sample_names, 0 ) {
	    	Ok(_) => (),
	    	Err(err) => panic!("Error in the data write: {err}" )
	    };

	    
	    let reads_genes = self.gex.n_reads( genes_idx , &self.gene_names );
	    let reads_ab = self.gex.n_reads( ab_idx , &self.ab_names );
	    let reads_samples = self.gex.n_reads( samples_idx , &self.sample_names );

	    println!( "{}",results.summary( reads_genes, reads_ab, reads_samples) );

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

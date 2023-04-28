use needletail::parse_fastx_file;
use kmers::naive_impl::Kmer;

use crate::cellids::CellIds;
use crate::singlecelldata::SingleCellData;
//use crate::geneids::GeneIds;
use crate::fast_mapper::FastMapper;

use crate::ofiles::Ofiles;
//use std::io::BufReader;
//use flate2::write::GzDecoder;

//use this::last5::Last5;

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
//use std::time::Duration;


use std::path::PathBuf;
use std::fs::File;
use std::path::Path;
use std::io::Write;


/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct MappingInfo{
	/// reads that did not pass the filters
	pub unknown:usize,
	/// reads that had no cell id
    pub no_sample:usize,
    /// reads that have no match in the geneIds object
    pub no_data:usize,
    /// reads with cell id and gene id
    pub ok_reads:usize,
    /// reads that are duplicates on the UMI level per cell and gene
    pub pcr_duplicates:usize,
    /// the amount of ok_reads after which to write a entry into the log file
   	pub split:usize,
   	/// the others are explained in the quantify_rhapsody.rs file.
    log_iter:usize,
    pub log_writer:File,
    pub min_quality:f32, 
    pub max_reads:usize, 
    pub ofile:Ofiles,
    pub local_dup:usize,
    pub total:usize,
}

impl MappingInfo{
	pub fn new(log_writer:File, min_quality:f32, max_reads:usize, ofile:Ofiles, ) -> Self{
		let unknown = 0;
		let no_sample = 0;
		let no_data = 0;
		let ok_reads = 0;
		let pcr_duplicates = 0;
		let split = 1_000_000;
		let log_iter = 0;
		let local_dup = 0;
		let total = 0;
		Self{
			unknown,
			no_sample,
			no_data,
			ok_reads,
			pcr_duplicates,
			split,
			log_iter,
			log_writer,
			min_quality,
			max_reads,
			ofile,
			local_dup,
			total,
		}
	}
	pub fn log( &mut self, pb:&ProgressBar ){
		if self.total % self.split == 0{
			self.log_iter+=1;
            let log_str = self.log_str();
            pb.set_message( log_str.clone() );
            pb.inc(1);
            match writeln!( self.log_writer, "{log_str}" ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}" );
                }
            };
            self.local_dup = 0;
            //std::thread::sleep(Duration::from_millis(100));
		}
	}
	pub fn log_str( &mut self ) -> String{
		format!("{:.2} mio reads ({:.2}% usable; {:.2}% PCR dupl. [usable] [{:.2}% for last batch])",
            self.total as f32 / self.split as f32,
            self.ok_reads as f32 / (self.ok_reads +self.no_sample+ self.unknown) as f32 * 100.0 , 
            self.pcr_duplicates as f32 / self.ok_reads as f32 * 100.0,
            self.local_dup as f32 / self.split as f32 * 100.0
         )
	}
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

/// the analysis calss is a wrapper around my 'old' quantify_rhapsody main funtion.
/// I started it to easily process multiple fastq files in a row.
pub struct Analysis<'a>{
	genes: FastMapper,
	cells: CellIds<'a>,
	gex: SingleCellData,
	sample_names:Vec<String>,
	gene_names:Vec<String>,
	ab_names:Vec<String>,
}


impl Analysis<'_>{


	pub fn new(gene_kmers:usize, version:String, expression:String, antibody:String, specie:String  ) -> Self{
		//let sub_len = 9;
	    //let mut cells = SampleIds::new( sub_len );// = Vec::with_capacity(12);
	    //cells.init_rhapsody( &opts.specie );

	    // let mut cell_umi:HashSet<u128> = HashSet::new();
	    //let mut genes :GeneIds = GeneIds::new(gene_kmers); // split them into 9 bp kmers
	    let mut genes :FastMapper = FastMapper::new( gene_kmers ); // split them into 9 bp kmers
	     
	    //  now we need to get a CellIDs object, too
	    let cells = CellIds::new( &version, 7);

	    // that is a class to strore gene expression data.
	    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
	    // and umi's are u64 again
	    // here I need the cell kmer site.
	    let gex = SingleCellData::new( );

	    

	    let mut sample_names:Vec<String> = Vec::with_capacity(12);
	    let mut gene_names:Vec<String> = Vec::with_capacity(600);
	    let mut ab_names:Vec<String> = Vec::with_capacity(30);

	    for i in 1..13{
	        sample_names.push( format!("Sample{i}") )
	    }

	    if  specie.eq("human") {
	        // get all the human sample IDs into this.
	        genes.add( &b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG".to_vec(), "Sample1".to_string() );
	        genes.add( &b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG".to_vec(), "Sample2".to_string() );
	        genes.add( &b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT".to_vec(), "Sample3".to_string() );
	        genes.add( &b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT".to_vec(), "Sample4".to_string() );
	        genes.add( &b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG".to_vec(), "Sample5".to_string() );
	        genes.add( &b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC".to_vec(), "Sample6".to_string() );
	        genes.add( &b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT".to_vec(), "Sample7".to_string() );
	        genes.add( &b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT".to_vec(), "Sample8".to_string() );
	        genes.add( &b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG".to_vec(), "Sample9".to_string() );
	        genes.add( &b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC".to_vec(), "Sample10".to_string() );
	        genes.add( &b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC".to_vec(), "Sample11".to_string() );
	        genes.add( &b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG".to_vec(), "Sample12".to_string() );
	    }
	    else if specie.eq("mouse") {
	        // and the mouse ones
	        genes.add( &b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG".to_vec(), "Sample1".to_string() );
	        genes.add( &b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC".to_vec(), "Sample2".to_string() );
	        genes.add( &b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC".to_vec(), "Sample3".to_string() );
	        genes.add( &b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT".to_vec(), "Sample4".to_string() );
	        genes.add( &b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT".to_vec(), "Sample5".to_string() );
	        genes.add( &b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC".to_vec(), "Sample6".to_string() );
	        genes.add( &b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC".to_vec(), "Sample7".to_string() );
	        genes.add( &b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG".to_vec(), "Sample8".to_string() );
	        genes.add( &b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC".to_vec(), "Sample9".to_string() );
	        genes.add( &b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC".to_vec(), "Sample10".to_string() );
	        genes.add( &b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT".to_vec(), "Sample11".to_string() );
	        genes.add( &b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC".to_vec(), "Sample12".to_string() );

	    } else {
	        println!("Sorry, but I have no primers for species {}", specie);
	        std::process::exit(1)
	    }
	    
	    if Path::new(&expression).exists(){

	    	let mut expr_file = parse_fastx_file(expression).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string() );
	                    	gene_names.push( id.to_string() );
	                    	//genes2.add_unchecked( &seqrec.seq(), id.to_string() );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	eprintln!("Expression file could not be read - ignoring")
	    }

	    if Path::new(&antibody).exists(){

	   		let mut ab_file = parse_fastx_file(antibody).expect("valid path/file");
	    	while let Some(ab_record) = ab_file.next() {
		        let seqrec = ab_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string() );
	                    	ab_names.push( id.to_string() );
	                    	//genes2.add_unchecked( &seqrec.seq(), id.to_string() );
	                	};
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	eprintln!("Antibody file could not be read - ignoring")
	    }
		Self{
			genes,
	//		genes2,
			cells,
			gex,
			sample_names,
			gene_names,
			ab_names,
		}
	}


	pub fn parse(&mut self,  f1:String, f2:String,  report:&mut MappingInfo,pos: &[usize;8], min_sizes: &[usize;2]  ){

        let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&f1).expect("valid path/file");
        let mut readefile = parse_fastx_file(&f2).expect("valid path/file");
        let m = MultiProgress::new();
        let pb = m.add(ProgressBar::new(5000));
        pb.set_style(spinner_style);
        //pb.set_prefix(format!("[{}/?]", i + 1));

        let report_gid = self.genes.get_id( "Cd3e".to_string() );

        'main: while let Some(record2) = readefile.next() {
            if let Some(record1) = readereads.next() {

            	report.total += 1;
                
                let seqrec = match record2{
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


                // the quality measurements I get here do not make sense. I assume I would need to add lots of more work here.
                //println!("The read1 mean quality == {}", mean_u8( seqrec.qual().unwrap() ));

                if mean_u8( seqrec.qual().unwrap() ) < report.min_quality {
                    report.unknown +=1;
                    //println!("filtered a read: {:?} ({})",std::str::from_utf8(&seqrec.seq()), mean_u8( seqrec.qual().unwrap()  ) );
                    continue 'main;
                }
                if mean_u8( seqrec1.qual().unwrap() ) < report.min_quality {
                    report.unknown +=1;
                    //println!("filtered a read: {:?} ({})",std::str::from_utf8(&seqrec1.seq()), mean_u8( seqrec1.qual().unwrap() ) );
                    continue 'main;
                }

                // totally unusable sequence
                // as described in the BD rhapsody Bioinformatic Handbook
                // GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf (google should find this)
                if seqrec1.seq().len() < min_sizes[0] {
                    report.unknown +=1;
                    continue;
                }
                if seqrec.seq().len() < min_sizes[1] {
                    report.unknown +=1;
                    continue;
                }
                //let seq = seqrec.seq().into_owned();
                for nuc in &seqrec1.seq()[pos[6]..pos[7]] {  
                    if *nuc ==b'N'{
                        report.unknown +=1;
                        continue 'main;
                    }
                }
                //let umi = Kmer::from( &seqrec1.seq()[52..60]).into_u64();
                let umi = Kmer::from( &seqrec1.seq()[pos[6]..pos[7]]).into_u64();

                // first match the cell id - if that does not work the read is unusable
                //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                match &self.cells.to_cellid( &seqrec1.seq(), vec![pos[0],pos[1]], vec![pos[2],pos[3]], vec![pos[4],pos[5]]){
                    Ok(cell_id) => {
                        // this is removing complexity from the data - in the test dataset 111 reads are ignored.
                        // let cell_id_umi:u128 = read_be_u128(  [ umi.to_be_bytes() , (cell_id as u64).to_be_bytes() ].concat().as_slice() );
                        // if ! cell_umi.insert( cell_id_umi ){
                        //     continue 'main;
                        // }
                        match &self.genes.get( &seqrec.seq() ){
                            Some(gene_id) =>{
                                self.gex.try_insert( 
                                	*cell_id as u64, 
                                	format!( "Cell{cell_id}",  ),
                                	gene_id,
                                	umi,
                                	report
                                );
                                if *gene_id == report_gid {
                                	match seqrec1.write(&mut report.ofile.buff1, None){
	                                    Ok(_) => (),
	                                    Err(err) => println!("{err}")
	                                };
	                                match seqrec.write( &mut report.ofile.buff2, None){
	                                    Ok(_) => (),
	                                    Err(err) => println!("{err}")
	                                };                              
                            		//println!("Cool I got a gene id: {gene_id}", );
                            	}
                            },
                            None => {
                                // match &self.genes2.get_unchecked( &seqrec.seq() ){
                                //     Some(gene_id) =>{ 
                                //         self.gex.try_insert( 
                                //         	*cell_id as u64, 
                                //         	format!( "Cell{cell_id}"  ),
                                //         	gene_id,
                                //         	umi, 
                                // 			report
                                //         	);
                                //         if *gene_id == report_gid {
		                                	match seqrec1.write(&mut report.ofile.buff1, None){
			                                    Ok(_) => (),
			                                    Err(err) => println!("{err}")
			                                };
			                                match seqrec.write( &mut report.ofile.buff2, None){
			                                    Ok(_) => (),
			                                    Err(err) => println!("{err}")
			                                };                              
		                        //     		println!("Cool I got a gene id: {gene_id}", );
		                        //     	}
	                            //     },
	                            //     None => {
	                            //     	// match seqrec1.write(&mut report.ofile.buff1, None){
		                        //         //     Ok(_) => (),
		                        //         //     Err(err) => println!("{err}")
		                        //         // };
		                        //         // match seqrec.write( &mut report.ofile.buff2, None){
		                        //         //     Ok(_) => (),
		                        //         //     Err(err) => println!("{err}")
		                        //         // };
	                            //      	report.no_data +=1;
	                            //     }
	                            // };
                                report.no_data +=1;

                                // all - samples genes and antibodies are classed as genes here.
                            }
                        };

                    },
                    Err(_err) => {
                        // cell id could not be recovered
                        
                        // println!("Cell ID could not be recovered from {:?}:\n{}\n{:?}, {:?}, {:?}", std::str::from_utf8(&seqrec1.seq()), _err, 
                        //      std::str::from_utf8( &seqrec1.seq()[pos[0]..pos[1]]),
                        //      std::str::from_utf8( &seqrec1.seq()[pos[2]..pos[3]]),
                        //      std::str::from_utf8( &seqrec1.seq()[pos[4]..pos[5]])
                        // );
                        
                        report.no_sample +=1;
                        continue
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

	pub fn write_data( &mut self, outpath:String, results:&MappingInfo, min_umi : usize ) {
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
    match self.gex.write_sparse_sub ( file_path_sp, &mut self.genes , &self.gene_names, min_umi ) {
        Ok(_) => (),
        Err(err) => panic!("Error in the data write: {err}")
    };

    let file_path_sp = PathBuf::from(&outpath).join(
        "BD_Rhapsody_antibodies"
    );

    match self.gex.write_sparse_sub ( file_path_sp, &mut self.genes, &self.ab_names, 1 ) {
        Ok(_) => (),
        Err(err) => panic!("Error in the data write: {err}")
    };


    match self.gex.write_sub ( file_path, &mut self.genes, &self.sample_names, 0 ) {
        Ok(_) => (),
        Err(err) => panic!("Error in the data write: {err}" )
    };

    

    let total = results.no_sample+ results.unknown + results.ok_reads;
    let reads_genes = self.gex.n_reads( &mut self.genes , &self.gene_names );
    let reads_ab = self.gex.n_reads( &mut self.genes , &self.ab_names );
    let reads_samples = self.gex.n_reads( &mut self.genes , &self.sample_names );

    println!( "\nSummary:");
    println!(     "total      reads  : {} reads", results.total );
    println!(     "no cell ID reads  : {} reads", results.no_sample );
    println!(     "no gene ID reads  : {} reads", results.no_data );
    println!(     "N's or too short  : {} reads", results.unknown );
    println!(     "cellular reads    : {} reads ({:.2}% of total)", results.ok_reads, (results.ok_reads as f32 / results.total as f32) * 100.0 );
    println!(     "expression reads  : {} reads ({:.2}% of total)", reads_genes, (reads_genes as f32 / results.total as f32) * 100.0 );
    println!(     "antibody reads    : {} reads ({:.2}% of total)", reads_ab, (reads_ab as f32 / total as f32) * 100.0 );
    println!(     "sample tag reads  : {} reads ({:.2}% of total)", reads_samples, (reads_samples as f32 / results.total as f32) * 100.0 );
    println!(     "pcr duplicates    : {} reads ({:.2}% of usable)", results.pcr_duplicates, ( results.pcr_duplicates as f32 / results.ok_reads as f32 ) * 100.0 );
    

    let file_path2 = format!("{}/SampleCounts.tsv", outpath );
    println!( "\nCell->Sample table written to {file_path2:?}\n" );
	}
}
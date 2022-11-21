use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
// use flate2::read::GzDecoder;

use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::io::BufWriter;
// use std::io::BufReader;
// use std::str;

use kmers::naive_impl::Kmer;

use this::cellids::CellIds;

use this::sampleids::SampleIds;
use this::cellids10x::CellIds10x;
use this::geneids::GeneIds;

use std::path::PathBuf;
use std::fs::File;
use std::fs;

use std::io::prelude::*;
use std::time::SystemTime;

// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl

/// Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    specie: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the fastq database containing the genes
    #[clap(short, long)]
    expression: String,
    /// the fastq database containing the antibody tags
    #[clap(short, long)]
    antybody: String,
}

struct Ofiles {
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufWriter<GzEncoder<File>>,
    pub buff2: BufWriter<GzEncoder<File>>
}

/// Ofiles encapsulates two BufWriter<GzEncoder<File>> objects to make handling of 20 of these more convenient.
impl Ofiles{
    pub fn new(id: usize, opts: &Opts )->Self {

        let fp1 = PathBuf::from(opts.reads.clone());
        let fp2 = PathBuf::from(opts.file.clone());
        let file1_path = PathBuf::from(&opts.outpath).join(format!("{}.{}", id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(&opts.outpath).join(format!("{}.{}", id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::create(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {} cound not be created: {}", format!("{}.{}", id, fp1.file_name().unwrap().to_str().unwrap() ), err)
        };
        let f2 = match File::create(file2_path){
            Ok(file) => file,
            Err(err) => panic!("The file {} cound not be created: {}", format!("{}.{}", id, fp2.file_name().unwrap().to_str().unwrap() ), err)
        };
        
        let file1 = GzEncoder::new(f1, Compression::default());
        let file2 = GzEncoder::new(f2, Compression::default());

        let buff1 = BufWriter::new( file1 );
        let buff2 = BufWriter::new( file2 );

        let count:u32 = 0;
        Self{
            count,
            //file1,
            //file2,
            buff1,
            buff2
        }
    }
    pub fn close( &mut self ){

        match self.buff1.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R1: {}",e),
        };
        match self.buff2.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R2: {}",e),
        };
    }
}


// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    match fs::create_dir_all(&opts.outpath){
        Ok(_) => (),
        Err(e) => panic!("I could not create the outpath: {}", e)
    };

    println!("starting to collect the expression data can take a LONG time");

    let sub_len = 9;
    let mut samples = SampleIds::new( sub_len );// = Vec::with_capacity(12);
    samples.init_rhapsody( &opts.specie );
 
    let mut cell2sample: Vec<HashSet<u32>>;
    cell2sample = vec![
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
        HashSet::with_capacity(10000),
    ];

    let mut genes:GeneIds = GeneIds::new(9); // split them into 9 bp kmers

    let mut expr_file = parse_fastx_file(&opts.expression).expect("valid path/file");

    if  opts.specie.eq("human") {
        // get all the human sample IDs into this.
        genes.add( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", "Sample1".to_string() );
        genes.add( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", "Sample2".to_string() );
        genes.add( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", "Sample3".to_string() );
        genes.add( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", "Sample4".to_string() );
        genes.add( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", "Sample5".to_string() );
        genes.add( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", "Sample6".to_string() );
        genes.add( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", "Sample7".to_string() );
        genes.add( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", "Sample8".to_string() );
        genes.add( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", "Sample9".to_string() );
        genes.add( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", "Sample10".to_string() );
        genes.add( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", "Sample11".to_string() );
        genes.add( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", "Sample12".to_string() );

    }
    else if opts.specie.eq("mouse") {
        // and the mouse ones
        genes.add( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", "Sample1".to_string() );
        genes.add( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", "Sample2".to_string() );
        genes.add( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", "Sample3".to_string() );
        genes.add( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", "Sample4".to_string() );
        genes.add( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", "Sample5".to_string() );
        genes.add( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", "Sample6".to_string() );
        genes.add( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", "Sample7".to_string() );
        genes.add( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", "Sample8".to_string() );
        genes.add( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", "Sample9".to_string() );
        genes.add( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", "Sample10".to_string() );
        genes.add( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", "Sample11".to_string() );
        genes.add( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", "Sample12".to_string() );

    } else {
        println!("Sorry, but I have no primers for species {}", &opts.specie);
        std::process::exit(1)
    }
    
    while let Some(e_record) = expr_file.next() {
        let seqrec = e_record.expect("invalid record");
        match std::str::from_utf8(seqrec.id()){
            Ok(st) => {
                for id in st.to_string().split("|"){
                    genes.add( &seqrec.seq(), id.to_string() );
                    break;
                }
            },
            Err(err) => eprintln!("The expression entry's id could not be read: {}", err),
        }
    }

    let mut ab_file = parse_fastx_file(&opts.antybody).expect("valid path/file");
    while let Some(ab_record) = ab_file.next() {
        let seqrec = ab_record.expect("invalid record");
        match std::str::from_utf8(seqrec.id()){
            Ok(st) => {
                for id in st.to_string().split("|"){
                    genes.add( &seqrec.seq(), id.to_string() );
                    break;
                }
            },
            Err(err) => eprintln!("The expression entry's id could not be read: {}", err),
        }
    }

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new();

    // that is a class to strore gene expression data.
    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
    // and umi's are u64 again
    let mut gex = CellIds10x::new( 9 );

    let mut unknown = 0;
    let mut no_sample = 0;

    {
        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
        let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

        'main: while let Some(record2) = readefile.next() {
            if let Some(record1) = readereads.next() {
                let seqrec = record2.expect("invalid record");
                let seqrec1 = record1.expect("invalid record");

                // totally unusable sequence
                // as described in the BD rhapsody Bioinformatic Handbook
                // GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf (google should find this)
                if seqrec1.seq().len() < 66 {
                    unknown +=1;
                    continue;
                }
                if seqrec.seq().len() < 64 {
                    unknown +=1;
                    continue;
                }
                //let seq = seqrec.seq().into_owned();
                for nuc in &seqrec1.seq()[52..60] {  
                    if *nuc ==b'N'{
                        unknown +=1;
                        continue 'main;
                    }
                }
                let umi = Kmer::from( &seqrec1.seq()[52..60]).into_u64();

                // first match the cell id - if that does not work the read is unusable
                match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                    Ok(cell_id) => {
                        match genes.get( &seqrec.seq() ){
                            Some(gene_id) =>{
                                match gex.get( cell_id as u64, format!( "Cell{}", cell_id ) ){
                                    Ok(cell_info) => {
                                        cell_info.add( gene_id, umi );
                                    },
                                    Err(err) => panic!("Could not add a gene expression: gene_id = {}, umi = {} and err= {}",gene_id, Kmer::from( &seqrec1.seq()[52..60]).into_u64(),  err),
                                }
                                //println!("Cool I got a gene id: {}", id);
                            },
                            None => {
                                unknown +=1;
                                eprintln!("I could not identify a gene in this read: {:?}", std::str::from_utf8( &seqrec.seq() ) );
                            }
                        };

                    },
                    Err(_err) => {
                        // cell id could not be recovered
                        // println!("Cell ID could not be recovered from {:?}:\n{}", std::str::from_utf8(&seqrec1.seq()), _err);
                        no_sample +=1;
                        continue
                    }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                };
            } else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        println!( "collected sample info:");


        let fp1 = PathBuf::from(opts.reads.clone());
        println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
        let file_path = PathBuf::from(&opts.outpath).join(
            format!("BD_Rhapsody_expression_and_sample_ids.{}.tsv", fp1.file_name().unwrap().to_str().unwrap() )
        );

        match gex.write ( file_path, &genes ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

        println!(     "no sample ID reads: {} reads", no_sample );
        println!(     "genomic reads     : {} reads", unknown );

        println!( "written to {:?}", file_path.into_os_string().into_string() );
    }


    match now.elapsed() {
       Ok(elapsed) => {println!("quantifyRhapsody finished in {} sec", elapsed.as_secs());},
       Err(e) => {println!("Error: {e:?}");}
    }

}

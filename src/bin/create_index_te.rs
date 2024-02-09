use clap::Parser;

use regex::Regex;

use rustody::fast_mapper::FastMapper;
use rustody::gene_family::GeneFamily;
use rustody::gene::Gene;

use needletail::parse_fastx_file;

//use rustody::sampleids::SampleIds;
//use rustody::analysis::

use std::path::PathBuf;
use std::fs;
//use std::path::Path;

//use std::time::SystemTime;

use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use flate2::read::GzDecoder;
use std::collections::HashMap;
//use std::collections::HashSet;

// use std::convert::TryInto;

//use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;
use indicatif::MultiProgress;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use rustody::mapping_info::MappingInfo;
use rustody::ofiles::Ofiles;

/// Create a binary indes for the quantify_rhapsody program.
/// In contrast to the not te version of this tool we here expect a lot of overlap between different genes.
/// So instead of actively removing sequences that mapp to multiple genes we here try to identify the family of the gene.
/// We still only report unique mapper combinations, but we here allow for family or class specific mappers, too. 
/// Make sure to not cat any gzipped files for Rust.
/// At least here we will only read from the first stream in a catted gz file!

#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the gff gene information table
    #[clap(default_value= "testData/TE_quant/chr14_KI270722v1_random_TEtranscripts.gtf.gz",short, long)]
    gtf: String,
    /// the fasta genome data
    #[clap( default_value=  "testData/TE_quant/chr14_KI270722v1_random.fa.gz",short, long)]
    file: String,
    /// the outpath
    #[clap(default_value=  "testData/TEmapperTest",short, long)]
    outpath: String,
    /// the mapping kmer length
    #[clap(default_value_t=32, long)]
    gene_kmers: usize,
    /// create text outfile instead of binary
    #[clap(default_value_t=false, long)]
    text: bool,
    /// the string to check for family level names
    #[clap(default_value="family_id", long)]
    family: String,
    /// the string to check for class level names
    #[clap(default_value="class_id", long)]
    class: String,
    /// the string to check for gene levels names
    #[clap(default_value="gene_id", long)]
    gene: String,
    /// the string to check for transcript levels names
    #[clap(default_value="transcript_id", long)]
    transcript: String,
    /// the maximum area of the transcript tend to be indexed (useful for single cell transcriptomes like 10x or Rhapsody)
    #[clap(default_value_t=200, long)]
    max_length: usize,
    /// the per processor chunk size (default 100_000)
    #[clap(default_value_t=100, long)]
    chunk_size: usize,
    /// processor cores to use (default all)
    #[clap( short)]
    n: Option<usize>,
}


/*
re_family_name =  Regex::new(r#"family_id "([\(\)/\w\d\-\._]*)""#).unwrap();
re_gene_id = Regex::new(r#"gene_id "([\(\)/\w\d\-\._]*)";"#).unwrap();
re_transcript_id = Regex::new(r#"transcript_id "([\(\)/\w\d\-\._]*)";"#).unwrap();
re_class_id = Regex::new(r#"class_id "([\(\)/\w\d\-\._\?]*)";"#).unwrap();

    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
*/


fn process_lines ( lines:&&[String], index: &mut FastMapper ,seq_records: &HashMap< String, Vec<u8>>, max_length:usize, re_class_id: &Regex, re_family_name: &Regex, re_gene_id: &Regex, re_transcript_id: &Regex){

    let mut families = HashMap::<String, GeneFamily>::new();

    for line in lines.iter() {

        let parts: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if parts.len() < 8{
            continue;
        }

        if parts[2] == "exon"{
            // For this example the exons info is all there is in the gtf file
            // so here I need to treat this as a whole gene
            let transcript_id = match re_transcript_id.captures(&parts[8].to_string()) {
                Some(captures) => {
                    let transcript_id = captures.get(1).unwrap().as_str().to_string();
                    transcript_id
                },
                None => {
                    eprintln!("I could not get a transcript id from this line!\n{:?}", &parts[8].to_string() ); //hope I nerver get to see that
                    continue;
                }
            };
            let family_name = match re_family_name.captures(&parts[8].to_string()) {
                Some(captures) => {
                    let transcript_id = captures.get(1).unwrap().as_str().to_string();
                    transcript_id
                },
                None => {
                    eprintln!("I could not get a family_name from this line!\n{:?}", &parts[8].to_string() ); //hope I nerver get to see that
                    continue;
                }
            };
            let gene_name = match re_gene_id.captures(&parts[8].to_string()) {
                Some(captures) => {
                    let transcript_id = captures.get(1).unwrap().as_str().to_string();
                    transcript_id
                },
                None => {
                    eprintln!("I could not get a gene_name from this line!\n{:?}", &parts[8].to_string() ); //hope I nerver get to see that
                    continue;
                }
            };
            let class_name = match re_class_id.captures(&parts[8].to_string()) {
                Some(captures) => {
                    let transcript_id = captures.get(1).unwrap().as_str().to_string();
                    transcript_id
                },
                None => {
                    eprintln!("I could not get a class_name from this line!\n{:?}", &parts[8].to_string() ); //hope I nerver get to see that
                    continue;
                }
            };
            // for this approach we need to use the family model?
            //eprintln!("I will add these gene names: {}, {}, {}, {}", transcript_id.to_string(), gene_name.to_string(), family_name.to_string(), class_name.to_string() );
            let mut gene = Gene::new(
                parts[0].to_string(),
                parts[3].to_string(),
                parts[4].to_string(),
                parts[6].to_string(),
                transcript_id.to_string(), // this will be the id we add to the index
                vec![transcript_id.to_string(), gene_name.to_string(), family_name.to_string(), class_name.to_string()],
            );
            gene.add_exon( parts[3].to_string(),parts[4].to_string());
            
            match families.get_mut( &family_name ){
                Some( family_object ) => {
                    family_object.push( gene );
                },
                None => {
                    let family_object : GeneFamily = 
                    GeneFamily::new( family_name.to_string(), gene );
                    families.insert( family_name.to_string(), family_object);
                }
            }
        }
    }

    for (_, family) in &families {
        // Do something with the gene, e.g. remove it
        family.index( index, max_length , &seq_records );
    }

}

/// One batch consists of several lines from the gtf or gff file.
/// They will be used to select sequences from the seq_records and
/// add the mapping indexes for these sequences to the index.
fn process_batch(
    lines: &Vec<String>, // No longer mutable
    seq_records: &HashMap< String, Vec<u8>>, // Assuming SeqRecords is a custom type    
    index: &mut FastMapper, // Assuming PartialIndex is a custom type
    report: &mut MappingInfo, // Assuming Report is a custom type
    num_threads: usize,
    opts: &Opts, // Assuming Opts is a custom type
    re_class_id: &Regex, // Assuming Regex is a custom type for regular expressions
    re_family_name: &Regex,
    re_gene_id: &Regex,
    re_transcript_id: &Regex,
) {
    report.stop_file_io_time();
    eprintln!("creating mapper");

    let results: Vec<FastMapper> = lines
        .par_chunks(lines.len() / num_threads + 1)
        .map(|data_split| {
            let mut idx = FastMapper::new(32, 10, 0 );
            let _res = process_lines(
                &data_split,
                &mut idx,
                &seq_records,
                opts.max_length,
                &re_class_id,
                &re_family_name,
                &re_gene_id,
                &re_transcript_id,
            );
            idx.make_index_te_ready_single();
            idx
        })
        .collect();

    report.stop_multi_processor_time();

    for idx in results {
        index.merge(idx);
    }

    report.stop_single_processor_time();
}


// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    //let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    if fs::metadata(&opts.outpath).is_ok() {
        if let Err(err) = fs::remove_dir_all(&opts.outpath) {
            eprintln!("Error old index directory: {}", err);
        } else {
            println!("Old index directory removed successfully!");
        }

        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory: {}", err);
        } else {
            println!("New index directory created successfully!");
        }
    }else {
        if let Err(err) = fs::create_dir_all(&opts.outpath) {
            eprintln!("Error creating directory: {}", err);
        } else {
            println!("New index directory created successfully!");
        }
    }

    let mut kmer_size = opts.gene_kmers;
    if opts.gene_kmers > 32{
        eprintln!("Sorry the max size of the kmers is 32 bp");
        kmer_size = 32;
    }

    //const COVERED_AREA:usize = 400; // cover 600 bp of the transcript
    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");

    let log_file_str = PathBuf::from(&opts.outpath).join(
        "index_log.txt"
    );

    println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );

    let log_file = match File::create( log_file_str ){
        Ok(file) => file,
        Err(err) => {
            panic!("Error: {err:#?}" );
        }
    };
    let ofile = Ofiles::new( 1, "NOT_USED", "A.gz", "B.gz",  opts.outpath.as_str() );

    let mut report = MappingInfo::new( Some(log_file), 32.0 , 0, Some(ofile) );

    report.start_counter();
    // read the fasta data in!
    let mut seq_records = HashMap::< String, Vec<u8>>::new();

    let mut expr_file = parse_fastx_file(&opts.file).expect("valid path/file");
    let mut id: &[u8];
    let delimiter: &[u8] = b"  ";  // The byte sequence you want to split by


    while let Some(e_record) = expr_file.next() {
        let seqrec = e_record.expect("invalid record");
        id = seqrec.id().split(|&x| x == delimiter[0]).collect::<Vec<&[u8]>>()[0];
        //id = seqrec.id().split(|&x| x == delimiter[0]).collect()[0];
        //eprintln!("'{}'", std::str::from_utf8(id).unwrap() );
        seq_records.insert( std::str::from_utf8( id ).unwrap().to_string() , seqrec.seq().to_vec());
    }

    /*
    The TE transcrip gtf I have looks like that:
    chr14_KI270722v1_random hg38_rmsk       exon    100070  100378  2251    +       .       gene_id "AluSp"; transcript_id "AluSp_dup53774"; family_id "Alu"; class_id "SINE";
    chr14_KI270722v1_random hg38_rmsk       exon    100379  101084  12116   -       .       gene_id "Tigger1"; transcript_id "Tigger1_dup12525"; family_id "TcMar-Tigger"; class_id "DNA";
    chr14_KI270722v1_random hg38_rmsk       exon    101085  101457  1599    -       .       gene_id "L1MCa"; transcript_id "L1MCa_dup7723"; family_id "L1"; class_id "LINE";
    */

    // and init the Index:

    let mut index = FastMapper::new( kmer_size, 900_000, 0 );

    // in short I need to get an internal model of a gene to work.
    // I want to know where the gene starts ans ends (likely transcripts)
    // I want to get intron/exon bounds and then add both a spliced version of the transcript 
    // as well as an unspliced version

    /*
    seqname   0
    source    1    
    Feature   2
    start     3
    end       4
    score     5
    strand    6
    frame     7
    attribute 8 
    */
    
    let gtf = Regex::new(r".*gtf.?g?z?$").unwrap();

    let re_family_name: Regex;
    let re_gene_id: Regex;
    let re_transcript_id : Regex;
    let re_class_id : Regex;

    println!("Starting to collect the data {}", report.now() );


    // gene_id "AluSp"; transcript_id "AluSp_dup53774"; family_id "Alu"; class_id "SINE";
    match gtf.is_match( &opts.gtf ){
        true => {
            eprintln!("gtf mode");
            let family_regex = format!(r#"{} "([\(\)/\w\d\-\._]*)""#, opts.family );
            re_family_name =  Regex::new(&family_regex).unwrap();

            let gene_regex = format!(r#"{} "([\(\)/\w\d\-\._]*)""#, opts.gene );
            re_gene_id = Regex::new(&gene_regex).unwrap();

            let transcript_regex = format!(r#"{} "([\(\)/\w\d\-\._]*)""#, opts.transcript );
            //panic!("{:?}", transcript_regex);
            re_transcript_id = Regex::new(&transcript_regex).unwrap();

            let class_regex = format!(r#"{} "([\(\)/\w\d\-\._\?]*)""#, opts.class );
            re_class_id = Regex::new(&class_regex).unwrap();
        },
        false => {
            eprintln!("gff mode.");
            let family_regex = format!(r#"{}=([\(\)/\w\d\-\._]*);"#, opts.family );
            re_family_name =  Regex::new(&family_regex).unwrap();

            let gene_regex = format!(r#"{}=([\(\)/\w\d\-\._]*);"#, opts.gene );
            re_gene_id = Regex::new(&gene_regex).unwrap();

            let transcript_regex = format!(r#"{}=([\(\)/\w\d\-\._]*);"#, opts.transcript );
            //panic!("{:?}", transcript_regex);
            re_transcript_id = Regex::new(&transcript_regex).unwrap();

            let class_regex = format!(r#"{}=([\(\)/\w\d\-\._\?]*);"#, opts.class );
            re_class_id = Regex::new(&class_regex).unwrap();
        },
    };
    //let mut gene_id:String;
    //let mut gene_name:String;
    //let mut transcript_id:String;

    // now we need to read the gtf info and get all families out of that
    // we specificly need the last 64bp of the end of the gene as that is what we are going to map
    // if there is a splice entry in that area we need to create a spliced and an unspliced entry

    if ! opts.gtf.ends_with(".gz") {
        panic!("Please gzip your gtf file - thank you! {}", opts.gtf.to_string());
    }

    let f1 = match File::open( &opts.gtf ){
        Ok(file) => file,
        Err(err) => panic!("The file {} does not exists: {err}", &opts.gtf ),
    };
    let file1 = GzDecoder::new(f1);
    let reader = BufReader::new( file1 );
    //let mut missing_chr:HashSet<String>  = HashSet::new();

    let num_threads = opts.n.unwrap_or(num_cpus::get());

    println!("I am using {num_threads} cores");

    let m = MultiProgress::new();
    let pb = m.add(ProgressBar::new(5000));
    let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");
    pb.set_style(spinner_style);

    let reads_per_chunk = opts.chunk_size;
    let mut good_read_count = 0;
    let max_dim = reads_per_chunk * num_threads;
    let mut lines = Vec::<String>::with_capacity( max_dim );
    eprintln!("In one batch I will analyze {} elemets {} cores x {} elements per batch", max_dim, num_threads, reads_per_chunk);
    let mut batch = 0;

    let mut last_chr = "".to_string();

    let mut partial_index = FastMapper::new( kmer_size, 900_000,0 );
    report.start_ticker();
    for line in reader.lines() {
        let rec = line.ok().expect("Error reading record.");
        let parts: Vec<&str> = rec.split('\t').collect();

        if last_chr != parts[0].to_string() {
            if good_read_count > 1 {
                
                process_batch(
                    &lines, &seq_records, &mut partial_index, &mut report, num_threads, &opts, 
                    &re_class_id, &re_family_name, &re_gene_id, &re_transcript_id,
                );
                
                batch +=1;
                let (_h, m, s, ms) = report.stop_ticker();
                println!(
                    "batch {}: For {} sequence regions we needed {m}min {s}sec {ms}msec",
                    batch, lines.len()
                );
                report.start_ticker();
                lines.clear();
                good_read_count = 0;
            }
            println!("Chromosome {last_chr} is finished - integrating it");

            partial_index.make_index_te_ready();
            partial_index.write_index( format!("{}/{}/",opts.outpath, last_chr )).unwrap();

            index.merge( partial_index );
            partial_index = FastMapper::new( kmer_size, 900_000, 0 );

            let (_h, m, s, ms) = report.stop_ticker();
            println!(
                "batch {}: For {} sequence regions we needed {m}min {s}sec {ms}msec",
                batch, lines.len()
            );
            report.start_ticker();
            last_chr = parts[0].to_string();
        }
        if good_read_count < max_dim{
            lines.push(rec);
            good_read_count+=1;
        }
        else {

            process_batch(
                &lines, &seq_records, &mut partial_index, &mut report, num_threads, &opts, 
                &re_class_id, &re_family_name, &re_gene_id, &re_transcript_id,
            );
            
            batch +=1;
            let (_h, m, s, ms) = report.stop_ticker();
            println!(
                "batch {}: For {} sequence regions we needed {m}min {s}sec {ms}msec",
                batch, lines.len()
            );
            report.start_ticker();
            lines.clear();
            good_read_count = 0;
        }
    }

    if good_read_count > 0{
        //good_read_count = 0;
        println!("Rounding up the last {good_read_count} values");

        process_batch(
            &lines, &seq_records, &mut partial_index, &mut report, num_threads, &opts, 
            &re_class_id, &re_family_name, &re_gene_id, &re_transcript_id,
        );
        
        batch +=1;
        let (_h,m, s, ms) = report.stop_ticker();
        println!(
            "batch {}: For {} sequence regions we needed {m}min {s}sec {ms}msec",
            batch, lines.len()
        );
        report.start_ticker();
        lines.clear();
        //good_read_count = 0;
    }

    partial_index.make_index_te_ready();
    index.merge( partial_index );
    let (h,m,s,_ms) = MappingInfo::split_duration( report.absolute_start.elapsed().unwrap() );

    println!("Collapsing the index down to a usable size? at {h}h {m}min {s}sec.");

    index.make_index_te_ready();


    eprintln!(" total first keys {}\n total second keys {}\n total single gene per second key {}\n total multimapper per second key {}", index.info()[0], index.info()[1], index.info()[2], index.info()[3] );

    index.write_index( opts.outpath.to_string() ).unwrap();

    if opts.text{
        index.write_index_txt( opts.outpath.to_string() ).unwrap();
    }

    //index.write_index_txt( opts.outpath.to_string() ).unwrap();
    //eprintln!("THIS IS STILL IN TEST MODE => TEXT INDEX WRITTEN!!! {}",opts.outpath.to_string() );
    
    eprintln!("{}", report.program_states_string() );

    eprintln!("Index stored in folder {}", opts.outpath );

    println!("finished at {}\n", report.now() );
    // match now.elapsed() {
    //     Ok(elapsed) => {
    //         let mut milli = elapsed.as_millis();

    //         let mil = milli % 1000;
    //         milli= (milli - mil) /1000;

    //         let sec = milli % 60;
    //         milli= (milli -sec) /60;

    //         let min = milli % 60;
    //         milli= (milli -min) /60;

    //         eprintln!("finished in {milli} h {min} min {sec} sec {mil} milli sec");
    //     },
    //     Err(e) => {println!("Error: {e:?}");}
    // }

}

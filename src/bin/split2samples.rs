use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
// use flate2::read::GzDecoder;

use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::io::BufWriter;
// use std::io::BufReader;
// use std::str;

use this::cellids::CellIds;
use std::collections::BTreeMap;

use this::sampleids::SampleIds;

use std::path::PathBuf;
use std::fs::File;
use std::fs;

use std::io::prelude::*;
use std::time::SystemTime;

//use crate::glob;
use glob::glob; // to search for files
use regex::Regex;

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
    /// the mode of the program [fastqSplit, cellIdent, sampleSplit]
    #[clap(short, long)]
    mode: String,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
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
        let file1_path = PathBuf::from(&opts.outpath).join(format!("{}.{}.{}", opts.mode, id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(&opts.outpath).join(format!("{}.{}.{}", opts.mode, id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::create(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {}.{id}.{} cound not be created: {err}", opts.mode,  fp1.file_name().unwrap().to_str().unwrap() )
        };
        let f2 = match File::create(file2_path){
            Ok(file) => file,
            Err(err) =>panic!("The file {}.{id}.{} cound not be created: {err}", opts.mode,  fp2.file_name().unwrap().to_str().unwrap() )
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
            Err(e) => eprintln!("Could not flush R1: {e}"),
        };
        match self.buff2.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R2: {e}"),
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
        Err(e) => panic!("I could not create the outpath: {e}")
    };

    match opts.mode.as_str() {
        "fastqSplit" => fastq_split( &opts), 
        "cellIdent" => cell_ident( &opts) , 
        "sampleSplit" => sample_split( &opts ),
        "analysis" => {
            cell_ident( &opts );

            match now.elapsed() {
                Ok(elapsed) => {println!("split2samples cellIdent finished in {} sec", elapsed.as_secs());},
                Err(e) => {println!("Error: {e:?}");}
            }
            let loc_now = SystemTime::now();
            sample_split( &opts );
            match loc_now.elapsed() {
                Ok(elapsed) => {println!("split2samples sampleSplit finished in {} sec", elapsed.as_secs());},
                Err(e) => {println!("Error: {e:?}");}
            }
        }
        &_ => eprintln!("mode {} is not supported:\nfastqSplit splits a fastq file into 10 subsets\ncellIdent captures the cells<->sample info\nsampleSplit splits the fastq file based on the collective results from the cellIdent runs", &opts.mode) 
    };

    match now.elapsed() {
       Ok(elapsed) => {println!("split2samples {} finished in {} sec", &opts.mode, elapsed.as_secs());},
       Err(e) => {println!("Error: {e:?}");}
    }

}


/// Split both fastq files into 10 smaller fastq files to process them in paralele.
/// This function is less perfcormant than the fastqsplit python package.
/// It is not developed further
fn fastq_split( opts:&Opts){
    // this will at the moment just split it into 10 files.

    println!("This is inefficient - please use another tool! >4 times slower than fastqsplit");

    let split:usize = 400;
    let mut id:usize = 0;
    let mut fileid:usize = 0;

    let mut ofiles: Vec<Ofiles>;
    ofiles = vec![
        Ofiles::new(1, opts),
        Ofiles::new(2, opts),
        Ofiles::new(3, opts),
        Ofiles::new(4, opts),
        Ofiles::new(5, opts),
        Ofiles::new(6, opts),
        Ofiles::new(7, opts),
        Ofiles::new(8, opts),
        Ofiles::new(9, opts),
        Ofiles::new(10, opts)
    ];

    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    while let Some(record2) = readefile.next() {
        if let Some(record1) = readereads.next() {
            let seqrec = record2.expect("invalid record");
            let seqrec1 = record1.expect("invalid record");
            match seqrec1.write(&mut ofiles[ fileid ].buff1, None){
                Ok(_) => (),
                Err(err) => println!("{err}")
            };
            match seqrec.write( &mut ofiles[ fileid ].buff2, None){
                Ok(_) => (),
                Err(err) => println!("{err}")
            };
            id += 1;
            if id > split{
                id = 0;
                fileid += 1;
                if fileid == ofiles.len(){
                    fileid = 0;
                }
            }
        }
        else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
    }
}


/// cell_ident idetifes the cells in a sample. It therfore checks the R2 reads for the sample strings and 
/// identifies the cellid for these reads only. The cell ids are written to simple human readable text files with one cellid per line.
fn cell_ident( opts:&Opts ){
    // for now just write some tests to the stdout!
    // thanks to Rob I can now skip the stdout part!
    //let stdout = std::io::stdout();
    //let lock = stdout.lock();
    //let buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    println!("starting to collect the cell ids per sample - can take a LONG time");

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

    // let mut file1_path = PathBuf::from(&opts.outpath).join("ambig.R1.fq.gz");
    // let mut file2_path = PathBuf::from(&opts.outpath).join("ambig.R2.fq.gz");

    //  now we need to get a CellIDs object, too
    let cells = CellIds::new( &opts.version, 9);

    let mut unknown = 0;
    let mut no_sample = 0;

    {
        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
        let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

        while let Some(record2) = readefile.next() {
            if let Some(record1) = readereads.next() {
                let seqrec = record2.expect("invalid record");
                let seqrec1 = record1.expect("invalid record");

                // totally unusable sequence
                if seqrec1.seq().len() < 53 {
                    unknown +=1;
                    continue;
                }
                //let seq = seqrec.seq().into_owned();

                // first match the cell id - if that does not work the read is unusable
                match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                    Ok(val) => {
                        // OK the read is usable - check if it is a sample
                        match samples.get( &seqrec.seq(), 9, 10 ){
                            Ok(id) => {
                                cell2sample[id as usize].insert( val ); // will never insert one element twice. Great!
                            },
                            Err(_err) => {
                                // so here we add the new gene detection!

                                //println!("{}",err);
                                //seqrec1.write(&mut file1_ambig_out, None)?;
                                //seqrec.write(&mut file2_ambig_out, None)?;
                                //seqrec1.write(&mut outBuff1, None)?;
                                //seqrec.write(&mut outBuff2, None)?;
                                unknown +=1;
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
        for (i, entry ) in cell2sample.iter().enumerate(){
            if ! entry.is_empty() {

                match samples.read_get( i as u32 )  {
                    Some(val) => println!( "    sample {}: {} reads and {} cells", i+1, val.total, cell2sample[i].len() ),
                    None => println!( "    sample {}: na reads and {} cells", i+1,  cell2sample[i].len() ),
                };

                let fp1 = PathBuf::from(opts.reads.clone());
                //println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
                let file_path = PathBuf::from(&opts.outpath).join(
                    format!("{}.sample{}.ints.txt", fp1.file_name().unwrap().to_str().unwrap(), i+1)
                );
                let file = match File::create( file_path ){
                    Ok(file) => file,
                    Err(err) => {
                        eprintln!("Error: {err:#?}" );
                        std::process::exit(1)
                    }
                };
                let mut writer = BufWriter::new(&file);
                for int in &cell2sample[i]{
                    //file.write( )
                    //BigEndian::write_u32(file, int).unwrap();
                    writeln!(writer, "{int}").unwrap();
                    //println!( "        with sample {}",int );
                }
            }
        }
        println!(     "no sample ID reads: {no_sample} reads" );
        println!(     "genomic reads     : {unknown} reads" );
    }

}

/// This function does the main work. It takes about three times the time to finish as the cell_ident function.
/// All in all the combination of the two functions needed about 4 to 5 hours to finish two fastq.gz files of 
fn sample_split( opts: &Opts){

    println!("Splitting the real data");

    let mut cells2sample =  BTreeMap::<u32, u32>::new(); //: BTreeMap<u32, u32>;

    //  now we need to get a CellIDs object, too
    let cells = CellIds::new(&opts.version, 9);

    // and we need to collect all cell2sample data
    let re = Regex::new(r"sample(\d*).int").unwrap();
    for entry in glob(&(opts.outpath.clone()+"/*.ints.txt")).expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => {
                // get the id of the sample we read in at the moment
                // let re = Regex::new(r"(\d{4})-(\d{2})-(\d{2})").unwrap();
                // let text = "2012-03-14, 2013-01-01 and 2014-07-05";
                // for cap in re.captures_iter(text) {
                //     println!("Month: {} Day: {} Year: {}", &cap[2], &cap[3], &cap[1]);
                // }
                let mut id_loc:u32 = 0;
                let file_name = match path.to_str(){
                    Some(file_name) => file_name,
                    None => {
                        continue
                    },
                };
                for cap in re.captures_iter( file_name ) {
                    match cap[1].to_string().parse(){
                        Ok(id) => {
                            id_loc = id;
                            let data = fs::read_to_string(path.clone());
                            for st in data.expect("std::io::Error").lines().map(|line| line.trim()){
                                cells2sample.insert( st.parse::<u32>().unwrap(), id );
                            }
                        },
                        Err(err) => {
                            eprintln!("file {file_name} does not contain the sample id! patternmatch threw this error {err}" );
                            std::process::exit(1)
                        }
                    }
                }
                                
                if id_loc == 0{
                    eprintln!("file {file_name} does not contain the sample id! patternmatch did not match" );
                    std::process::exit(1)
                }
            },
            Err(e) => println!("{e:?}"),
        }
    }

    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    let mut unknown = 0;

    let mut ofiles: Vec<Ofiles>;
    ofiles = vec![
        Ofiles::new(1, opts),
        Ofiles::new(2, opts),
        Ofiles::new(3, opts),
        Ofiles::new(4, opts),
        Ofiles::new(5, opts),
        Ofiles::new(6, opts),
        Ofiles::new(7, opts),
        Ofiles::new(8, opts),
        Ofiles::new(9, opts),
        Ofiles::new(10, opts),
        Ofiles::new(11, opts),
        Ofiles::new(12, opts)
    ];

    let mut missing:u32 = 0;
    while let Some(record2) = readefile.next() {
        if let Some(record1) = readereads.next() {
            let seqrec = record2.expect("invalid record");
            let seqrec1 = record1.expect("invalid record");
            //let seq = seqrec.seq().into_owned();
            if seqrec1.seq().len() < 53 {
                unknown +=1;
                continue;
            }

            match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                Ok(id) => {
                    // check if the cell has been detected before or print error
                    match cells2sample.get( &id ){
                        Some(cell_id) => {
                            // print record to sample specific out file
                            //println!("Cell {} is written to file", id);
                            let loc_id:usize = *cell_id as usize -1;
                            ofiles[loc_id].count += 1;
                            // here we could possibly 'fix' the sequence. Look into the splitp source on how to do that!
                            ofiles[loc_id].buff1.write_all(b"@").unwrap();
                            ofiles[loc_id].buff1.write_all(seqrec1.id()).unwrap();
                            ofiles[loc_id].buff1.write_all(b"\n").unwrap();
                            // here comes the seq
                            let mut seq = seqrec1.seq().into_owned();
                            let expected = cells.to_sequence( id );
                            seq[..9].copy_from_slice(&expected[0][..9]);
                            seq[21..30].copy_from_slice(&expected[1][..9]);
                            seq[43..52].copy_from_slice(&expected[2][..9]);

                            ofiles[loc_id].buff1.write_all( &seq ).unwrap();
                            ofiles[loc_id].buff1.write_all(b"\n+\n").unwrap();
                            ofiles[loc_id].buff1.write_all(seqrec1.qual().unwrap()).unwrap();


                            //seqrec1.write(&mut ofiles[*cell_id as usize].file1, None)?;
                            match seqrec.write( &mut ofiles[loc_id].buff2, None){
                                Ok(_) => (),
                                Err(e) => eprintln!("sequence could not be written {e}")
                            };
                        },
                        None => {
                            missing += 1;
                            //let tmp = seqrec1.seq();
                            //let seq = str::from_utf8( &tmp );
                            //println!( "cell with id {} has not been detected before:\n{:?}", id, seq )
                        }
                    }
                }
                Err(_err) => {
                    unknown += 1;
                    //println!("trying to print the data, but: {}",err);
                    continue
                }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
            };
        } else {
            println!("file 2 had reads remaining, but file 1 ran out of reads!");
        }
    }

    let mut i = 1;
    for mut f in ofiles{
        println!( "{} reads in sample {i}", f.count );
        f.close();
        i += 1;
    }
    println!("{missing} reads did not match a sample");
    println!("{unknown} reads had no cell id" );

    //Ok(())
}




use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;

use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::io::BufWriter;
use std::str;

mod cellids;
use crate::cellids::CellIds;
use std::collections::BTreeMap;

mod sampleids;
use crate::sampleids::SampleIds;

use std::path::PathBuf;
use std::fs::File;
use std::fs;

use std::io::prelude::*;
use std::time::SystemTime;

crate glob;
use glob::glob; // to search for files

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
}

struct Ofiles {
    pub id: usize,
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufWriter<GzEncoder<File>>,
    pub buff2: BufWriter<GzEncoder<File>>
}

impl Ofiles{
    pub fn new(id: usize, opts: &Opts )->Self {
        let file1_path = PathBuf::from(opts.outpath).join(format!("{}.{}.{}", opts.mode, id, opts.reads ));
        let file2_path = PathBuf::from(opts.outpath).join(format!("{}.{}.{}", opts.mode, id, opts.file ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let file1 = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
        let file2 = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());

        let mut buff1 = BufWriter::new( file1 );
        let mut buff2 = BufWriter::new( file2 );

        let count:u32 = 0;
        Self{
            id,
            count,
            //file1,
            //file2,
            buff1,
            buff2
        }
    }
    pub fn close( &mut self ){

        self.buff1.flush();
        self.buff2.flush();

        // let mut file1 = self.buff1.into_inner().unwrap();
        // let mut file2 = self.buff2.into_inner().unwrap();
         
        // match file1.try_finish(){
        //     Ok(_val) => (), //println!("closed sample {} R1 file with {} reads",self.id, self.count),
        //     Err(err) => println!("I could not close the sample {} R1 file with {} reads\nErr {}", self.id, self.count, err)
        // };
        // match file2.try_finish(){
        //     Ok(_val) => (),//println!("closed sample {} R2 file with {} reads",self.id, self.count),
        //     Err(err) => println!("I could not close the sample {} R2 file with {} reads\nErr {}", self.id, self.count, err)
        // };
    }
}


fn main() -> anyhow::Result<()> {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    match mode.as_str() {
        "fastqSplit" => fastqSplit( &opts), 
        "cellIdent" => cellIdent( &opts) , 
        "sampleSplit" => sampleSplit ( &opts )
    };

    match now.elapsed() {
       Ok(elapsed) => {
           // it prints '2'
           println!("split2samples {} finished in {} sec", &opts.mode, elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {e:?}");
       }
   }

}


fn fastqSplit( opts:Opts){
    // this will at the moment just split it into 10 files.
    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");
    let split:usize = 100;
    let id:usize = 0;
    let fileid:usize = 0;

    let mut ofiles: Vec<Ofiles>;
    ofiles = vec![
        Ofiles::new(1, &opts),
        Ofiles::new(2, &opts),
        Ofiles::new(3, &opts),
        Ofiles::new(4, &opts),
        Ofiles::new(5, &opts),
        Ofiles::new(6, &opts),
        Ofiles::new(7, &opts),
        Ofiles::new(8, &opts),
        Ofiles::new(9, &opts),
        Ofiles::new(10, &opts)
    ];
    while let Some(record2) = readefile.next() {
        if let Some(record1) = readereads.next() {
            let seqrec = record2.expect("invalid record");
            let seqrec1 = record1.expect("invalid record");
            seqrec1.write(&mut ofiles[ fileid ].buff1, None)?;
            seqrec.write( &mut ofiles[ fileid ].buff2, None)?;
            id += 1;
            if id > split{
                id = 0;
                fileid += 1;
                if file_id > ofiles.len(){
                    fileid = 0;
                }
            }
        }
        else {
                anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
    }
}

fn cellIdent( opts:Opts ){
    // for now just write some tests to the stdout!
    // thanks to Rob I can now skip the stdout part!
    //let stdout = std::io::stdout();
    //let lock = stdout.lock();
    //let buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    println!("starting to collect the cell ids per sample - can take a LONG time");

    let sub_len = 9;
    let mut samples = SampleIds::new( sub_len );// = Vec::with_capacity(12);
 
    fs::create_dir_all(&opts.outpath)?;

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

    if  opts.specie.eq("human") {
        // get all the human sample IDs into this.
        samples.add( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG" );
        samples.add( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG" );
        samples.add( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT" );
        samples.add( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT" );
        samples.add( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG" );
        samples.add( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC" );
        samples.add( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT" );
        samples.add( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT" );
        samples.add( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG" );
        samples.add( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC" );
        samples.add( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC" );
        samples.add( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG" );

    }
    else if opts.specie.eq("mouse") {
        // and the mouse ones
        samples.add( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG" );
        samples.add( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC" );
        samples.add( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC" );
        samples.add( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT" );
        samples.add( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT" );
        samples.add( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC" );
        samples.add( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC" );
        samples.add( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG" );
        samples.add( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC" );
        samples.add( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC" );
        samples.add( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT" );
        samples.add( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC" );

    } else {
        println!("Sorry, but I have no primers for species {}", &opts.specie);
        std::process::exit(1)
    }

    // let mut file1_path = PathBuf::from(&opts.outpath).join("ambig.R1.fq.gz");
    // let mut file2_path = PathBuf::from(&opts.outpath).join("ambig.R2.fq.gz");

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new();

    let mut unknown = 0;

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

                match samples.get( &seqrec.seq(), 9, 10 ){
                    Ok(id) => {
                        match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                            Ok(val) => {
                                cell2sample[id as usize].insert( val ); // will never insert one element twice. Great!
                            },
                            Err(_err) => {
                                //println!("{}",err);
                                continue
                            }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                        };
                    },
                    Err(_err) => {
                        //println!("{}",err);
                        //seqrec1.write(&mut file1_ambig_out, None)?;
                        //seqrec.write(&mut file2_ambig_out, None)?;
                        //seqrec1.write(&mut outBuff1, None)?;
                        //seqrec.write(&mut outBuff2, None)?;
                        unknown +=1;
                    }
                }
            } else {
                anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        println!( "collected sample info:");
        for i in 0..cell2sample.len(){
            if cell2sample[i].len() > 0 {

                match samples.read_get( i as u32 )  {
                    Some(val) => println!( "    sample {}: {} reads and {} cells", i+1, val.total, cell2sample[i].len() ),
                    None => println!( "    sample {}: {} reads and {} cells", i+1, "na", cell2sample[i].len() ),
                };

                let file_path = PathBuf::from(&opts.outpath).join(format!("{}sample{}.ints.txt", opts.reads, i+1));
                let file = File::create( file_path )?;
                let mut writer = BufWriter::new(&file);
                for int in &cell2sample[i]{
                    //file.write( )
                    //BigEndian::write_u32(file, int).unwrap();
                    writeln!(writer, "{}", int)?;
                    //println!( "        with sample {}",int );
                }
            }
        }
        println!(     "genomic reads: {} reads", unknown );
    }

}


fn sampleSplit( opts: Opts){

    println!("Splitting the real data");

    let now = SystemTime::now();

    let sub_len = 9;
    let mut samples = SampleIds::new( sub_len );// = Vec::with_capacity(12);
 
    fs::create_dir_all(&opts.outpath)?;

    let mut cells2sample =  BTreeMap::<u32, u32>::new(); //: BTreeMap<u32, u32>;

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new();

    // and we need to collect all cell2sample data
    let re = Regex::new(r"sample(\d*).int").unwrap();
    for entry in glob(opts.outpath+"/*.ints.txt").expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => {
                // get the id of the sample we read in at the moment
                // let re = Regex::new(r"(\d{4})-(\d{2})-(\d{2})").unwrap();
                // let text = "2012-03-14, 2013-01-01 and 2014-07-05";
                // for cap in re.captures_iter(text) {
                //     println!("Month: {} Day: {} Year: {}", &cap[2], &cap[3], &cap[1]);
                // }
                let id:u32;
                for cap in re.captures_iter(path) {
                    id = id as u32;
                }

                let data = fs::read_to_string(path);
                for st in data.lines().map(|line| line.trim()){
                    cells2sample.insert( st as u32, id );
                }
            },
            Err(e) => println!("{:?}", e),
        }
    }

    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    unknown = 0;

    let mut ofiles: Vec<Ofiles>;
    ofiles = vec![
        Ofiles::new(1, &opts),
        Ofiles::new(2, &opts),
        Ofiles::new(3, &opts),
        Ofiles::new(4, &opts),
        Ofiles::new(5, &opts),
        Ofiles::new(6, &opts),
        Ofiles::new(7, &opts),
        Ofiles::new(8, &opts),
        Ofiles::new(9, &opts),
        Ofiles::new(10, &opts),
        Ofiles::new(11, &opts),
        Ofiles::new(12, &opts)
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
                            ofiles[*cell_id as usize].count += 1;
                            // here we could possibly 'fix' the sequence. Look into the splitp source on how to do that!
                            ofiles[*cell_id as usize].buff1.write_all(b"@").unwrap();
                            ofiles[*cell_id as usize].buff1.write_all(seqrec1.id()).unwrap();
                            ofiles[*cell_id as usize].buff1.write_all(b"\n").unwrap();
                            // here comes the seq
                            let mut seq = seqrec1.seq().into_owned();
                            let expected = cells.to_sequence( id );
                            for i in 0..9{
                                seq[i] = expected[0][i];
                            }
                            for i in 21..30{
                                seq[i] = expected[0][i-21];
                            }
                            for i in 43..52{
                                seq[i] = expected[0][i-43];
                            }
                            ofiles[*cell_id as usize].buff1.write_all( &seq ).unwrap();
                            ofiles[*cell_id as usize].buff1.write_all(b"\n+\n").unwrap();
                            ofiles[*cell_id as usize].buff1.write_all(seqrec1.qual().unwrap()).unwrap();


                            //seqrec1.write(&mut ofiles[*cell_id as usize].file1, None)?;
                            seqrec.write( &mut ofiles[*cell_id as usize].buff2, None)?;
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
            anyhow::bail!("file 2 had reads remaining, but file 1 ran out of reads!");
        }
    }

    let mut i = 1;
    for mut f in ofiles{
        println!( "{} reads in sample {}", f.count, i);
        f.close();
        i += 1;
    }
    println!("{} reads did not match a sample",missing);
    println!("{} reads had no cell id", unknown);

    Ok(())
}




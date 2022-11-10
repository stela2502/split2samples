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
    /// the mode of the program [fastqSplit, cellIdent, toSamples]
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
    pub fn new(id: usize, outpath: &str )->Self {
        let file1_path = PathBuf::from(outpath).join(format!("{}.R1.fq.gz", id));
        let file2_path = PathBuf::from(outpath).join(format!("{}.R2.fq.gz", id));
        
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

    println!("starting to collect the cell ids per sample - can take a LONG time");
    let now = SystemTime::now();
    println!("Start time {:#?}", now );

    let opts: Opts = Opts::parse();

    match mode.as_str() {
        "fastqSplit" => (), 
        "cellIdent" => (), 
        "toSamples" => ()
    };


    // for now just write some tests to the stdout!
    // thanks to Rob I can now skip the stdout part!
    //let stdout = std::io::stdout();
    //let lock = stdout.lock();
    //let buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    let sub_len = 9;
    let mut samples = SampleIds::new( sub_len );// = Vec::with_capacity(12);
    // //let File1 = Path::new(outpath).join( Path::new(reads).file_name());
    // let mut File1 = Path::new(&opts.outpath);

    // let mut File2 = File1.join( Path::new(&opts.reads).file_name().unwrap());

    // let file1 = get_writer( &File2 );
    // File2 = File1.join( Path::new(&opts.file).file_name().unwrap() );
    // let file2 = get_writer( &File2);
 
    fs::create_dir_all(&opts.outpath)?;

    let mut cells2sample =  BTreeMap::<u32, u32>::new(); //: BTreeMap<u32, u32>;
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

    let mut file1_path = PathBuf::from(&opts.outpath).join("ambig.R1.fq.gz");
    let mut file2_path = PathBuf::from(&opts.outpath).join("ambig.R2.fq.gz");

    match now.elapsed() {
       Ok(elapsed) => {
           // it prints '2'
           println!("time to prepare the search modules {} sec", elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {e:?}");
       }
   }
    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new();

    let mut unknown = 0;

    {
        // need better error handling here too
        //println!( "why does this file break? {}", file1_path.display() );
        let mut file1_ambig_out = GzEncoder::new(File::create(file1_path).unwrap(), Compression::default());
        let mut file2_ambig_out = GzEncoder::new(File::create(file2_path).unwrap(), Compression::default());


        let mut outBuff1 = BufWriter::new( file1_ambig_out );
        let mut outBuff2 = BufWriter::new( file2_ambig_out );
        
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
                                cells2sample.insert( val, id);
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

                let file_path = PathBuf::from(&opts.outpath).join(format!("sample{}.ints.txt",i+1));
                let file = File::create( file_path )?;
                let mut writer = BufWriter::new(&file);
                for int in &cell2sample[i]{
                    //file.write( )
                    //BigEndian::write_u32(file, int).unwrap();
                    writeln!(writer, "{}", int)?;
                    //println!( "        with sample {}",int );
                }
            }
            // this is not working - I need to first make sure the i32 does work for me here....
            //let file_path = PathBuf::from(&opts.outpath).join("sample{i}.ints.dat");
            //let mut file = File::create( file_path )?;
            //for int in samples[i].hash_set.drain(){
            //    file.write( int )?;
            //}
            
            //samples[i].file1.try_finish()?;
            //samples[i].file2.try_finish()?;
        }
        println!(     "genomic reads: {} reads", unknown );
        outBuff1.into_inner().unwrap().try_finish()?;
        outBuff2.into_inner().unwrap().try_finish()?;
    }
    // now lets rewind the fastq files and actually process the 
    match now.elapsed() {
       Ok(elapsed) => {
           // it prints '2'
           println!("time to identify cells -> sample connections {} sec", elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {e:?}");
       }
   }


    println!("Splitting the real data");

    file1_path = PathBuf::from(&opts.outpath).join("ambig.R1.fq.gz");
    file2_path = PathBuf::from(&opts.outpath).join("ambig.R2.fq.gz");

    //let mut readereads = parse_fastx_file(file1_path).expect("valid path/file");
    //let mut readefile = parse_fastx_file(file2_path).expect("valid path/file");

    let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
    let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

    unknown = 0;

    let mut ofiles: Vec<Ofiles>;
    ofiles = vec![
        Ofiles::new(1, &opts.outpath),
        Ofiles::new(2, &opts.outpath),
        Ofiles::new(3, &opts.outpath),
        Ofiles::new(4, &opts.outpath),
        Ofiles::new(5, &opts.outpath),
        Ofiles::new(6, &opts.outpath),
        Ofiles::new(7, &opts.outpath),
        Ofiles::new(8, &opts.outpath),
        Ofiles::new(9, &opts.outpath),
        Ofiles::new(10, &opts.outpath),
        Ofiles::new(11, &opts.outpath),
        Ofiles::new(12, &opts.outpath)
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

    match now.elapsed() {
       Ok(elapsed) => {
           // it prints '2'
           println!("overall run time {} sec", elapsed.as_secs());
       }
       Err(e) => {
           // an error occurred!
           println!("Error: {e:?}");
       }
   }

    Ok(())
}



// macro_rules! to_cellid {
//    ($r1: expr) => {
//       //1-9,22-30,44-52
//       to_cellid($r1, [0,8], [21,29], [43,51]  )
//    };
// }

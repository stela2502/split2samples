use clap::Parser;
use needletail::parse_fastx_file;
use this::cellids::CellIds;
use this::ofiles::Ofiles;

use std::collections::HashSet;

//use this::sampleids::SampleIds;

//use this::last5::Last5;

use std::fs;

use std::time::SystemTime;

use serde::Deserialize;

// use std::collections::HashSet;
// use std::convert::TryInto;

/// subsets a fastq file set to only those reads containing a cell read.

#[derive(Parser)]
#[clap(version = "0.3.6", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
    /// Optional: end the analysis after processing <max_reads> cell fastq entries 
    #[clap(default_value_t=usize::MAX, short, long)]
    max_reads: usize,
    /// Optional: CellIDs you want to get into the outfile 
    #[clap(default_value="no file", short, long)]
    cell_ids: String,
}

#[derive(Debug, Deserialize)]
struct BCMapRecord {
    id: String,
    cellid: String,
}


fn parse_bc_map(bc_map: &str ) -> HashSet<u32> {

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(true)
        .from_path(bc_map)
        .expect("cannot open cellIDs");

    //let sub_len = 9;
    let mut hash = HashSet::new();

    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BCMapRecord = result.expect("could not deserialize barcode map record.");

        // now I need to get the seq into a [&u8] vector....
        // from https://gist.github.com/jimmychu0807/strinf-conversion.rs 
        // let byte1: Vec<u8> = src1.iter().map(|c| *c as u8).collect::<Vec<_>>();
        //println!("what do I have here: {}", record.bc);
        let cellid: u32 =  record.cellid.parse().unwrap();
        println!("{}", record.id );
        hash.insert( cellid );
    }

    hash
}

// the main function nowadays just calls the other data handling functions
fn main() {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();


    // let max_reads = match &opts.max_reads {
    //     Some(val) => val,
    //     None => &usize::MAX,
    // };

    println!("Analysis will stop after having processed {} fastq entries containing a cell info\n", opts.max_reads);


    println!("init models");    

    fs::create_dir_all(&opts.outpath).expect("AlreadyExists");
    
    let pos:Vec<usize>;
    let min_sizes:Vec<usize>;

    if &opts.version == "v1"{
        pos = vec![0,9, 21,30, 43,52, 52,60 ];
        min_sizes = vec![ 66, 60 ];
    }
    else {
        pos = vec![0,9, 12,22, 26,35 , 36,42 ];
        min_sizes = vec![ 51, 51 ];
    }


    //let sub_len = 9;
    //let mut cells = SampleIds::new( sub_len );// = Vec::with_capacity(12);
    //cells.init_rhapsody( &opts.specie );

    let cell_umi = if opts.cell_ids != "no file" {
        parse_bc_map( &opts.cell_ids )
    }   
    else{
        HashSet::new()
    };

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new(&opts.version, 7);

    // let mut gex = SingleCellData::new( );

    // let mut expr_file = parse_fastx_file(&opts.expression).expect("valid path/file");

    // let mut sample_names:Vec<String> = Vec::with_capacity(12);
    // let mut gene_names:Vec<String> = Vec::with_capacity(600);
    // let mut ab_names:Vec<String> = Vec::with_capacity(30);

    // for i in 1..13{
    //     sample_names.push( format!("Sample{}",i) )
    // }

    // if  opts.specie.eq("human") {
    //     // get all the human sample IDs into this.
    //     genes.add( b"ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG", "Sample1".to_string() );
    //     genes.add( b"TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG", "Sample2".to_string() );
    //     genes.add( b"CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT", "Sample3".to_string() );
    //     genes.add( b"ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT", "Sample4".to_string() );
    //     genes.add( b"CTCCCTGGTGTTCAATACCCGATGTGGTGGGCAGAATGTGGCTGG", "Sample5".to_string() );
    //     genes.add( b"TTACCCGCAGGAAGACGTATACCCCTCGTGCCAGGCGACCAATGC", "Sample6".to_string() );
    //     genes.add( b"TGTCTACGTCGGACCGCAAGAAGTGAGTCAGAGGCTGCACGCTGT", "Sample7".to_string() );
    //     genes.add( b"CCCCACCAGGTTGCTTTGTCGGACGAGCCCGCACAGCGCTAGGAT", "Sample8".to_string() );
    //     genes.add( b"GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG", "Sample9".to_string() );
    //     genes.add( b"GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC", "Sample10".to_string() );
    //     genes.add( b"CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC", "Sample11".to_string() );
    //     genes.add( b"GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG", "Sample12".to_string() );
    // }
    // else if opts.specie.eq("mouse") {
    //     // and the mouse ones
    //     genes.add( b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", "Sample1".to_string() );
    //     genes.add( b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", "Sample2".to_string() );
    //     genes.add( b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", "Sample3".to_string() );
    //     genes.add( b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT", "Sample4".to_string() );
    //     genes.add( b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", "Sample5".to_string() );
    //     genes.add( b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC", "Sample6".to_string() );
    //     genes.add( b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", "Sample7".to_string() );
    //     genes.add( b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG", "Sample8".to_string() );
    //     genes.add( b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", "Sample9".to_string() );
    //     genes.add( b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC", "Sample10".to_string() );
    //     genes.add( b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", "Sample11".to_string() );
    //     genes.add( b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC", "Sample12".to_string() );

    // } else {
    //     println!("Sorry, but I have no primers for species {}", &opts.specie);
    //     std::process::exit(1)
    // }

    // that is a class to strore gene expression data.
    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
    // and umi's are u64 again
    // here I need the cell kmer site.

    let mut ofile:Ofiles = Ofiles::new( 1, "cells", &opts.reads.to_string(), &opts.file.to_string(), &opts.outpath.to_string() );

    let mut ok_reads = 0;

    {

        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
        let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

        'main: while let Some(record2) = readefile.next() {
            if let Some(record1) = readereads.next() {

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
                // println!("The read1 mean quality == {}", mean_u8( seqrec.qual().unwrap() ));

                // if mean_u8( seqrec.qual().unwrap() ) < 50.0 {
                //     unknown +=1;
                //     println!("filtered a read");
                //     continue 'main;
                // }
                // if mean_u8( seqrec1.qual().unwrap() ) < 50.0 {
                //     unknown +=1;
                //     println!("filtered a read");
                //     continue 'main;
                // }

                // totally unusable sequence
                // as described in the BD rhapsody Bioinformatic Handbook
                // GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf (google should find this)
                if seqrec1.seq().len() < min_sizes[0] {
                    continue;
                }
                if seqrec.seq().len() < min_sizes[1] {
                    continue;
                }
                //let seq = seqrec.seq().into_owned();
                for nuc in &seqrec1.seq()[pos[6]..pos[7]] {  
                    if *nuc ==b'N'{
                        continue 'main;
                    }
                }

                //let umi = Kmer::from( &seqrec1.seq()[52..60]).into_u64();

                // first match the cell id - if that does not work the read is unusable
                //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                match cells.to_cellid( &seqrec1.seq(), vec![pos[0],pos[1]], vec![pos[2],pos[3]], vec![pos[4],pos[5]]){
                    Ok(cell_id) => {
                        // this is removing complexity from the data - in the test dataset 111 reads are ignored.
                        // let cell_id_umi:u128 = read_be_u128(  [ umi.to_be_bytes() , (cell_id as u64).to_be_bytes() ].concat().as_slice() );
                        // if ! cell_umi.insert( cell_id_umi ){
                        //     continue 'main;
                        // }

                        if opts.cell_ids != "no file" && ! cell_umi.contains( &cell_id ) {
                            continue 'main;
                        }

                        match seqrec1.write(&mut ofile.buff1, None){
                            Ok(_) => (),
                            Err(err) => println!("{err}")
                        };
                        match seqrec.write( &mut ofile.buff2, None){
                            Ok(_) => (),
                            Err(err) => println!("{err}")
                        };
                        ok_reads += 1;
                        if ok_reads == opts.max_reads{
                            break 'main;
                        }
                    },
                    Err(_err) => {
                        // cell id could not be recovered
                        
                        // println!("Cell ID could not be recovered from {:?}:\n{}\n{:?}, {:?}, {:?}", std::str::from_utf8(&seqrec1.seq()), _err, 
                        //     std::str::from_utf8( &seqrec1.seq()[pos[0]..pos[1]]),
                        //     std::str::from_utf8( &seqrec1.seq()[pos[2]..pos[3]]),
                        //     std::str::from_utf8( &seqrec1.seq()[pos[4]..pos[5]])
                        // );
                        continue
                    }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                };
            } else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }


        println!( "\nSummary:");
        println!(     "total      reads  : {ok_reads} reads" );

    }


    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            let sec = milli % 60;
            milli= (milli -sec) /60;

            let min = milli % 60;
            milli= (milli -min) /60;

            println!("quantify_rhapsody finished in {milli}h {min}min {sec} sec {mil}milli sec\n" );},
       Err(e) => {println!("Error: {e:?}");}
    }

}

use clap::Parser;
use needletail::parse_fastx_file;
use this::cellids::CellIds;
use this::ofiles::Ofiles;

use std::collections::HashSet;

//use this::sampleids::SampleIds;

//use this::last5::Last5;

use std::fs;

use std::time::SystemTime;


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
    #[clap(default_value_t=usize::MAX, long)]
    max_reads: usize,
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
    // match fs::create_dir_all(&opts.outpath).expect("AlreadyExists"){
    //     Ok(_) => (),
    //     Err(e) => {
    //         panic!("I could not create the outpath: {}", e);
    //     }
    // };

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

    // let mut cell_umi:HashSet<u128> = HashSet::new();

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new(&opts.version, 7);

    // that is a class to strore gene expression data.
    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
    // and umi's are u64 again
    // here I need the cell kmer site.

    let mut ofile:Ofiles = Ofiles::new( 1, &"cells".to_string(), &opts.reads.to_string(), &opts.file.to_string(), &opts.outpath.to_string() );

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
                        eprintln!("could not read from R2:\n{}", err);
                        continue 'main;
                    }
                };
                let seqrec1 = match record1{
                    Ok(res) => res,
                    Err(err) => {
                        eprintln!("could not read from R1:\n{}", err);
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
                    Ok(_cell_id) => {
                        // this is removing complexity from the data - in the test dataset 111 reads are ignored.
                        // let cell_id_umi:u128 = read_be_u128(  [ umi.to_be_bytes() , (cell_id as u64).to_be_bytes() ].concat().as_slice() );
                        // if ! cell_umi.insert( cell_id_umi ){
                        //     continue 'main;
                        // }

                        match seqrec1.write(&mut ofile.buff1, None){
                            Ok(_) => (),
                            Err(err) => println!("{}",err)
                        };
                        match seqrec.write( &mut ofile.buff2, None){
                            Ok(_) => (),
                            Err(err) => println!("{}",err)
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
        println!(     "total      reads  : {} reads", ok_reads );

    }


    match now.elapsed() {
        Ok(elapsed) => {
            let mut milli = elapsed.as_millis();
            let sec:u128;
            let min:u128;

            let mil = milli % 1000;
            milli= (milli - mil) /1000;

            sec = milli % 60;
            milli= (milli -sec) /60;

            min = milli % 60;
            milli= (milli -min) /60;

            println!("quantify_rhapsody finished in {}h {}min {} sec {}milli sec\n", milli, min, sec, mil );},
       Err(e) => {println!("Error: {e:?}");}
    }

}

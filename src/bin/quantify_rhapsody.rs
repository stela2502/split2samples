use clap::Parser;
use needletail::parse_fastx_file;
use kmers::naive_impl::Kmer;
use this::cellids::CellIds;

//use this::sampleids::SampleIds;
use this::singlecelldata::SingleCellData;
use this::geneids::GeneIds;

//use this::last5::Last5;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
//use std::time::Duration;


use std::io::BufWriter;
use std::fs::File;
use std::io::Write;

// use std::collections::HashSet;
// use std::convert::TryInto;

/// Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles.
/// You need quite long R1 and R2 reads for this! (>70R1 and >70R2 [v1] and 52 bp reads for v2.96 and v2.384)

#[derive(Parser)]
#[clap(version = "0.3.4", author = "Stefan L. <stefan.lang@med.lu.se>")]
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
    antibody: String,
    /// the minimum reads per cell (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,
    /// UMI min count - use every umi (per gene; 1) or only reoccuring ones (>1)
    #[clap(default_value_t=1,short, long)]
    umi_count: u8,
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

    

    //let sub_len = 9;
    //let mut cells = SampleIds::new( sub_len );// = Vec::with_capacity(12);
    //cells.init_rhapsody( &opts.specie );

    // let mut cell_umi:HashSet<u128> = HashSet::new();
    let mut genes:GeneIds = GeneIds::new(12); // split them into 9 bp kmers

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new(&opts.version, 7);

    // that is a class to strore gene expression data.
    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
    // and umi's are u64 again
    // here I need the cell kmer site.
    let mut gex = SingleCellData::new( );

    let mut expr_file = parse_fastx_file(&opts.expression).expect("valid path/file");

    let mut sample_names:Vec<String> = Vec::with_capacity(12);
    let mut gene_names:Vec<String> = Vec::with_capacity(600);
    let mut ab_names:Vec<String> = Vec::with_capacity(30);

    for i in 1..13{
        sample_names.push( format!("Sample{}",i) )
    }

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
    

    let mut i = 1;

    while let Some(e_record) = expr_file.next() {
        let seqrec = e_record.expect("invalid record");
        match std::str::from_utf8(seqrec.id()){
            Ok(st) => {
                for id in st.to_string().split("|"){
                    i += 1;
                    //println!("Adding gene #{} {}",i, id );
                    genes.add( &seqrec.seq(), id.to_string() );
                    gene_names.push( id.to_string() );
                    break;
                }
            },
            Err(err) => eprintln!("The expression entry's id could not be read: {}", err),
        }
    }

    let mut ab_file = parse_fastx_file(&opts.antibody).expect("valid path/file");
    while let Some(ab_record) = ab_file.next() {
        let seqrec = ab_record.expect("invalid record");
        match std::str::from_utf8(seqrec.id()){
            Ok(st) => {
                for id in st.to_string().split("|"){
                    genes.add( &seqrec.seq(), id.to_string() );
                    ab_names.push( id.to_string() );
                    break;
                }
            },
            Err(err) => eprintln!("The expression entry's id could not be read: {}", err),
        }
    }

    let mut unknown = 0;
    let mut no_sample = 0;
    let mut ok_reads = 0;
    let mut pcr_duplicates = 0;
    let mut local_dup = 0;
    let split:usize = 1000*1000;
    let mut log_iter = 0;
    //let split:usize = 1000;

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

    let log_file_str = PathBuf::from(&opts.outpath).join(
        format!("Mapping_log.txt" )
    );

    println!( "the log file: {}", log_file_str.file_name().unwrap().to_str().unwrap() );
    
    let log_file = match File::create( log_file_str ){
        Ok(file) => file,
        Err(err) => {
            panic!("Error: {:#?}", err);
        }
    };
    let mut log_writer = BufWriter::new(&log_file);

    {
        let spinner_style = ProgressStyle::with_template("{prefix:.bold.dim} {spinner} {wide_msg}")
            .unwrap()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
        let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");
        let m = MultiProgress::new();
        let pb = m.add(ProgressBar::new(5000));
        pb.set_style(spinner_style.clone());
        pb.set_prefix(format!("[{}/?]", i + 1));

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

                // totally unusable sequence
                // as described in the BD rhapsody Bioinformatic Handbook
                // GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf (google should find this)
                if seqrec1.seq().len() < min_sizes[0] {
                    unknown +=1;
                    continue;
                }
                if seqrec.seq().len() < min_sizes[1] {
                    unknown +=1;
                    continue;
                }
                //let seq = seqrec.seq().into_owned();
                for nuc in &seqrec1.seq()[pos[6]..pos[7]] {  
                    if *nuc ==b'N'{
                        unknown +=1;
                        continue 'main;
                    }
                }
                //let umi = Kmer::from( &seqrec1.seq()[52..60]).into_u64();
                let umi = Kmer::from( &seqrec1.seq()[pos[6]..pos[7]]).into_u64();

                // first match the cell id - if that does not work the read is unusable
                //match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                match cells.to_cellid( &seqrec1.seq(), vec![pos[0],pos[1]], vec![pos[2],pos[3]], vec![pos[4],pos[5]]){
                    Ok(cell_id) => {
                        // this is removing complexity from the data - in the test dataset 111 reads are ignored.
                        // let cell_id_umi:u128 = read_be_u128(  [ umi.to_be_bytes() , (cell_id as u64).to_be_bytes() ].concat().as_slice() );
                        // if ! cell_umi.insert( cell_id_umi ){
                        //     continue 'main;
                        // }
                        match genes.get( &seqrec.seq() ){
                            Some(gene_id) =>{
                                match gex.get( cell_id as u64, format!( "Cell{}", cell_id ) ){
                                    Ok(cell_info) => {
                                        //println!("I got gene {} for cell {}", gene_id , cell_info.name);
                                        ok_reads += 1;
                                        if cell_info.add( gene_id, umi ) == false {
                                            //println!("  -> pcr duplicate");
                                            pcr_duplicates += 1;
                                            local_dup += 1;
                                        }
                                        if ok_reads == opts.max_reads{
                                            break 'main;
                                        }
                                        if ok_reads % split == 0{
                                            log_iter+=1;

                                            let log_str = format!("{} mio usable ({:.2}% total; {:.2}% PCR dupl. [{:.2}% for last batch]): {:?} -> gene id: {}, cell id: {}",
                                                log_iter , 
                                                ok_reads as f32 / (ok_reads +no_sample+ unknown) as f32 * 100.0 , 
                                                pcr_duplicates as f32 / ok_reads as f32 * 100.0,
                                                local_dup as f32 / split as f32 * 100.0,
                                                std::str::from_utf8(&seqrec1.seq()), gene_id , cell_id 
                                            );
                                            pb.set_message( log_str.clone() );
                                            pb.inc(1);
                                            match writeln!( log_writer, "{}", log_str ){
                                                Ok(_) => (),
                                                Err(err) => {
                                                    eprintln!("write error: {}", err);
                                                }
                                            };
                                            local_dup = 0;
                                            //std::thread::sleep(Duration::from_millis(100));
                                        }
                                    },
                                    Err(err) => panic!("Could not add a gene expression: gene_id = {}, umi = {} and err= {}",gene_id, Kmer::from( &seqrec1.seq()[52..60]).into_u64(),  err),
                                }
                                //println!("Cool I got a gene id: {}", id);
                            },
                            None => {
                                unknown +=1;
                                // all - samples genes and antibodies are classed as genes here.
                                //eprintln!("I could not identify a gene in this read: {:?}", std::str::from_utf8( &seqrec.seq() ) );
                            }
                        };

                    },
                    Err(_err) => {
                        // cell id could not be recovered
                        
                        // println!("Cell ID could not be recovered from {:?}:\n{}\n{:?}, {:?}, {:?}", std::str::from_utf8(&seqrec1.seq()), _err, 
                        //     std::str::from_utf8( &seqrec1.seq()[pos[0]..pos[1]]),
                        //     std::str::from_utf8( &seqrec1.seq()[pos[2]..pos[3]]),
                        //     std::str::from_utf8( &seqrec1.seq()[pos[4]..pos[5]])
                        // );
                        
                        no_sample +=1;
                        continue
                    }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                };
            } else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        let log_str = format!("{} mio usable ({:.2}% total; {:.2}% PCR dupl. [{:.2}% for last batch])",
            log_iter , 
            ok_reads as f32 / (ok_reads +no_sample+ unknown) as f32 * 100.0 , 
            pcr_duplicates as f32 / ok_reads as f32 * 100.0,
            local_dup as f32 / split as f32 * 100.0
        );

        pb.finish_with_message( log_str );
        println!("\n\nWriting outfiles ...");

        // calculating a little bit wrong - why? no idea...
        //println!( "collected sample info:i {}", gex.mtx_counts( &mut genes, &gene_names, opts.min_umi , opts.umi_count) );


        //let fp1 = PathBuf::from(opts.reads.clone());
        //println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
        let file_path = PathBuf::from(&opts.outpath).join(
            format!("SampleCounts.tsv" )
        );

        let file_path_sp = PathBuf::from(&opts.outpath).join(
            format!("BD_Rhapsody_expression" )
        );

        // this always first as this will decide which cells are OK ones!
        match gex.write_sparse_sub ( file_path_sp, &mut genes , &gene_names, opts.min_umi, opts.umi_count ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

        let file_path_sp = PathBuf::from(&opts.outpath).join(
            format!("BD_Rhapsody_antibodies" )
        );

        match gex.write_sparse_sub ( file_path_sp, &mut genes, &ab_names, 1, 1 ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

    
        match gex.write_sub ( file_path, &mut genes, &sample_names, 0, 0 ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

        

        let total = no_sample+ unknown + ok_reads;
        let reads_genes = gex.n_reads( &mut genes , &gene_names );
        let reads_ab = gex.n_reads( &mut genes , &ab_names );
        let reads_samples = gex.n_reads( &mut genes , &sample_names );

        println!( "\nSummary:");
        println!(     "total      reads  : {} reads", total );
        println!(     "no cell ID reads  : {} reads", no_sample );
        println!(     "N's or too short  : {} reads", unknown );
        println!(     "usable reads      : {} reads ({:.2}% of total)", ok_reads, (ok_reads as f32 / total as f32) * 100.0 );
        println!(     "expression reads  : {} reads ({:.2}% of total)", reads_genes, (reads_genes as f32 / total as f32) * 100.0 );
        println!(     "antibody reads    : {} reads ({:.2}% of total)", reads_ab, (reads_ab as f32 / total as f32) * 100.0 );
        println!(     "sample tag reads  : {} reads ({:.2}% of total)", reads_samples, (reads_samples as f32 / total as f32) * 100.0 );
        println!(     "pcr duplicates    : {} reads ({:.2}% of usable)", pcr_duplicates, ( pcr_duplicates as f32 / ok_reads as f32 ) * 100.0 );
        

        let file_path2 = format!("{}/SampleCounts.tsv", opts.outpath );
        println!( "\nCell->Sample table written to {:?}\n", file_path2);

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

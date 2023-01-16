use clap::Parser;
use needletail::parse_fastx_file;
use kmers::naive_impl::Kmer;
use this::cellids::CellIds;

use this::sampleids::SampleIds;
use this::singlecelldata::SingleCellData;
use this::geneids::GeneIds;

use std::path::PathBuf;
use std::fs;

use std::time::SystemTime;

/// Quantifies a DB Rhapsody experiment and creates sparse matrix outfiles.
/// You need quite long R1 and R2 reads for this! (>70R1 and >70R2)

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
    antibody: String,
    /// the minimum reads (sample + genes + antibody combined)
    #[clap(short, long)]
    min_umi: usize,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    version: String,
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
    let mut cells = SampleIds::new( sub_len );// = Vec::with_capacity(12);
    cells.init_rhapsody( &opts.specie );


    let mut genes:GeneIds = GeneIds::new(9); // split them into 9 bp kmers

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

    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new(&opts.version);

    // that is a class to strore gene expression data.
    // sample ids are meant to be u64, gene ids usize (as in the GeneIds package)
    // and umi's are u64 again
    let mut gex = SingleCellData::new( 9 );

    let mut unknown = 0;
    let mut no_sample = 0;
    let mut ok_reads = 0;

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
                        match genes.get( &seqrec.seq() ){
                            Some(gene_id) =>{
                                match gex.get( cell_id as u64, format!( "Cell{}", cell_id ) ){
                                    Ok(cell_info) => {
                                        ok_reads += 1;
                                        cell_info.add( gene_id, umi );
                                    },
                                    Err(err) => panic!("Could not add a gene expression: gene_id = {}, umi = {} and err= {}",gene_id, Kmer::from( &seqrec1.seq()[52..60]).into_u64(),  err),
                                }
                                //println!("Cool I got a gene id: {}", id);
                            },
                            None => {
                                unknown +=1;
                                //eprintln!("I could not identify a gene in this read: {:?}", std::str::from_utf8( &seqrec.seq() ) );
                            }
                        };

                    },
                    Err(_err) => {
                        // cell id could not be recovered
                        //println!("Cell ID could not be recovered from {:?}:\n{}", std::str::from_utf8(&seqrec1.seq()), _err);
                        no_sample +=1;
                        continue
                    }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                };
            } else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        // calculating a little bit wrong - why? no idea...
        println!( "collected sample info:i {}", gex.mtx_counts( &mut genes, &gene_names, opts.min_umi ) );


        let fp1 = PathBuf::from(opts.reads.clone());
        println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
        let file_path = PathBuf::from(&opts.outpath).join(
            format!("SampleCounts.tsv" )
        );

        let file_path_sp = PathBuf::from(&opts.outpath).join(
            format!("BD_Rhapsody_expression" )
        );
        // this always first as this will decide which cells are OK ones!
        match gex.write_sparse_sub ( file_path_sp, &mut genes , &gene_names, opts.min_umi ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

        let file_path_sp = PathBuf::from(&opts.outpath).join(
            format!("BD_Rhapsody_antibodies" )
        );
        // this always first as this will decide which cells are OK ones!
        
        match gex.write_sparse_sub ( file_path_sp, &mut genes, &ab_names, opts.min_umi / 2 ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

    
        match gex.write_sub ( file_path, &mut genes, &sample_names, 0 ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };

        

        let total = no_sample+ unknown + ok_reads;
        println!(     "no sample ID reads: {} reads", no_sample );
        println!(     "N's or too short  : {}", unknown );
        println!(     "usable reads      : {} ({:.2}%)", ok_reads, (ok_reads as f32 / total as f32) );

        let file_path2 = format!("{}/SampleCounts.tsv", opts.outpath );
        println!( "written to {:?}", file_path2);
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

            println!("quantifyRhapsody finished in {}h {}min {} sec {}milli sec", milli, min, sec, mil );},
       Err(e) => {println!("Error: {e:?}");}
    }

}

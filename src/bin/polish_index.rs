use clap::Parser;

use rustody::fast_mapper::FastMapper;
use rustody::mapping_info::MappingInfo;


/// The indices from create_index_te might have genes in the genes.txt file that are actually not covered by the index.
/// This is no problem per se, but might increase the size of the index significantly.
/// This tool removes uncovered genes from the index.


#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the gff/gtf gene information table (ONE gzipped stream - do NOT cat these files!)
    #[clap(short, long)]
    index: String,
    /// the outpath
    #[clap(default_value=  "testData/mapperTest",short, long)]
    outpath: String,
}


fn main() {
    // parse the options
    let mut report = MappingInfo::new(None, 32.0 , 0, None );
    report.start_counter();

    
    //// create the report object /////////////////////////////////////
    let opts: Opts = Opts::parse();

    let mut index = FastMapper::new( 32, 10 );

    match index.load_index( opts.index ){
    	Ok(_) => {},
    	Err(e) => {
    		panic!("polish_index main: Load index hit an error: {e:?}" )
    	}
    }

    report.stop_file_io_time();

    index.print();

    let mut fixed_index = FastMapper::new(32, 10);

    fixed_index.merge(index); // this should re-calculate the genes vector!

    report.stop_single_processor_time();

    match fixed_index.write_index( opts.outpath ){
    	Ok(_) => {},
    	Err(e) => {
    		panic!("polish_index main: Write index hit an error: {e:?}" )
    	}
    }

    report.stop_file_io_time();

    println!("polish_index has finished:\n{}", report.program_states_string() );

}





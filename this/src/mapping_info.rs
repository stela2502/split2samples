use crate::ofiles::Ofiles;
use std::fs::File;
use indicatif::{ProgressBar};

use std::io::Write;

/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct MappingInfo{
	/// reads that did not pass the filters
	pub unknown:usize,
	/// reads that had no cell id
    pub no_sample:usize,
    /// reads that have no match in the geneIds object
    pub no_data:usize,
    /// reads with cell id and gene id
    pub ok_reads:usize,
    /// reads that are duplicates on the UMI level per cell and gene
    pub pcr_duplicates:usize,
    /// the amount of ok_reads after which to write a entry into the log file
   	pub split:usize,
   	/// the others are explained in the quantify_rhapsody.rs file.
    log_iter:usize,
    pub log_writer:File,
    pub min_quality:f32, 
    pub max_reads:usize, 
    pub ofile:Ofiles,
    pub local_dup:usize,
    pub total:usize,
}

impl MappingInfo{
	pub fn new(log_writer:File, min_quality:f32, max_reads:usize, ofile:Ofiles, ) -> Self{
		let unknown = 0;
		let no_sample = 0;
		let no_data = 0;
		let ok_reads = 0;
		let pcr_duplicates = 0;
		let split = 1_000_000;
		let log_iter = 0;
		let local_dup = 0;
		let total = 0;
		Self{
			unknown,
			no_sample,
			no_data,
			ok_reads,
			pcr_duplicates,
			split,
			log_iter,
			log_writer,
			min_quality,
			max_reads,
			ofile,
			local_dup,
			total,
		}
	}

	pub fn write_to_log ( &mut self, text:String ){
		match writeln!( self.log_writer, "{text}" ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {err}" );
            }
        };
	}
	pub fn log( &mut self, pb:&ProgressBar ){
		if self.total % self.split == 0{
			self.log_iter+=1;
            let log_str = self.log_str();
            pb.set_message( log_str.clone() );
            pb.inc(1);
            match writeln!( self.log_writer, "{log_str}" ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}" );
                }
            };
            self.local_dup = 0;
            //std::thread::sleep(Duration::from_millis(100));
		}
	}
	pub fn log_str( &mut self ) -> String{
		format!("{:.2} mio reads ({:.2}% usable; {:.2}% PCR dupl. [usable] [{:.2}% for last batch])",
            self.total as f32 / self.split as f32,
            self.ok_reads as f32 / (self.ok_reads +self.no_sample+ self.unknown) as f32 * 100.0 , 
            self.pcr_duplicates as f32 / self.ok_reads as f32 * 100.0,
            self.local_dup as f32 / self.split as f32 * 100.0
         )
	}
}

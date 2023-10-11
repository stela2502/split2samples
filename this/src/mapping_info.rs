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
    /// reads with cell_id - gene_id is not checked
    pub cellular_reads: usize,
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
		let cellular_reads = 0;
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
			cellular_reads,
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


	pub fn merge(&mut self, other:&MappingInfo ){
		self.no_sample += other.no_sample;
		self.no_data += other.no_data;
		//unknown is defined without multiprocessor support
		self.ok_reads += other.ok_reads;
		self.pcr_duplicates += other.pcr_duplicates;
		self.cellular_reads += other.cellular_reads;
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
		format!("{:.2} mio reads ({:.2}% with cell info, {:.2}% with gene match",
            self.total as f32 / self.split as f32,
            self.cellular_reads as f32 / (self.total) as f32 * 100.0 , 
            self.ok_reads as f32 / (self.total) as f32 * 100.0 
         )
	}

	pub fn summary( &mut self, reads_genes:usize, reads_ab:usize, reads_samples:usize ) -> String{
	    let result = "\nSummary:\n".to_owned()
	    	+format!(     "total      reads  : {} reads\n", self.total ).as_str()
	    	+format!(     "no cell ID reads  : {} reads ({:.2}% of total)\n", self.no_sample, (self.no_sample as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "no gene ID reads  : {} reads ({:.2}% of total)\n", self.no_data, (self.no_data as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "N's or too short  : {} reads ({:.2}% of total)\n", self.unknown, (self.unknown as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "cellular reads    : {} reads ({:.2}% of total)\n", self.cellular_reads, (self.cellular_reads as f32 / self.total as f32) * 100.0 ).as_str()
	    	+format!(     "expression reads  : {} reads ({:.2}% of cellular)\n", reads_genes, (reads_genes as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "antibody reads    : {} reads ({:.2}% of cellular)\n", reads_ab, (reads_ab as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "sample   reads    : {} reads ({:.2}% of cellular)\n", reads_samples, (reads_samples as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "mapped reads      : {} reads ({:.2}% of cellular)\n", reads_genes + reads_ab + reads_samples, ( (reads_genes + reads_ab + reads_samples) as f32 / self.cellular_reads as f32) * 100.0 ).as_str();
	   		//+format!(     "pcr duplicates    : {} reads ({:.2}% of cellular)\n\n", results.pcr_duplicates, ( self.pcr_duplicates as f32 / self.ok_reads as f32 ) * 100.0 ).as_str();
	   	match writeln!( self.log_writer, "{result}" ){
                Ok(_) => (),
                Err(err) => {
                    eprintln!("write error: {err}" );
                }
            };
        return result
	    
	}
}

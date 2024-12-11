use crate::ofiles::{Ofiles, Fspot};
use std::fs::File;
use indicatif::{ProgressBar};
use std::collections::BTreeMap;

use std::io::Write;
use std::time::{Duration, SystemTime};

use chrono::{DateTime, Utc};
use std::collections::HashMap;

/// MappingInfo captures all mapping data and is a way to easily copy this data over multiple analysis runs.
pub struct MappingInfo{
	/// reads that did not pass the filters
	pub quality:usize,
	pub length:usize,
	analyzed:usize,
	pub n_s:usize,
	pub poly_a:usize,
	/// reads that had no cell id
    pub no_sample:usize,
    /// reads that have no match in the geneIds object
    pub no_data:usize,
    /// reads with cell id and gene id
    pub ok_reads:usize,
    /// reads with cell_id - gene_id is not checked
    pub cellular_reads: usize,
    pub multimapper: usize,
    /// reads that are duplicates on the UMI level per cell and gene
    pub pcr_duplicates:usize,
    /// the amount of ok_reads after which to write a entry into the log file
   	pub split:usize,
   	/// the others are explained in the quantify_rhapsody.rs file.
    log_iter:usize,
    pub log_writer: Option<File>,
    pub min_quality:f32, 
    pub max_reads:usize, 
    pub ofile: Option<Ofiles>,
    pub local_dup:usize,
    pub total:usize,
    pub absolute_start: SystemTime,
    realtive_start: Option<SystemTime>,
    tmp_counter: Option<SystemTime>,
    pub single_processor_time: Duration,
    pub multi_processor_time: Duration,
    pub file_io_time: Duration,
    pub subprocess_time: Duration,
    pub reads_log: BTreeMap<String, usize >,
    pub error_counts: HashMap<String, usize>,  // To store error types and their counts

}

impl MappingInfo{
	pub fn new(log_writer:Option<File>, min_quality:f32, max_reads:usize, ofile:Option<Ofiles>, ) -> Self{
		let absolute_start = SystemTime::now();
		let single_processor_time = Duration::new(0,0);
		let multi_processor_time = Duration::new(0,0);
		let file_io_time = Duration::new(0,0);
		let reads_log = BTreeMap::new();
		let subprocess_time = Duration::new(0,0);
		let mut this = Self{
			quality: 0,
		    length: 0,
		    analyzed: 1,
		    n_s: 0,
		    poly_a: 0,
			no_sample: 0,
			no_data: 0,
			ok_reads: 0,
			cellular_reads: 0,
			multimapper: 0,
			pcr_duplicates: 0,
			split: 1_000_000,
			log_iter: 0,
			log_writer,
			min_quality,
			max_reads,
			ofile,
			local_dup: 0,
			total: 0,
			absolute_start,
			realtive_start: None,
			tmp_counter: None,
			single_processor_time,
			multi_processor_time,
			file_io_time,
			subprocess_time,
			reads_log,
			error_counts: HashMap::new(),  // Initialize the HashMap
		};
		this.start_counter();
		this
	}

	pub fn start_counter ( &mut self ){
		self.realtive_start = Some( SystemTime::now() );
	}

	pub fn start_ticker ( &mut self )  {
		self.tmp_counter = Some( SystemTime::now() );
	}
	pub fn stop_ticker ( &mut self ) -> ( u128, u128, u128, u128 ) {
		let ret = MappingInfo::split_duration( self.tmp_counter.unwrap_or( SystemTime::now() ).elapsed().unwrap() );
		self.tmp_counter = Some( SystemTime::now() );
		ret
	}

	pub fn stop_single_processor_time ( &mut self ) {
		self.single_processor_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn stop_multi_processor_time ( &mut self ) {
		self.multi_processor_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn subprocess_time ( &mut self ) {
		self.subprocess_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn stop_file_io_time ( &mut self ) {
		self.file_io_time += self.realtive_start.unwrap().elapsed().unwrap();
		self.start_counter();
	}

	pub fn elapsed_time_split ( &self ) -> ( u128, u128, u128, u128 ){
		let elapsed = self.absolute_start.elapsed().unwrap();
		MappingInfo::split_duration( elapsed )
	}

	pub fn now (&self) -> String{

		let now: DateTime<Utc> = Utc::now();
		format!("{}", now)
	    
	}

	// Unified reporting method that logs errors into the HashMap
    pub fn report(&mut self, issue: &str) {
        // Increment the count for the issue type
        *self.error_counts.entry(issue.to_string()).or_insert(0) += 1;
        //println!("Issue reported: {}", issue);
    }

    // Optionally, add a method to retrieve counts for a specific issue
    pub fn get_issue_count(&self, issue: &str) -> usize {
        *self.error_counts.get(issue).unwrap_or(&0)  // Return count or 0 if not present
    }

    // Method to export error_counts to a CSV file
    pub fn report_to_csv(&self, file_path: &str) {
        let mut file = File::create(file_path).unwrap();  // Create a file for writing
        writeln!(file, "Error Type,Count").unwrap();  // Write CSV header

        // Iterate over the error_counts and write each as a row in the CSV
        for (error_type, count) in &self.error_counts {
            writeln!(file, "{},{}", error_type, count);  // Write each error type and count
        }

        //Ok(())  // Return Ok if successful
    }

    // Method to export error_counts to a CSV-formatted String
	pub fn report_to_string(&self) -> String {
	    // Start with the header
	    let mut output = String::from("Error Type,Count\n");

	    // Iterate over the error_counts and append each as a row in the CSV format
	    for (error_type, count) in &self.error_counts {
	        // Append each error type and count, followed by a newline
	        output.push_str(&format!("{},{}\n", error_type, count));
	    }

	    output // Return the content as a String
	}

	pub fn split_duration( elapsed:Duration ) -> ( u128, u128, u128, u128 ){

        let mut milli = elapsed.as_millis();

        let mil = milli % 1000;
        milli= (milli - mil) /1000;

        let sec = milli % 60;
        milli= (milli -sec) /60;

        let min = milli % 60;
        milli= (milli -min) /60;

        (milli, min, sec, mil )

    }

    pub fn iter_read_type(&mut self, name:&str ){
    	*self.reads_log.entry(name.to_string()).or_insert(0) += 1; 
    }

    fn read_types_to_string(&self, names:Vec<&str> ) -> String {
        let mut formatted_entries = String::new();

        for name in &names {
            let value = self.reads_log.get(*name).unwrap_or(&0);
        	let formatted_name = format!("{:<18}:", name); 
            formatted_entries.push_str(&format!("{} {} reads ({:.2}% of cellular)\n", formatted_name, value, *value as f32 / self.cellular_reads as f32 *100_f32 ));
        }
        formatted_entries
    }


	pub fn merge(&mut self, other:&MappingInfo ){
		self.no_sample += other.no_sample;
		self.no_data += other.no_data;
		//unknown is defined without multiprocessor support
		self.ok_reads += other.ok_reads;
		self.pcr_duplicates += other.pcr_duplicates;
		self.cellular_reads += other.cellular_reads;
		for (name, value) in &other.reads_log {
			*self.reads_log.entry(name.to_string()).or_insert(0) += value;
		}
		self.analyzed = self.total;
		for (error_type, count) in &other.error_counts {
            // For each error type in `other`, increment the value in `self`
            *self.error_counts.entry(error_type.clone()).or_insert(0) += count;
        }
	}



	pub fn write_to_log ( &mut self, text:String ){

		match &mut self.log_writer{
			Some(file) => {
				match writeln!( file , "{text}" ){
		            Ok(_) => (),
		            Err(err) => {
		                eprintln!("write error: {err}" );
		            }
		        };
			},
			None => {},
		}
		
	}

	pub fn write_to_ofile ( &mut self, which: Fspot, text:String ){

		match &mut self.ofile{
			Some(file) => {
				file.write_to_oufile( which, text );
			},
			None => {},
		}
		
	}

	pub fn log_report( &mut self ) {
		let log_str = self.report_to_string();
		self.write_to_log( log_str );
	}

	pub fn log( &mut self, pb:&ProgressBar ){
		if self.total % self.split == 0{
			self.log_iter+=1;
            let log_str = self.log_str();
            pb.set_message( log_str.clone() );
            pb.inc(1);
            self.write_to_log( log_str );
            self.local_dup = 0;
		}
	}

	pub fn log_str( &mut self ) -> String{
		format!("{:.2} mio reads ({:.2}% with cell_id, {:.2}% with gene_id {:.2}% multimapper)",
            self.total as f32 / self.split as f32,
            self.cellular_reads as f32 / (self.analyzed) as f32 * 100.0 , 
            self.ok_reads as f32 / (self.analyzed) as f32 * 100.0,
            self.multimapper as f32 / (self.analyzed) as f32 * 100.0,
         )
	}
	pub fn program_states_string( &self ) -> String{
		let mut result = String::from("");
		let  (mut hours,mut min,mut sec ,mut mulli ) = Self::split_duration( self.absolute_start.elapsed().unwrap() );
	   	result += format!("   overall run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	( hours, min, sec , mulli ) = Self::split_duration( self.file_io_time);
	   	result += format!("   file-io run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	( hours, min, sec , mulli ) = Self::split_duration( self.single_processor_time);
	   	result += format!("single-cpu run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	( hours, min, sec , mulli ) = Self::split_duration( self.multi_processor_time);
	   	result += format!(" multi-cpu run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	if self.subprocess_time != Duration::new(0,0) {
	   		( hours, min, sec , mulli ) = Self::split_duration( self.subprocess_time);
	    	result += format!("subprocess run time {} h {} min {} sec {} millisec\n", hours, min, sec , mulli ).as_str();
	   	}
	   	result
	}

	pub fn summary( &mut self, reads_genes:usize, reads_ab :usize, reads_samples:usize ) -> String{

		let pcr_duplicates = self.cellular_reads - reads_genes - reads_ab - reads_samples;

		let unknown = self.quality + self.length + self.n_s + self.poly_a;
	    let mut result = "\nSummary:\n".to_owned()
	    	+format!(     "cellular   reads  : {} reads ({:.2}% of total)\n", self.cellular_reads, (self.cellular_reads as f32 / self.total as f32) * 100.0 ).as_str()
	    	+format!(     "no cell ID reads  : {} reads ({:.2}% of total)\n", self.no_sample, (self.no_sample as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "no gene ID reads  : {} reads ({:.2}% of total)\n", self.no_data.saturating_sub(self.no_sample), ( self.no_data.saturating_sub( self.no_sample) as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     "filtered   reads  : {} reads ({:.2}% of total)\n", unknown, (unknown as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " ->  multimapper  : {} reads ({:.2}% of total)\n", self.multimapper, ( self.multimapper as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " -> bad qualiity  : {} reads ({:.2}% of total)\n", self.quality, ( self.quality as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " ->    too short  : {} reads ({:.2}% of total)\n", self.length, ( self.length as f32 / self.total as f32) * 100.0).as_str()
	    	+format!(     " ->          N's  : {} reads ({:.2}% of total)\n", self.n_s, ( self.n_s as f32 / self.total as f32) * 100.0).as_str()
	    	+"\n"
	    	+format!(     "total      reads  : {} reads\n", self.total ).as_str()
	    	+"\ncollected read counts:\n"
	    	+self.read_types_to_string(vec!["expression reads", "antibody reads", "sample reads"]).as_str()
	    	+"\nreported UMI counts:\n"
	    	+format!(     "expression reads  : {} UMIs ({:.2}% of cellular)\n", reads_genes, (reads_genes as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "antibody reads    : {} UMIs ({:.2}% of cellular)\n", reads_ab, (reads_ab as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "sample reads      : {} UMIs ({:.2}% of cellular)\n", reads_samples, (reads_samples as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	//+format!(     "unique reads      : {} reads ({:.2}% of cellular)\n", reads_genes + reads_ab + reads_samples, ( (reads_genes + reads_ab + reads_samples) as f32 / self.cellular_reads as f32) * 100.0 ).as_str()
	    	+format!(     "\nPCR duplicates or bad cells: {} reads ({:.2}% of cellular)\n\n", pcr_duplicates, ( pcr_duplicates as f32 / self.cellular_reads as f32 ) * 100.0 ).as_str()
	   		+"timings:\n";
	   	result += &self.program_states_string();
	   	self.write_to_log( result.clone() );
        result
	}
	
}

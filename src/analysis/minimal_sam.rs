
use crate::genes_mapper::{ GenesMapper, MapperResult, SeqRec };

#[derive(Debug, Clone)]
pub struct MinimalSam {

}

impl MinimalSam {

	pub fn new() -> Self{
		Self{}
	}
	pub fn to_sam_line (&self, read:&SeqRec, gene_id:&Vec<MapperResult>, cell_id:&SeqRec, umi:&SeqRec, index:&GenesMapper ) -> Option<String>{

		let mut read2 = read.clone();
    	let ( start, cigar_str ) = match gene_id[0].cigar() {
    		Some(cigar) => 
    		{

    			let (mine, _other) = cigar.calculate_covered_nucleotides( &cigar.to_sam_string() );
    			if mine < read2.len(){
    				#[cfg(debug_assertions)]
    				eprintln!("build_sam_record tries to slice the {}bp long seq \n{read2}\n as the cigar {cigar} is shorter.", read2.len());
    				read2=read2.slice(0, mine ).unwrap();

    			}
    			if mine + gene_id[0].start() > gene_id[0].db_length() {
    				panic!("Cigar suggest longer match than db_length allows!");
    			}
    			let dels_at_start = cigar.length_del_start();
    			if dels_at_start > 0 {
    				( 
    					gene_id[0].start() +1 + dels_at_start ,
    					cigar.to_sam_string()
    				)
    			}else {
    				(
    					gene_id[0].start() +1,
    					cigar.to_sam_string()
    				)
    			}

    		},
    		None=> {
    			panic!("An alignement needs to always have a cigar attached!")
    		},
    	};


    	let mut record = "".to_string();
    	// starting
    	record += &String::from_utf8_lossy(read2.id());
    	record += "\t";
    	// the ID
    	record += &format!("{}\t", 0b00000000000);
    	/*
            0x1 (1): PAIRED - Read is paired in sequencing.
		    0x2 (2): PROPER_PAIR - Read is mapped in a proper pair.
		    0x4 (4): UNMAPPED - Read is unmapped.
		    0x8 (8): MATE_UNMAPPED - Mate is unmapped.
		    0x10 (16): REVERSE_STRAND - Read is mapped to the reverse strand.
		    0x20 (32): MATE_REVERSE_STRAND - Mate is mapped to the reverse strand.
		    0x40 (64): READ1 - Read is the first read in a pair.
		    0x80 (128): READ2 - Read is the second read in a pair.
		    0x100 (256): SECONDARY - Secondary alignment.
		    0x200 (512): QC_FAIL - Read fails quality checks.
		    0x400 (1024): DUPLICATE - PCR or optical duplicate.
		 */
		let gene_name = match index.get_gene(gene_id[0].gene_id()){
		 	Some(gene_data) => {
		 		gene_data.get_unique_name()
		 	},
		 	None => {
		 		return None
		 	}
		};
		record += &format!("{}\t", gene_name ); // needs to be the transcript id!
		// the name of the database sequence
    	record += &format!("{}\t", start );
    	// the start position of the read
    	record += &format!("{}\t", gene_id[0].mapq() );
    	// the map quality
		record += &cigar_str;
    	record +="\t";
		// the Cigar string
		record += "*\t";
		// the sequence name of the mate/next read
		record += "0\t";
		//1-based position of the mate/next read
		record += "0\t";
		//Template length
		record += &String::from_utf8_lossy( &read2.seq() ) ;
    	record +="\t";
		//Read sequence
		record += &String::from_utf8_lossy( &read2.qual() ) ;
    	record +="\t";
		//Phred quality scores

		// and here come the 10x specific tags - might be needed later on
		record += &format!("NH:i:{}\t", gene_id.len());
	    /*NH:i:9: Number of reported alignments.
	        NH: Tag identifier, stands for "Number of reported alignments."
	        i: Tag data type, denotes a 32-bit signed integer.
	        9: Value associated with the tag, indicating the number of reported alignments for the read.*/
	    record += "HI:i:1\t";
		/*HI:i:2: Alignment hit index.
	        HI: Tag identifier, stands for "Alignment hit index."
	        i: Tag data type, denotes a 32-bit signed integer.
	        2: Value associated with the tag, representing the alignment hit index for the read. This index is typically used to distinguish between multiple hits for the same query in the SAM/BAM file.
		*/
	    record += &format!("AS:i:{}\t", gene_id[0].mapq() ); // AS:i:89
	    /*AS:i:89: Alignment score.
	        AS: Tag identifier, stands for "Alignment score."
	        i: Tag data type, denotes a 32-bit signed integer.
	        89: Value associated with the tag, indicating the alignment score for the read. The alignment score represents the quality or confidence of the alignment.
		*/
	    record += &format!("nM:i:{}\t", gene_id[0].edit_dist() ); // nM:i:0
	    /*nM:i:0: Edit distance.
	        nM: Tag identifier, stands for "Edit distance."
	        i: Tag data type, denotes a 32-bit signed integer.
	        0: Value associated with the tag, representing the edit distance between the read and the reference sequence. The edit distance is the number of mismatches and gaps (insertions and deletions) between the read and the reference sequence.
		*/
	    record += "RE:A:I\t"; // RE:A:I - I do not report anything else
	    /*
	    RE:A:I: Read extension.
	        RE: Tag identifier, stands for "Read extension."
	        A: Tag data type, denotes a single character (ASCII).
	        I: Value associated with the tag, indicating the type of read extension. In this case, I represents an insert read extension.
	    Possible values:
	    1. Insert Read Extension (RE:A:I): Indicates that the read is an insert read, meaning it is derived from the original DNA fragment without any modifications or rearrangements.
	    2. Long Fragment Read Extension (RE:A:L): Indicates that the read is a long fragment read, suggesting that it spans a longer portion of the original DNA fragment compared to other reads.
	    3. Linked Read Extension (RE:A:Q): Indicates that the read is a linked read, meaning it is part of a set of reads that are linked together, typically by sharing a common barcode or unique molecular identifier (UMI).
	    4. Chimeric Read Extension (RE:A:C): Indicates that the read is a chimeric read, meaning it contains sequences from two or more distinct DNA molecules or genomic regions.
	    */
	    //record += &format!("xf:i:{}\t",gene_id[0].edit_dist()); // li:i:0 // not linked

	    record += "li:i:0\t"; // li:i:0 // not linked
	    /*li:i:0: Linked-read identifier.
	        li: Tag identifier, stands for "Linked-read identifier."
	        i: Tag data type, denotes a 32-bit signed integer.
	        0: Value associated with the tag, representing the linked-read identifier for the read. Linked-reads are reads that are part of the same molecule or DNA fragment.
		*/
	    record += &format!("BC:Z:{}\t", umi.as_dna_string() ); // BC:Z:GCCATTCC
	    /*BC:Z:GCCATTCC: Barcode sequence.
	        BC: Tag identifier, stands for "Barcode sequence."
	        Z: Tag data type, denotes a string.
	        GCCATTCC: Value associated with the tag, representing the barcode sequence assigned to the read.
		*/
	    record += &format!("QT:Z:{}\t", &String::from_utf8_lossy(umi.qual()) ); // QT:Z:AAAAAEEE
	    /*QT:Z:AAAAAEEE: Quality tags.
	        QT: Tag identifier, stands for "Quality tags."
	        Z: Tag data type, denotes a string.
	        AAAAAEEE: Value associated with the tag, representing the quality scores of the read bases.
		*/
		record += &format!("CR:Z:{}\t",cell_id.as_dna_string() ); // CR:Z:CTCCTTTCATACTGTG
	    /*CR:Z:CTCCTTTCATACTGTG: Chromium cellular barcode sequence.
	        CR: Tag identifier, stands for "Chromium cellular barcode sequence."
	        Z: Tag data type, denotes a string.
	        CTCCTTTCATACTGTG: Value associated with the tag, representing the cellular barcode sequence assigned by the Chromium platform.
		*/
	    record += &format!("CY:Z:{}\t",  String::from_utf8_lossy(cell_id.qual()) ); // CY:Z:AAAAAEEEEEEEEEEE
	    /*CY:Z:AAAAAEEEEEEEEEEE: Chromium cellular quality scores.
	        CY: Tag identifier, stands for "Chromium cellular quality scores."
	        Z: Tag data type, denotes a string.
	        AAAAAEEEEEEEEEEE: Value associated with the tag, representing the quality scores of the cellular barcode bases.
		*/
	    record += &format!("CB:Z:{}-1\t", cell_id.as_dna_string() ); // CB:Z:CTCCTTTCATACTGTG-1
	    /*CB:Z:CTCCTTTCATACTGTG-1: Chromium gem group barcode sequence.
	        CB: Tag identifier, stands for "Chromium gem group barcode sequence."
	        Z: Tag data type, denotes a string.
	        CTCCTTTCATACTGTG-1: Value associated with the tag, representing the gem group barcode sequence assigned by the Chromium platform.
		*/
	    record += &format!("UR:Z:{}\t", umi.as_dna_string() ); // UR:Z:TTTAGGTCTTGG
	    /*UR:Z:TTTAGGTCTTGG: Chromium UMI sequence.
	        UR: Tag identifier, stands for "Chromium UMI sequence."
	        Z: Tag data type, denotes a string.
	        TTTAGGTCTTGG: Value associated with the tag, representing the unique molecular identifier (UMI) sequence assigned by the Chromium platform.
		*/
	    record += &format!("UZ:Z:{}\t", String::from_utf8_lossy(umi.qual()) ); // UY:Z:EEEEEEEEEEEE
	    /*UY:Z:EEEEEEEEEEEE: Chromium UMI quality scores.
	        UY: Tag identifier, stands for "Chromium UMI quality scores."
	        Z: Tag data type, denotes a string.
	        EEEEEEEEEEEE: Value associated with the tag, representing the quality scores of the UMI bases.
		*/
	    record += &format!("UB:Z:{}\t", umi.as_dna_string()); // UB:Z:TTTAGGTCTTGG
	    /*UB:Z:TTTAGGTCTTGG: Chromium corrected UMI sequence.
	        UB: Tag identifier, stands for "Chromium corrected UMI sequence."
	        Z: Tag data type, denotes a string.

		*/	
	    record += &format!("RG:Z:{}", "Sample4:0:1:HN2CKBGX9:1"); // RG:Z:Sample4:0:1:HN2CKBGX9:1
	    #[cfg(debug_assertions)]
	    println!("the bam line is this:\n{record}");
	    Some(record)
	}

}

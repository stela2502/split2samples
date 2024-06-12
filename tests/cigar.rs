// /tests/cigar.rs


#[cfg(test)]
mod tests {
	use rustody::genes_mapper::Cigar;
	use rustody::genes_mapper::cigar::CigarEnum;
	use rustody::genes_mapper::CigarEndFix;



	#[test]
	fn test_cvalculate_covered_1(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("5D30M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (30, 35), "corect sizes" );

		obj.restart_from_cigar("5I30M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (35, 30), "corect sizes" );

		obj.restart_from_cigar("5I15M3D15M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (35, 33), "corect sizes" );

		obj.restart_from_cigar("5D15M3I15M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (33, 35), "corect sizes" );
	}

	#[test]
	fn test_to_sam_string(){
		let mut obj = Cigar::default();

		obj.restart_from_cigar("1D15M1I30M");
		let fixed = obj.to_sam_string();

		assert_eq!( fixed, "15M1I30M", "internal insert overhanging D fixed at start");
		
		obj.restart_from_cigar("15M1I30M1D");
		let fixed = obj.to_sam_string();
		assert_eq!( fixed, "15M1I30M", "internal insert overhanging D fixed at end");
	}

	#[test]
	fn test_soft_start(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M");
		obj.soft_clip_start_end( );
		assert_eq!( obj.cigar, "24X65M");
		assert_eq!( obj.fixed, Some(CigarEndFix::Start), "start fixed");
	}

	#[test]
	fn test_soft_end(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("59M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		obj.soft_clip_start_end( );
		assert_eq!( obj.cigar, "59M30X");
		assert_eq!( obj.fixed, Some(CigarEndFix::End), "end fixed");
	}

	#[test]
	fn test_soft_both(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		obj.soft_clip_start_end( );
		assert_eq!( obj.cigar, "24X65M30X");
		assert_eq!( obj.fixed, Some(CigarEndFix::Both), "both fixed");
	}

	#[test]
	fn test_quality(){
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		assert_eq!( obj.mapping_quality(), 26 );
	}

	#[test]
	fn test_calculate_covered_nucleotides() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("14M1X39M1X7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (75, 75), "corect sizes" )
	}

	#[test]
	fn test_calculate_covered_nucleotides_deletions() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("14M1X39M1D7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (74, 75), "corect sizes" )
	}
	#[test]
	fn test_calculate_covered_nucleotides_insertions() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("14M1X39M1I7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (75, 74), "corect sizes" )
	}
	#[test]
	fn test_calculate_covered_nucleotides_real() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X8M1I39M2X16M7X");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (74, 73), "corect sizes" )
	}

	#[test]
	fn test_calculate_len_real() {
		let mut obj = Cigar::default();
		obj.restart_from_cigar("1X8M1I39M2X16M7X");
		assert_eq!( obj.len(), 74, "1X8M1I39M2X16M7X - corect sizes" )
	}


	#[test]
	fn test_default_is_worst() {
		let mut obj1 = Cigar::default();
		let mut obj2 = Cigar::default();

		obj1.restart_from_cigar( "32M" );

		assert!( obj1.better_as(&obj2), "{obj1:?}\nis better than \n{obj2:?} ({})", obj1.better_as(&obj2) );
	}

	#[test]
	fn test_compare() {
		let mut obj1 = Cigar::new("32M");
		let mut obj2 = Cigar::new("36M");
		obj1.convert_to_cigar(&vec![CigarEnum::Match; 32]);
		obj2.convert_to_cigar(&vec![CigarEnum::Match; 36]);

		assert!( obj2.better_as(&obj1), "{obj2} is better than {obj1}" );

		obj1.convert_to_cigar( &[CigarEnum::Insertion, CigarEnum::Match, CigarEnum::Match, 
			CigarEnum::Match, CigarEnum::Match, CigarEnum::Match, CigarEnum::Match, 
			CigarEnum::Deletion, CigarEnum::Match, CigarEnum::Match, CigarEnum::Match, 
			CigarEnum::Match, ] );
		obj2.convert_to_cigar( &[CigarEnum::Match, CigarEnum::Match, CigarEnum::Deletion, 
			CigarEnum::Match,CigarEnum::Insertion, CigarEnum::Match, CigarEnum::Match ]);
		assert!( obj1.better_as(&obj2), "{obj1:?} is better than {obj2:?} ({})", obj2.better_as(&obj1) );

		obj1.restart_from_cigar( "21M1D49M1I3M" ); // len 74
		obj2.restart_from_cigar( "1I20M1D52M" ); // len 73

		assert!( obj1.better_as(&obj2), "{obj1:?}\nis better than \n{obj2:?} ({})", obj1.better_as(&obj2) );
	}

	



}
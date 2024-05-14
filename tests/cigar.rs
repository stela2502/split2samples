// /tests/cigar.rs


#[cfg(test)]
mod tests {
	use rustody::genes_mapper::Cigar;

	#[test]
	fn test_soft_start(){
		let mut obj = Cigar::new("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M");
		obj.soft_clip_start_end( );
		assert_eq!( obj.cigar, "24S65M");
	}

	#[test]
	fn test_soft_end(){
		let mut obj = Cigar::new("59M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		obj.soft_clip_start_end( );
		assert_eq!( obj.cigar, "59M30S");
	}

	#[test]
	fn test_soft_both(){
		let mut obj = Cigar::new("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		obj.soft_clip_start_end( );
		assert_eq!( obj.cigar, "24S65M30S");
	}

	#[test]
	fn test_quality(){
		let obj = Cigar::new("1X1M2X1I2X1M1X1M1X1I2M1X1M1X1M1X2M1X2D2X65M1X1M2X1M1X1M1X1M1X1I2X1I2X1I1M1D3X1M1X1D1X2M2D2X1M1I");
		assert_eq!( obj.mapping_quality(), 26 )
	}

	#[test]
	fn test_calculate_covered_nucleotides() {
		let obj = Cigar::new("14M1X39M1X7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (75, 75), "corect sizes" )
	}

	#[test]
	fn test_calculate_covered_nucleotides_deletions() {
		let obj = Cigar::new("14M1X39M1D7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (74, 75), "corect sizes" )
	}
	#[test]
	fn test_calculate_covered_nucleotides_insertions() {
		let obj = Cigar::new("14M1X39M1I7M1X12M");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (75, 74), "corect sizes" )
	}
	#[test]
	fn test_calculate_covered_nucleotides_real() {
		let obj = Cigar::new("1X8M1I39M2X16M7S");
		assert_eq!( obj.calculate_covered_nucleotides( &obj.to_string() ), (74, 66), "corect sizes" )
	}

}
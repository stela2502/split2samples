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
		assert_eq!( obj.mapping_quality(), 13 )
	}
}
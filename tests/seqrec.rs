//test/seqrec.rs

#[cfg(test)]
mod tests {

    use rustody::genes_mapper::SeqRec;
	use rustody::int_to_str::IntToStr;
    #[test]
    fn test_simple(){
    	let obj = SeqRec::new(
    		b"read1",
    		b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG",
    		b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
    	);

    	let tool = IntToStr::new(b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG".to_vec(), 32);
    	assert_eq!( obj.to_u32(), tool.into_u32(),"I got what I expected {:b}, {:b}", obj.to_u32() ,tool.into_u32() );

    	assert_eq!( obj.to_u16(), tool.into_u16(),"I got what I expected {:b}, {:b}", obj.to_u16() ,tool.into_u16() );

    	assert_eq!( obj.to_u64(), tool.into_u64(),"I got what I expected {:b}, {:b}", obj.to_u64() ,tool.into_u64() );
    }

}
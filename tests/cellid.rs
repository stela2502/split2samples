// tests/cellid.rs

#[cfg(test)]
mod tests {
    use rustody::cellids10x::CellId10x;
    use rustody::traits::BinaryMatcher;

    #[test]
    fn test_get_nucleotide_2bit(){
    	// 0123456789012345
    	// CTGGGCTGTCACTAGT
    	let obj = CellId10x( 0b011100011010001111011011010101101 );
    	for (id, ch) in "CTGGGCTGTCACTAGT".chars().enumerate() {
    		if let nuc = match obj.get_nucleotide_2bit(id) {
	            Some(0) => 'A', //0b00
	            Some(1) => 'C', // 0b01
	            Some(2) => 'G', // 0b10
	            Some(3) => 'T', // 0b11
	            Some(_) => 'N',
	            None => panic!("position not reachable {id}"),
	       	}{
	       		assert_eq!( nuc, ch, "get_nucleotide_2bit({id}) should return {ch}" );
	       	}else {
	       		panic!("position not reachable {id}");
	       	}
    	}
    }
}
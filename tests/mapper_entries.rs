#[cfg(test)]
mod tests {

    use rustody::fast_mapper::mapper_entries::MapperEntry;

    #[test]
    fn check_geneids() {
        let mut mapper = MapperEntry::new( 4 );

        //let tool = IntToStr::new(b"AGCTGTGAGACTCTTCACACTATCATCATTATTCGGAGG".to_vec(), 16);
        mapper.add(12, (4,0), vec![16]);
        mapper.add(45, (3,0), vec![16] );

        assert_eq!( mapper.get(&12).unwrap().get(), vec![(4,0)] );
        assert_eq!( mapper.get(&45).unwrap().get(), vec![(3,0)] );

		mapper.add(12, (14,0), vec![16] );
		assert_eq!( mapper.get(&12).unwrap().get(), vec![(4,0), (14,0)] );


        assert_eq!( mapper.get(&14), None );

        assert_eq!( mapper.info(), [2,1,1] );
    }

    #[test]
    fn hamming_distance() {
        //let mapper = MapperEntry::new(4);

        let a:u64 = 0b11010101;
        let b:u64 = 0b10110001;

        assert_eq!( MapperEntry::hamming_distance(a, b), 3 as u32 );
    }


}

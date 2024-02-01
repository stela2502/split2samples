#[cfg(test)]
mod tests {

    use rustody::fast_mapper::mapper_entries::MapperEntry;
    use rustody::fast_mapper::mapper_entries::second_seq::SecondSeq;

    #[test]
    fn check_geneids() {
        let mut mapper = MapperEntry::new( );

        //let tool = IntToStr::new(b"AGCTGTGAGACTCTTCACACTATCATCATTATTCGGAGG".to_vec(), 16);
        mapper.add( SecondSeq(12_u64,20_u8), (4,0), vec![16]);
        mapper.add( SecondSeq(45_u64,20_u8), (3,0), vec![16] );

        assert_eq!( mapper.get(&SecondSeq(12_u64,20_u8)).unwrap()[0].get(), vec![(4,0)] );
        assert_eq!( mapper.get(&SecondSeq(45_u64,20_u8)).unwrap()[0].get(), vec![(3,0)] );

		mapper.add(SecondSeq(12_u64,20_u8), (14,0), vec![16] );
        assert_eq!( mapper.get(&SecondSeq(12_u64,20_u8)).unwrap()[0].get(), vec![ (4,0), (14,0)] );
		//assert_eq!( mapper.get(&SecondSeq(12_u64,20_u8)).unwrap()[0].get(), vec![ (14,0), (4,0)] );


        assert_eq!( mapper.get(&SecondSeq(14_u64,20_u8)), None );

        assert_eq!( mapper.info(), [2,1,1] );
    }

}


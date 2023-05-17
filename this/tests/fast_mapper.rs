#[cfg(test)]
mod tests {

    use this::fast_mapper::FastMapper;
    //use std::path::Path;
    use kmers::naive_impl::Kmer;
    #[test]
    fn check_geneids() {
        let mut mapper = FastMapper::new( 16 );

        let mut geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        //ATCCCATCCTTCATTGTTCGCCTGGA #0
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #1
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #2
        //...........................................................................

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene2".to_string() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 10 );

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC" ), Some(0) );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC" ), Some(1));

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene3".to_string() );

        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC" ), None );

        let mut gnames = Vec::<String>::with_capacity(3);
        gnames.push( "Gene1".to_string() );
        gnames.push( "Gene2".to_string() );
        gnames.push( "Gene3".to_string() );
        assert_eq!( mapper.names_store, gnames );
        mapper.print();

    }

    #[test]
    fn check_write_index() {
        let mut mapper = FastMapper::new( 16 );

        let mut geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene2".to_string() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 10);

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC" ), Some(0) );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC" ), Some(1) );

        let opath = "testData/output_index_test";
        mapper.write_index( opath.to_string() ).unwrap();
        mapper.print();

        let mut mapper2 = FastMapper::new( 16 );
        mapper2.load_index( opath.to_string() ).unwrap();

        assert_eq!(  mapper.with_data, mapper2.with_data );

        assert_eq!(  mapper.kmer_len, mapper2.kmer_len );

        assert_eq!(  mapper.kmer_len, mapper2.kmer_len );

        let key = b"ATCCCATC";
        let kmer = needletail::kmer::Kmers::new(key, 8 as u8).next();
        let _idx = Kmer::from(kmer.unwrap()).into_u64() as usize;

        assert_eq!(  mapper.names_store, mapper2.names_store );

    }

}

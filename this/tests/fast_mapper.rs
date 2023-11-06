#[cfg(test)]
mod tests {

    use this::fast_mapper::FastMapper;
    //use std::path::Path;
    use kmers::naive_impl::Kmer;
    use this::mapping_info::MappingInfo;
    use this::ofiles::Ofiles;
    use std::fs::File;
    use std::path::PathBuf;

    #[test]
    fn check_geneids() {
        let mut mapper = FastMapper::new( 16, 10 );

        let mut geneid = 0;
        let outpath = String::from("testData/");
        let log_file_str = PathBuf::from(&outpath).join(
            "Mapping_log.txt"
        );
        let log_file = match File::create( log_file_str ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {err:#?}" );
            }
        };
        let ofile = Ofiles::new( 1, "Umapped_with_cellID", "R2.fastq.gz", "R1.fastq.gz",  outpath.as_str() );
        let mut results = MappingInfo::new( log_file, 0.20, 1232423523523, ofile );

        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        //ATCCCATCCTTCATTGTTCGCCTGGA #0
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #1
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #2
        //...........................................................................

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene2".to_string() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 31 );

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC", &mut results ), Some(0) );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut results ), Some(1));

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene3".to_string() );

        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut results ), None );

        let mut gnames = Vec::<String>::with_capacity(3);
        gnames.push( "Gene1".to_string() );
        gnames.push( "Gene2".to_string() );
        gnames.push( "Gene3".to_string() );
        assert_eq!( mapper.names_store, gnames );
        mapper.print();

    }
    #[test]
    fn check_changed_start() {
        let mut mapper = FastMapper::new( 16, 10 );
        mapper.last_count = 10;
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string() );
        assert_eq!( mapper.last_count, 11 );
        assert_eq!(mapper.names_store[11] == "Gene1".to_string());
    }

    #[test]
    fn check_write_index() {
        let mut mapper = FastMapper::new( 16, 10 );

        //log_writer:File, min_quality:f32, max_reads:usize, ofile:Ofiles,
        let outpath = String::from("testData/");
        let log_file_str = PathBuf::from(&outpath).join(
            "Mapping_log.txt"
        );
        let log_file = match File::create( log_file_str ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {err:#?}" );
            }
        };
        let ofile = Ofiles::new( 1, "Umapped_with_cellID", "R2.fastq.gz", "R1.fastq.gz",  outpath.as_str() );
        let mut results = MappingInfo::new( log_file, 0.20 as f32, 1232423523523, ofile );

        let mut geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene2".to_string() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 31);

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC", &mut results ), Some(0) );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut results ), Some(1) );

        let opath = "testData/output_index_test";
        mapper.write_index( opath.to_string() ).unwrap();
        mapper.print();

        let mut mapper2 = FastMapper::new( 16, 10 );
        mapper2.load_index( opath.to_string() ).unwrap();

        assert_eq!(  mapper.with_data, mapper2.with_data );

        assert_eq!(  mapper.kmer_len, mapper2.kmer_len );

        assert_eq!(  mapper.kmer_len, mapper2.kmer_len );

        let key = b"ATCCCATC";
        let kmer = needletail::kmer::Kmers::new(key, 8 as u8).next();
        let _idx = Kmer::from(kmer.unwrap()).into_u64() as usize;

        assert_eq!(  mapper.names_store, mapper2.names_store );

    }

    #[test]
    fn check_merge() {

        let mut mapper = FastMapper::new( 16, 10 );


        let mut geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        let mut other = FastMapper::new( 16, 10 );

        other.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene2".to_string() );

        assert!( mapper.names_store.len() = 1, "Not exactly one gene in mapper 1" );
        assert!( other.names_store.len() = 1, "Not exactly one gene in mapper 2" );
        mapper.merge( other );

        assert!( mapper.names_store.len() = 2, "Not exactly two genes in the merged mapper" );

    }

}

#[cfg(test)]
mod tests {

    use rustody::fast_mapper::FastMapper;
    //use std::path::Path;
    use kmers::naive_impl::Kmer;
    //use rustody::mapping_info::MappingInfo;
    use rustody::int_to_str::IntToStr;
    static EMPTY_VEC: Vec<String> = Vec::new();


    #[test]
    fn check_geneids() {
        let mut mapper = FastMapper::new( 16, 10 );

        let mut geneid = 0;
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);



        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string(), EMPTY_VEC.clone() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        //ATCCCATCCTTCATTGTTCGCCTGGA #0
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #1
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #2
        //...........................................................................

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), 
            "Gene2".to_string(),EMPTY_VEC.clone() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        //assert_eq!( mapper.with_data, 51 );

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC", &mut tool ), Some(vec![0]) );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut tool ), Some(vec![1]) );

        //adding the same sequence as Gene2 with the name Gene3 should not make the next search return None
        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene3".to_string(),EMPTY_VEC.clone() );

        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut tool ), None );

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
        mapper.change_start_id( 10 );
        assert_eq!( mapper.last_count, 10);
        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        assert_eq!( mapper.last_count, 11 );
        assert_eq!(mapper.names_store[10] , "Gene1".to_string(), "gene was 'Gene1'");
    }

    #[test]
    fn check_write_index() {
        let mut mapper = FastMapper::new( 16, 10 );

        //log_writer:File, min_quality:f32, max_reads:usize, ofile:Ofiles

        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);
        let mut geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(),
         "Gene2".to_string(),EMPTY_VEC.clone() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        assert_eq!( mapper.with_data, 54);

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC", &mut tool ), Some(vec![0]) );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut tool ), Some(vec![1]) );

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


        let geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        let mut other = FastMapper::new( 16, 10 );

        other.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), 
            "Gene2".to_string(),EMPTY_VEC.clone() );

        assert!( mapper.names_store.len() == 1, "Not exactly one gene in mapper 1" );
        assert!( other.names_store.len() == 1, "Not exactly one gene in mapper 2" );
        mapper.merge( other );

        assert!( mapper.names_store.len() == 2, "Not exactly two genes in the merged mapper" );

    }


    #[test]
    fn check_samples() {
        let mut mapper = FastMapper::new( 16, 10 );
        let sample2 = b"GTTGTCAAGATGCTACCGTTCAGAGGGCAAGGTGTCACATTGGGCTACCGCGGGAAGTCGACCAGATCCTA";
        //let sample  = b"GTTGTCAAGATGCTACCGTTCTGAGGGCAAGGTGTCACTTTGGGCTACCGCGGGAAGTCGACCAGATCCTA";
        let sample_real = b"GTTGTCAAGATGCTACCGTTCAGAGAAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAGAAAA";
        let sequences = [
        b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 
        b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
        b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
        b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
        b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
        b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC"];

        let mut id = 1;
        for seq in sequences{
            //seq.reverse();
            mapper.add( &seq.to_vec(), format!("Sample{id}"),EMPTY_VEC.clone() );
            id +=1;
        }
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

        assert_eq!( mapper.get_strict( sequences[0], &mut tool ), Some(vec![0]) );
        assert_eq!( mapper.get_strict( sequences[1], &mut tool ), Some(vec![1]) );
        assert_eq!( mapper.get_strict( sample2, &mut tool ), None );
        assert_eq!( mapper.get_strict( &sequences[0][7..], &mut tool ), Some(vec![0]) );
        assert_eq!( mapper.get_strict( &sequences[11][7..], &mut tool ), Some(vec![11]) );
        assert_eq!( mapper.get_strict( sample_real, &mut tool ), Some(vec![0]) );
    }



    #[test]
    fn check_samples_shifted() {
        let mut mapper = FastMapper::new( 16, 10 );
        mapper.change_start_id( 10 );

        let sample2 = b"GTTGTCAAGATGCTACCGTTCAGAGGGCAAGGTGTCACATTGGGCTACCGCGGGAAGTCGACCAGATCCTA";
        //let sample  = b"GTTGTCAAGATGCTACCGTTCTGAGGGCAAGGTGTCACTTTGGGCTACCGCGGGAAGTCGACCAGATCCTA";
        let sample_real = b"GTTGTCAAGATGCTACCGTTCAGAGAAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAGAAAA";
        let sequences = [
        b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 
        b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
        b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
        b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
        b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
        b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC"];

        let mut id = 1;
        for seq in sequences{
            //seq.reverse();
            mapper.add( &seq.to_vec(), format!("Sample{id}"),EMPTY_VEC.clone() );
            id +=1;
        }
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

        assert_eq!( mapper.get_strict( sequences[0], &mut tool ), Some(vec![10]) );
        assert_eq!( mapper.get_strict( sequences[1], &mut tool ), Some(vec![11]) );
        assert_eq!( mapper.get_strict( sample2, &mut tool ), None );
        assert_eq!( mapper.get_strict( &sequences[0][7..], &mut tool ), Some(vec![10]) );
        assert_eq!( mapper.get_strict( &sequences[11][7..], &mut tool ), Some(vec![21]) );
        assert_eq!( mapper.get_strict( sample_real, &mut tool ), Some(vec![10]) );
    }

}

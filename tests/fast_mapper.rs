#[cfg(test)]
mod tests {

    use rustody::fast_mapper::FastMapper;
    //use std::path::Path;
    use kmers::naive_impl::Kmer;
    //use rustody::mapping_info::MappingInfo;
    use rustody::int_to_str::IntToStr;
    use rustody::errors::MappingError;
    static EMPTY_VEC: Vec<String> = Vec::new();


    #[test]
    fn check_samples() {
        let mut mapper = FastMapper::new( 32, 10, 0 );
        let sample2 =  b"GTTGTCAAGATGCTACCGTTCAGAGGGCAAGGTGTCACATTGGGCTACCGCGGGAAGTCGACCAGATCCTA";
        //let sample_0 =                           b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG";
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
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 27);

        let expected_results = (0..12).map(|i| vec![i]).collect::<Vec<_>>();

        for (sequence, expected_result) in sequences.iter().zip(expected_results.iter()) {
            match mapper.get(*sequence, &mut tool) {
                Ok(result) => {
                    assert_eq!(result, *expected_result, "Unexpected result for sequence: {:?}", String::from_utf8_lossy(*sequence));
                }
                Err(err) => {
                    panic!("Error occurred for sequence {:?}: {:?}", String::from_utf8_lossy(*sequence), err);
                }
            }
            
        }

        match mapper.get_strict( sample2, &mut tool ){
            Err(MappingError::NoMatch) => {
                assert_eq!( "NoMatch Error", "NoMatch Error", "expected error detected" )
            },
            Err(e) => {
                panic!("The expected error MappingError::NoMatch was not thrown! {e:?}")
            },
            Ok(v) => {
                panic!("This search should have given a MappingError::NoMatch and I got {v:?}");
            }
        }

        assert_eq!( mapper.get( &sequences[0][7..], &mut tool ).unwrap(), vec![0] );
        assert_eq!( mapper.get( &sequences[11][7..], &mut tool ).unwrap(), vec![11] );
        assert_eq!( mapper.get( sample_real, &mut tool ).unwrap(), vec![0] );
    }


    #[test]
    fn check_integrate_repeats1() {
        let mut mapper = FastMapper::new( 32, 10, 0 );
        let classes1 = vec![ "Ulk2".to_string(), "famA".to_string(), "clusterA".to_string(), "regionA".to_string()];

        let seq = &b"CCAAGAATGGTTCCTGTGTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTGAATAAACCATTCACCGGTTTTAATCCTTTTACTTCAAAACTTACACATACTGACCTAC";
        mapper.add( &seq.to_vec(), "Ulk2".to_string(), classes1 );
        let classes2 = vec![ "Ulk1".to_string(), "famA".to_string(), "clusterA".to_string(), "regionA".to_string()];
        mapper.add( &seq.to_vec(), "Ulk1".to_string(), classes2 );
        
        let mut tool = IntToStr::new(seq.to_vec(), 32);
        if let Some((first, second)) = tool.next(){
            let mapper_obj = &mapper.mapper[ first as usize ];
            let matched = mapper_obj.get( &second );
            assert!(matched.is_some(), "found the object" );
            if let Some(name_entry) = matched{
                assert_eq!( name_entry[0].key, second, "not the right SecondSeq? {:?} != {:?}", name_entry[0].key, second );
                assert_eq!( name_entry[0].get().len(), 2, "I have two gene names in one name_entry" );
            }
        }else {
            panic!("could not get a single next from IntToSeq");
        }

        mapper.make_index_te_ready();

        if let Some((first, second)) = tool.next(){
            let mapper_obj = &mapper.mapper[ first as usize ];
            let matched = mapper_obj.get( &second );
            assert!( matched.is_some(), "found a match using second next element" );
            if let Some(name_entry) = matched{
                assert_eq!( name_entry[0].key, second, "not the right SecondSeq? {:?} != {:?}", name_entry[0].key, second );
                assert_eq!( name_entry[0].get().len(), 1, "I have two gene names in one name_entry" );
                assert_eq!( name_entry[0].get(), vec![(1,2)], "The famA is the tag for the collapsed name entry" );
            }else {
                panic!("No match using the second next sequenes {first:?} and {second:?}")
            }
        }else {
            panic!("IntToStr did not give me a second (first, second) touple!");
        }
    }

    #[test]
    fn check_integrate_repeats2() {
        let mut mapper = FastMapper::new( 32, 10, 0 );
        let classes1 = vec![ "Ulk2".to_string(), "famA".to_string(), "clusterA".to_string(), "regionA".to_string()];

        let seq = &b"CCAAGAATGGTTCCTGTGTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTGAATAAACCATTCACCGGTTTTAATCCTTTTACTTCAAAACTTACACATACTGACCTAC";
        mapper.add( &seq.to_vec(), "Ulk2".to_string(), classes1 );
        let classes2 = vec![ "Ulk1".to_string(), "famA".to_string(), "clusterA".to_string(), "regionA".to_string()];
        mapper.add( &seq.to_vec(), "Ulk1".to_string(), classes2 );

        let mut mapper2 = FastMapper::new( 32, 10, 0 );
        let classes3 = vec![ "Ulk2".to_string(), "famB".to_string(), "clusterA".to_string(), "regionA".to_string()];

        let seq = &b"CCAAGAATGGTTCCTGTGTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTGAATAAACCATTCACCGGTTTTAATCCTTTTACTTCAAAACTTACACATACTGACCTAC";
        mapper2.add( &seq.to_vec(), "Ulk2".to_string(), classes3 );
        let classes4 = vec![ "Ulk1".to_string(), "famB".to_string(), "clusterA".to_string(), "regionA".to_string()];
        mapper2.add( &seq.to_vec(), "Ulk1".to_string(), classes4 );

        mapper.merge( mapper2 );

        let mut tool = IntToStr::new(seq.to_vec(), 32);
        if let Some((first, second)) = tool.next(){
            let mapper_entry = &mapper.mapper[ first as usize ];
            let name_entry = match mapper_entry.get( &second ) {
                Some(obj) => {
                     assert_eq!( 1,1, "found the object");
                     obj
                },
                None => {
                    panic!("mapper_entry.get( second ) failed!");
                }
            };
            assert_eq!( name_entry[0].key, second, "not the right SecondSeq? {:?} != {:?}", name_entry[0].key, second );
            assert_eq!( name_entry[0].get().len(), 4, "I have four gene names in one name_entry" );
        }else {
            panic!("IntToStr did not give me a single (first, second) touple!");
        }

        mapper.make_index_te_ready();
    }


    #[test]
    fn check_geneids() {
        let mut mapper = FastMapper::new( 32, 10, 0 );

        let mut geneid = 0;
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);



        mapper.add( &b"CCAAGAATGGTTCCTGTGTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTGAATAAACCATTCACCGGTTTTAATCCTTTTACTTCAAAACTTACACATACTGACCTAC".to_vec(), "Ulk2".to_string(), EMPTY_VEC.clone() );
        mapper.names4sparse.insert( "Ulk2".to_string(), geneid );

        //ATCCCATCCTTCATTGTTCGCCTGGA #0
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #1
        //CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC #2
        //...........................................................................

        mapper.add( &b"CGAGGCTGTGTATTACTGTGCCGTGGGGCTCCGGAGCCAGGAAAAGAAGAGGAtggagagggagtgggaaggaGAAAAGTCGTATACAGATTTGGGATCTTAGGCTCTGGAGACATTCAG".to_vec(), 
            "Vpreb1".to_string(),EMPTY_VEC.clone() );
        geneid +=1;
        mapper.names4sparse.insert( "Vpreb1".to_string(), geneid );

        //assert_eq!( mapper.with_data, 51 );

        assert_eq!(  mapper.get( b"GTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTG", &mut tool ).unwrap(), vec![0] );
        assert_eq!(  mapper.get( b"GGGCTCCGGAGCCAGGAAAAGAAGAGGAtggagagggagt", &mut tool ).unwrap(), vec![1] );

        //adding the same sequence as Gene2 with the name Gene3 should not make the next search return None
        mapper.add( &b"CCAAGAATGGTTCCTGTGTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTGAATAAACCATTCACCGGTTTTAATCCTTTTACTTCAAAACTTACACATACTGACCTAC".to_vec(), "Gene3".to_string(),EMPTY_VEC.clone() );

        match mapper.get( b"GTTGTATATTATTTGGTATCTTTTACTTACCTGCTTGAATACTTG", &mut tool ){
            Err(MappingError::MultiMatch) => {
                assert_eq!( "NoMatch Error", "NoMatch Error", "expected error detected" )
            },
            Err(e) => {
                panic!("The expected error MappingError::NoMatch was not thrown! {e:?}")
            },
            Ok(v) => {
                panic!("This search should have given a MappingError::NoMatch and I got {v:?}");
            }
        }

    }

    
    #[test]
    fn check_changed_start() {
        let mut mapper = FastMapper::new( 16, 10, 0 );
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);
        mapper.change_start_id( 10 );
        assert_eq!( mapper.last_count, 0);
        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        assert_eq!( mapper.last_count, 1 );
        assert_eq!(mapper.names_store[0] , "Gene1".to_string(), "gene was 'Gene1'");
        let res = mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut tool );
        assert!( res.is_ok() , "Got a result from the search" );
        if let Ok(gene_id) = res {
            assert_eq!( gene_id, vec![10] , "obtained th shifted gene id");
        }
    }
    
    #[test]
    fn check_offset() {
        let mut mapper = FastMapper::new( 16, 10, 10 );
        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);
        assert_eq!( mapper.last_count, 0);
        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        assert_eq!( mapper.last_count, 1 );
        assert_eq!(mapper.names_store[0] , "Gene1".to_string(), "gene was 'Gene1'");
        let res = mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut tool );
        assert!( res.is_ok() , "Got a result from the search" );
        if let Ok(gene_id) = res {
            assert_eq!( gene_id, vec![10] , "obtained th shifted gene id");
        }
    }

    
    #[test]
    fn check_write_index() {
        let mut mapper = FastMapper::new( 16, 10, 0 );

        //log_writer:File, min_quality:f32, max_reads:usize, ofile:Ofiles

        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);
        let mut geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        mapper.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(),
         "Gene2".to_string(),EMPTY_VEC.clone() );
        geneid +=1;
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        //assert_eq!( mapper.with_data, 26);

        assert_eq!(  mapper.get( b"ATCCCATCCTTCATTGTTCGCCTGGACTCTCAGAAGCACATCGACTTCTCCCTCCGTTCTCCTTATGGCGGCGGC", &mut tool ).unwrap(), vec![0] );
        assert_eq!(  mapper.get( b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC", &mut tool ).unwrap(), vec![1] );

        let opath = "testData/output_index_test";
        mapper.write_index( opath.to_string() ).unwrap();
        mapper.print();

        let mut mapper2 = FastMapper::new( 16, 10, 0 );
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

        let mut mapper = FastMapper::new( 16, 10, 0 );


        let geneid = 0;
        
        mapper.add( &b"ATCCCATCCTTCATTGTTCGCCTGGA".to_vec(), "Gene1".to_string(),EMPTY_VEC.clone() );
        mapper.names4sparse.insert( "Gene1".to_string(), geneid );

        let mut other = FastMapper::new( 16, 10, 0 );

        other.add( &b"CGATTACTTCTGTTCCATCGCCCACACCTTTGAACCCTAGGGCTGGGTTGAACATCTTCTGTCTCCTAGGTCTGC".to_vec(), 
            "Gene2".to_string(),EMPTY_VEC.clone() );

        assert!( mapper.names_store.len() == 1, "Not exactly one gene in mapper 1" );
        assert!( other.names_store.len() == 1, "Not exactly one gene in mapper 2" );
        mapper.merge( other );

        assert!( mapper.names_store.len() == 2, "Not exactly two genes in the merged mapper" );

    }

    #[test]
    fn check_duplicate_seqs() {
        let mut index = FastMapper::new( 32, 10, 0 );

        //             "AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC"
        let seq = b"GTTGTCAAGATGCTACCGTTCAGAGGGCAAGGTGTCACAT";
        index.add( &seq.to_vec(), "Transcript0".to_string(),vec!("Transcript0".to_string(), "Gene0".to_string(), "Family0".to_string(),  "Class0".to_string(), ) );
        index.add( &seq.to_vec(), "Transcript1".to_string(),vec!( "Transcript1".to_string(), "Gene1".to_string(), "Family0".to_string(),  "Class0".to_string(), ) );

        let mut tool = IntToStr::new( b"AAGGCCTT".to_vec(), 32);

        match index.get( seq,  &mut tool ){
            Err(MappingError::MultiMatch) => {
                assert_eq!( "MultiMatch Error", "MultiMatch Error", "expected error detected" )
            },
            Err(e) => {
                panic!("The expected error MappingError::MultiMatch was not thrown! {e:?}")
            },
            Ok(v) => {
                panic!("This search should have given a MappingError::MultiMatch and I got {v:?}");
            }
        }

        index.make_index_te_ready();

        assert_eq!( index.get( seq,  &mut tool ).unwrap(), vec![2_usize] ); // family0


    }

}

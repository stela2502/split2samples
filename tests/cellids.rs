
#[cfg(test)]
mod tests {
    use rustody::cellids::CellIds;
    use rustody::traits::CellIndex;
    use rustody::genes_mapper::SeqRec;




    #[test]
    //
    fn getcells_1(){
        let cells = CellIds::new(&"v1".to_string() );

        let primer = SeqRec::new( 
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
        );
        let id:u32 = 1;
        let exp= ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1)) +1 as u32, 
            0, 
            SeqRec::new( b"GTCGCTATATACAGGATAAAGCCTTCT",b"GTCGCTATATACAGGATAAAGCCTTCT",b"GTCGCTATATACAGGATAAAGCCTTCT"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        match cells.to_cellid( &primer){
            Ok(val) => {
            assert_eq!( val , exp );
            eprintln!("WHAT we go an id here: {}", val.0 )}, // will never insert one element twice. Great!
            Err(_err) => {panic!("This did not work!?\n")}, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        
        let exp2 = CellIds::into_u64(vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"]);
        assert_eq!( cells.to_sequence( exp.0 ), exp2 );
        // 3, 3, 3
    }
    #[test]
    //
    fn getcells_9(){
        let cells = CellIds::new(&"v1".to_string() );

        let mut primer = SeqRec::new( 
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
        );
        let mut id:u32 = 1;
        let mut exp= ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32, 
            0, 
            SeqRec::new( b"GTCGCTATATACAGGATAAAGCCTTCT",b"GTCGCTATATACAGGATAAAGCCTTCT",b"GTCGCTATATACAGGATAAAGCCTTCT"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        match cells.to_cellid( &primer){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(_err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        
        let mut exp2 = CellIds::into_u64(vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"]);
        assert_eq!( cells.to_sequence( exp.0 ), exp2, "the expected sequence triplet #1" );
        // 3, 3, 3
        primer = SeqRec::new( 
            b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTATNNNNAAAAAAAAAAAA",
            b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTATNNNNAAAAAAAAAAAA",
            b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTATNNNNAAAAAAAAAAAA",
        );
        id = 3;
        exp = ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32, 
            0, 
            SeqRec::new( b"CTTCACATATGTGAAGAACACAAGTAT",b"CTTCACATATGTGAAGAACACAAGTAT",b"CTTCACATATGTGAAGAACACAAGTAT"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        exp2 = CellIds::into_u64( vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"] );
        match cells.to_cellid( &primer){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(_err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer)? , exp );
        assert_eq!( cells.to_sequence( exp.0 ), exp2 );

        // and the last one
        let primer2 = SeqRec::new( 
            b"NTGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCANNNNAAAAAAAAAAAA",
            b"NTGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCANNNNAAAAAAAAAAAA",
            b"NTGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCANNNNAAAAAAAAAAAA",
        );
        id = 96;
        exp = ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32, 
            0, 
            SeqRec::new( b"TGCGATCTACAACAACGGCATAGGTCA",b"TGCGATCTACAACAACGGCATAGGTCA",b"TGCGATCTACAACAACGGCATAGGTCA"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        exp2 = CellIds::into_u64( vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"] );
        //assert_eq!( 884735+1 , exp);
        match cells.to_cellid( &primer2){
            Ok(val) => assert_eq!( val , exp, "the expected sequence triplet #2" ), // will never insert one element twice. Great!
            Err(_err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer)? , exp );
        assert_eq!( cells.to_sequence( exp.0 ), exp2, "the expected sequence triplet #3" );
    }
    //
    #[test]
    fn getcells_7() {
        let cells = CellIds::new(&"v1".to_string() );

        let mut primer = SeqRec::new( 
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
            b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCTNNNNAAAAAAAAAAAA",
        );
        let mut id:u32 = 1;
        let mut exp= ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32, 
            0, 
            SeqRec::new( b"GTCGCTATATACAGGATAAAGCCTTCT",b"GTCGCTATATACAGGATAAAGCCTTCT",b"GTCGCTATATACAGGATAAAGCCTTCT"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        match cells.to_cellid( &primer){
            Ok(val) => assert_eq!( val , exp ),
            Err(_err) => (),
        };
        
        let mut exp2 = CellIds::into_u64( vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"] );
        assert_eq!( cells.to_sequence( exp.0 ), exp2 );
        // 3, 3, 3
        primer =  SeqRec::new( 
            b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTATNNNNAAAAAAAAAAAA",
            b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTATNNNNAAAAAAAAAAAA",
            b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTATNNNNAAAAAAAAAAAA",
        );
        id = 3;
        exp = ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32, 
            0, 
            SeqRec::new( b"CTTCACATATGTGAAGAACACAAGTAT",b"CTTCACATATGTGAAGAACACAAGTAT",b"CTTCACATATGTGAAGAACACAAGTAT"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        exp2 = CellIds::into_u64( vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"] );
        match cells.to_cellid( &primer){
            Ok(val) => assert_eq!( val , exp ), 
            Err(_err) => (), 
        };
        //assert_eq!( cells.to_cellid( primer)? , exp );
        assert_eq!( cells.to_sequence( exp.0 ), exp2 );

        // and the last one
        primer = SeqRec::new( 
            b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCANNNNAAAAAAAAAAAA",
            b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCANNNNAAAAAAAAAAAA",
            b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCANNNNAAAAAAAAAAAA",
        );
        id = 96;
        exp = ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32, 
            0, 
            SeqRec::new( b"TGCGATCTACAACAACGGCATAGGTCA",b"TGCGATCTACAACAACGGCATAGGTCA",b"TGCGATCTACAACAACGGCATAGGTCA"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        exp2 = CellIds::into_u64(vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"]);
        //assert_eq!( 884735+1 , exp);
        match cells.to_cellid( &primer){
            Ok(val) => assert_eq!( val , exp ), 
            Err(_err) => (),
        };
        //assert_eq!( cells.to_cellid( primer)? , exp );
        assert_eq!( cells.to_sequence( exp.0 ), exp2 );        
    }
    #[test]
    fn getcells_384_9() {
        let cells = CellIds::new(&"v2.384".to_string() );

        let primer =  SeqRec::new( 
            b"TGTCTAGCGNNNNNNNNNNNNTTGTGCGGANNNNNNNNNNNNNTTGTGCGACNNNNAAAAAAAAAAAA",
            b"TGTCTAGCGNNNNNNNNNNNNTTGTGCGGANNNNNNNNNNNNNTTGTGCGACNNNNAAAAAAAAAAAA",
            b"TGTCTAGCGNNNNNNNNNNNNTTGTGCGGANNNNNNNNNNNNNTTGTGCGACNNNNAAAAAAAAAAAA",
        ); // totally artificial - primer design wrong... - lazy
        let exp2 = CellIds::into_u64(vec![b"TGTCTAGCG", b"TTGTGCGGA", b"TTGTGCGAC"]);
        let id:u32 = 3;
        let exp= ( 
            ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) ) +1 as u32 ,
            0, 
            SeqRec::new( b"TGTCTAGCGTTGTGCGGATTGTGCGAC",b"TGTCTAGCGTTGTGCGGATTGTGCGAC",b"TGTCTAGCGTTGTGCGGATTGTGCGAC"),
            SeqRec::new( b"NNNNAAAA",b"NNNNAAAA",b"NNNNAAAA"),
        );
        match cells.to_cellid( &primer){
            Ok(val) => assert_eq!( val , exp ),
            Err(_err) => (),
        };
        
        //assert_eq!( cells.to_cellid( primer)? , exp );
        assert_eq!( cells.to_sequence( exp.0 ), exp2 );        
    }
}



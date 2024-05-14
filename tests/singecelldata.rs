#[cfg(test)]
mod tests {
	use rustody::singlecelldata::SingleCellData;
	use rustody::singlecelldata::cell_data::GeneUmiHash;
	use rustody::fast_mapper::FastMapper;
	use rustody::mapping_info::MappingInfo;
    //use rustody::singlecelldata::IndexedGenes;

	static EMPTY_VEC: Vec<String> = Vec::new();


    #[test]
    fn singlecelldata_to_sparse() {
        let mut celldata = SingleCellData::new( 1 ); // only one thread here


        let mut report = MappingInfo::new(None, 56.0, 10000, None );

        // add one cell with two genes and each 20 umi counts

        for umi in 0..20 {
        	assert!( celldata.try_insert( &1_u64, GeneUmiHash(0, umi as u64), &mut report), "I add Gene1 (0) umi {umi}" );
        }
        assert!( ! celldata.try_insert( &1_u64, GeneUmiHash(0, 0 as u64), &mut report ),"I add Gene1 (0) umi 0 - again - and failed" );
        for umi in 20..30 {
            assert!( celldata.try_insert( &1_u64, GeneUmiHash(3, umi as u64), &mut report ),"I add Gene4 (3) umi {umi} ");
        }

        // add one other cell with one gene and 20 umi counts

        for umi in 0..20 {
        	assert!( celldata.try_insert( &13452355_u64, GeneUmiHash(2, umi as u64), &mut report), "I add Gene3 (2) umi {umi}" );
        	assert!( celldata.try_insert( &13452355_u64, GeneUmiHash(0, umi as u64), &mut report), "I add Gene1 (0) umi {umi}" );
        }
        assert!( ! celldata.try_insert( &13452355_u64, GeneUmiHash(0, 0 as u64), &mut report ),"I add Gene1 umi 0 - again - and failed" );


	    let sequences = [
        b"AAGAGTCGACTGCCATGTCCCCTCCGCGGGTCCGTGCCCCCCAAG", b"ACCGATTAGGTGCGAGGCGCTATAGTCGTACGTCGTTGCCGTGCC", 
        b"AGGAGGCCCCGCGTGAGAGTGATCAATCCAGGATACATTCCCGTC", b"TTAACCGAGGCGTGAGTTTGGAGCGTACCGGCTTTGCGCAGGGCT",
        b"GGCAAGGTGTCACATTGGGCTACCGCGGGAGGTCGACCAGATCCT", b"GCGGGCACAGCGGCTAGGGTGTTCCGGGTGGACCATGGTTCAGGC",
        b"ACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC", b"TCGCTGCCGTGCTTCATTGTCGCCGTTCTAACCTCCGATGTCTCG",
        b"GCCTACCCGCTATGCTCGTCGGCTGGTTAGAGTTTACTGCACGCC", b"TCCCATTCGAATCACGAGGCCGGGTGCGTTCTCCTATGCAATCCC",
        b"GGTTGGCTCAGAGGCCCCAGGCTGCGGACGTCGTCGGACTCGCGT", b"CTGGGTGCCTGGTCGGGTTACGTCGGCCCTCGGGTCGCGAAGGTC"];
        let mut mapper = FastMapper::new( 32, 10, 0 );
        let mut id = 1;
        for seq in sequences{
            //seq.reverse();
            mapper.add( &seq.to_vec(), format!("Gene{id}"),EMPTY_VEC.clone() );
            id +=1;
        }

        //to_str<'live>(&mut self, gene_info:&GeneIds, names: &Vec<String> ) 
        let  names= vec!("Gene1".to_string(), "Gene3".to_string(), "Gene4".to_string() );
        // this string counts: genes, cell, lines
        let  exp2:String = "3 2 4".to_string();
        celldata.update_genes_to_print( &mapper.as_indexed_genes() , &names);
        let  val = celldata.mtx_counts( &mapper.as_indexed_genes(), 1, 1 );

        assert_eq!( val,  exp2 ); 

    }

//     #[test]
//     fn singlecelldata_to_sparse2cells() {
//         let mut celldata = SingleCellData::new( );

//         let mut cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 7 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         println!("Gene id == {}", geneid);
//         for umi in 0..20 {
//             println!("I add Gene1 ({}) umi {}", geneid, umi );
//             cell1.add( geneid, umi as u64);
//         }
//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }

//         cell1 = match celldata.get( 2 , "Cell2".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }

//         let names = vec!("Gene3".to_string());
//         let val = celldata.mtx_counts( &mut genes, &names, 1 as usize );
//         let exp2 = "1 2 2".to_string();
        
//         assert_eq!( val,  exp2 ); 

//     }
//     use std::fs;
//     use std::fs::File;
//     use std::path::PathBuf;
//     use flate2::read::GzDecoder;
//     use std::io::Read;

//     #[test]
//     fn singlecelldata_to_sparse2cells_outfiles() {
//         let mut celldata = SingleCellData::new( );

//         let mut cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 12 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         println!("Gene id == {}", geneid);

//         for umi in 0..20 {
//             //println!("I add Gene1 ({}) umi {}", geneid, umi );
//             assert_eq!( cell1.add( geneid, umi as u64), true);
//         }
//         assert_eq!( cell1.add( geneid, 1 as u64), false);

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..20 {
//             cell1.add( geneid, umi as u64);
//         }

//         cell1 = match celldata.get( 2 , "Cell2".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..10 {
//             cell1.add( geneid, umi as u64);
//         }

//         let names = vec!("Gene3".to_string());

//         let file_path_sp = PathBuf::from( "../testData/output_sparse");

//         match celldata.write_sparse_sub ( file_path_sp, &mut genes, &names, 1) {
//             Ok(_) => (),
//             Err(err) => panic!("Error in the data write: {}", err)
//         };

//         match fs::create_dir("../testData/output_sparse/"){
//             Ok(_) => (),
//             Err(err) => println!("{}", err),
//         };

//         let file = match File::open("../testData/output_sparse/matrix.mtx.gz"){
//             Ok(f) => f,
//             Err(_err) => panic!("expected outfile is missing"),
//         };
//         let mut gz = GzDecoder::new(file);

//         let val = "%%MatrixMarket matrix coordinate integer general\n1 2 2\n1 1 20\n1 2 10\n".to_string();

//         let mut buffer = String::new();

//         match gz.read_to_string(&mut buffer) {
//             Ok(_string) => {},
//             Err(_err) => {},
//         };

//         assert_eq!( val,  buffer ); 

//     }


// #[test]
//     fn singlecelldata_only_one_call_to_cellid_add() {
//         let mut celldata = SingleCellData::new( );

//         let mut cell1 = match celldata.get( 1 , "Cell1".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         let mut genes = GeneIds::new( 12 );

//         let mut geneid = 0;
        
//         genes.add( b"AGCTGCTAGCCGATATT", "Gene1".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         genes.add( b"CTGTGTAGATACTATAGATAA", "Gene2".to_string() );
//         genes.names4sparse.insert( "Gene1".to_string(), geneid );
//         // the next two should not be in the output
//         genes.add( b"CGCGATCGGATAGCTAGATAGG", "Gene3".to_string() );
//         genes.add( b"CATACAACTACGATCGAATCG", "Gene4".to_string() );

//         geneid = genes.get_id( "Gene1".to_string() );
//         assert_eq!( cell1.add( geneid, 1 as u64), true);

//         geneid = genes.get_id( "Gene3".to_string() );
//         assert_eq!( cell1.add( geneid, 1 as u64), true);

//         cell1 = match celldata.get( 2 , "Cell2".to_string() ){
//             Ok( cell) => cell,
//             Err(err) => panic!("{}", err ),
//         };

//         geneid = genes.get_id( "Gene3".to_string() );
//         for umi in 0..10 {
//             cell1.add( geneid, umi as u64);
//         }

//         let names = vec!("Gene1".to_string(), "Gene3".to_string());

//         let file_path_sp = PathBuf::from( "../testData/output_sparse2");

//         match celldata.write_sparse_sub ( file_path_sp, &mut genes, &names, 1) {
//             Ok(_) => (),
//             Err(err) => panic!("Error in the data write: {}", err)
//         };

//         match fs::create_dir("../testData/output_sparse2/"){
//             Ok(_) => (),
//             Err(err) => println!("{}", err),
//         };

//         let file = match File::open("../testData/output_sparse2/matrix.mtx.gz"){
//             Ok(f) => f,
//             Err(_err) => panic!("expected outfile is missing"),
//         };
//         let mut gz = GzDecoder::new(file);

//         let val = "%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 1\n2 1 1\n2 2 10\n".to_string();

//         let mut buffer = String::new();

//         match gz.read_to_string(&mut buffer) {
//             Ok(_string) => {},
//             Err(_err) => {},
//         };

//         assert_eq!( val,  buffer ); 

//     }
 }
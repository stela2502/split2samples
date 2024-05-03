// tests/genes_mapper.rs

#[cfg(test)]
mod tests {

	use rustody::genes_mapper::GenesMapper;
	use needletail::parse_fastx_file;
	use std::path::Path;

	#[test]
	fn test_create(){
		let mut genes = GenesMapper::new(0);

		let ex = "testData/genes.fasta";
		let mut i=0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

	    assert_eq!( genes.get_gene_count() , i, "Got the expecetd gene counts out" );
	    assert_eq!( genes.get_max_gene_id(), i, "Got the last id out" );
	}

	#[test]
	fn test_create_offset(){
		let mut genes = GenesMapper::new(100);

		let ex = "testData/genes.fasta";
		let mut i=0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

	    assert_eq!( genes.get_gene_count() , i, "Got the expecetd gene counts out" );

	    assert_eq!( genes.get_max_gene_id(), 100 +i, "Got the last id out" );
	}

	#[test]
	fn test_indexed_genes(){
		let mut genes = GenesMapper::new(100);

		let ex = "testData/genes.fasta";
		#[allow(unused_variables)] //otherwise the compiler falsly brags about this.
		let mut i = 0;
		if Path::new(ex).exists(){
	    	let mut expr_file = parse_fastx_file(ex).expect("valid path/file");

	    	while let Some(e_record) = expr_file.next() {
	    		i+=1;
		        let seqrec = e_record.expect("invalid record");
	        	match std::str::from_utf8(seqrec.id()){
		            Ok(st) => {
	                	if let Some(id) = st.to_string().split('|').next(){
		                    genes.add( &seqrec.seq().to_vec(), id.to_string(), id.to_string(), 0 );
	                	}
	            	},
	            	Err(err) => eprintln!("The expression entry's id could not be read: {err}"),
	        	}
	        }
	    }else {
	    	panic!("Expression file could not be read - ignoring")
	    }

	    let indexed = genes.as_indexed_genes();

	    let indexed_names = indexed.get_all_gene_names();

	    for i in 0..5{
	    	assert_eq!( indexed.ids_for_gene_names( &vec![indexed_names[i].to_string()])[0], genes.extern_id_for_gname( &indexed_names[i] ).unwrap(), "the id was as expoected" );
	    }

	}
}

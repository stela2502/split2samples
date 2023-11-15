/// here I plann to define all my traits - or the one I have planned for now :-D


pub trait Index{

	fn new( kmer_len:usize, allocate:usize )-> Self;
	fn get(&self, seq: &[u8]) ->  Option< usize >;
	fn add(&mut self, seq: &Vec<u8>, name: std::string::String, class_ids: Vec<String> ) -> usize;
	fn get_id( &self, name: String ) -> usize;
	fn print( &self );
	fn change_start_id ( &mut self, new_start :usize );
	fn names_len(&self) -> usize;
	fn names(&self) -> Vec<String>;
	fn to_header_n( &self, names: &Vec<String> ) -> std::string::String;
	fn names4sparse( &self ) -> Vec<String>;
	fn reset_names4sparse( &mut self );
	fn add_2_names4sparse( &mut self, name:&str );
	fn max_id( &self ) -> usize; // return the max:id foir the sparse export of the data
}
/// here I plann to define all my traits - or the one I have planned for now :-D


pub trait Index : Sync{

	fn get(&self, seq: &[u8]) ->  Option< usize >;
	fn add(&mut self, seq: &Vec<u8>, name: std::string::String, class_ids: Vec<String> ) -> usize;
	fn get_id( &self, name: String ) -> usize;
	fn print( &self );
	fn change_start_id ( &mut self, new_start :usize );
	fn names_len(&self) -> usize;
	fn names(&self) -> Vec<String>;
	fn to_header_n( &self, names: &Vec<String> ) -> std::string::String;
	fn max_id( &self ) -> usize; // return the max:id foir the sparse export of the data
	fn write_index( &mut self, path: String ) -> Result< (), &str>;
	fn load_index( &mut self, path: String ) -> Result< (), &str>;
}

pub trait CellIndex: Sync{
	fn to_cellid (&self, r1: &[u8]  )-> Result<( u32, u64), &str>;
}
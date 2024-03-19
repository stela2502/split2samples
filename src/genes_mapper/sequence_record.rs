// sequence_record.rs

/// A really simple copy of needletails SequenceRecord to allow a hard copy of the data
/// This would be necessary for me to used

#[derive(Debug, Clone)]
pub struct SeqRec{
	id: Vec<u8>,
	seq: Vec<u8>,
	qual:Vec<u8>,
}

impl SeqRec{
	pub fn new(id:&[u8], seq:&[u8], qual:&[u8])->Self{
		Self{
			id: id.to_vec(),
			seq: seq.to_vec(),
			qual: qual.to_vec(),
		}

	}
	pub fn id(&self) -> &[u8] {
        &self.id
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    pub fn qual(&self) -> &[u8] {
        &self.qual
    }
}
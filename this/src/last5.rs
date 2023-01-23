

pub struct Last5 {
	data : Vec::<f32>,
	pos: usize,
	max: usize
}

impl Last5 {
	pub fn new( size:usize ) -> Self{
		let max = size -1;
		let data: Vec::<f32> = vec!(0.0; size);
		let pos = 0;
        Self{
            data,
            pos,
            max
        }
    }

    pub fn push(&mut self, val: f32 ) {
    	self.data[self.pos] = val;
    	self.pos += 1;
    	if self.pos > self.max {
    		self.pos = 0;
    	}
    }

    pub fn all(&mut self, val: f32 ) -> bool {
    	for dat in &self.data {
    		if *dat != val {
    			return false
    		}
    	}
    	return true
    }
}



/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use crate::traits::CellIndex;
use crate::int_to_str::IntToStr;


/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellIds10x{
    pub ver: String,
}

impl CellIndex for CellIds10x{

    fn to_cellid (&self, r1: &[u8]  )-> Result<( u32, u64 ), &str>{

        let mut tool = IntToStr::new( r1.to_vec(), 9);
        let cell:u32 = tool.into_u32();
        tool.drop_n(4); // drop the 4 u8 that defined the cellid

        let umi:u64  = tool.into_u64();
        //println!("We found a cellid {} and an umi {}", cell, umi);
        return Ok( (cell, umi) )

    }
}


impl  CellIds10x{
    pub fn new(ver:&String )-> Self {
        //println!("Yes you initialized a cellIds10x object!");
        Self {
            ver: ver.to_string(),
        }
    }
}

// macro_rules! to_cellid {
//    ($r1: expr, $r2: expr) => {
//       //1-9,22-30,44-52
//       $r1.to_cellid($r1, $r2, [0,8], [21,29], [43,51]  )
//    };
// }

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn getsamples() {
//         let cells = crate::cell_ids::new();

//         let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
//         let mut id:u32 = 1;
//         let mut exp= ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
        
//         let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
//         assert_eq!( cells.to_sequence( exp ), exp2 );
//         // 3, 3, 3
//         primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
//         id = 3;
//         exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
//         //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
//         assert_eq!( cells.to_sequence( exp ), exp2 );

//         // and the last one
//         primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
//         id = 96;
//         exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
//         assert_eq!( 884735+1 , exp);
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
//         //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
//         assert_eq!( cells.to_sequence( exp ), exp2 );        
//     }
// }



/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use kmers::naive_impl::Kmer;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;


// and here the data
pub struct CellIds<'a>{    
    c1s: Vec<&'a [u8; 9]>,
    c2s: Vec<&'a [u8; 9]>,
    c3s: Vec<&'a [u8; 9]>,
    csl1: BTreeMap<u64, u32>,
    csl2: BTreeMap<u64, u32>,
    csl3: BTreeMap<u64, u32>,
}


// here the functions
impl CellIds<'_>{

    pub fn new()-> Self {
        let  c1s: Vec<&[u8; 9]>;
        let  c2s: Vec<&[u8; 9]>;
        let  c3s: Vec<&[u8; 9]>;

        let mut csl1 = BTreeMap::<u64, u32>::new();
        let mut csl2 = BTreeMap::<u64, u32>::new();
        let mut csl3 = BTreeMap::<u64, u32>::new();

        c1s = vec![
              b"GTCGCTATA", b"CTTGTACTA", b"CTTCACATA", b"ACACGCCGG",
              b"CGGTCCAGG", b"AATCGAATG", b"CCTAGTATA", b"ATTGGCTAA",
              b"AAGACATGC", b"AAGGCGATC", b"GTGTCCTTA", b"GGATTAGGA",
              b"ATGGATCCA", b"ACATAAGCG", b"AACTGTATT", b"ACCTTGCGG",
              b"CAGGTGTAG", b"AGGAGATTA", b"GCGATTACA", b"ACCGGATAG",
              b"CCACTTGGA", b"AGAGAAGTT", b"TAAGTTCGA", b"ACGGATATT",
              b"TGGCTCAGA", b"GAATCTGTA", b"ACCAAGGAC", b"AGTATCTGT",
              b"CACACACTA", b"ATTAAGTGC", b"AAGTAACCC", b"AAATCCTGT",
              b"CACATTGCA", b"GCACTGTCA", b"ATACTTAGG", b"GCAATCCGA",
              b"ACGCAATCA", b"GAGTATTAG", b"GACGGATTA", b"CAGCTGACA",
              b"CAACATATT", b"AACTTCTCC", b"CTATGAAAT", b"ATTATTACC",
              b"TACCGAGCA", b"TCTCTTCAA", b"TAAGCGTTA", b"GCCTTACAA",
              b"AGCACACAG", b"ACAGTTCCG", b"AGTAAAGCC", b"CAGTTTCAC",
              b"CGTTACTAA", b"TTGTTCCAA", b"AGAAGCACT", b"CAGCAAGAT",
              b"CAAACCGCC", b"CTAACTCGC", b"AATATTGGG", b"AGAACTTCC",
              b"CAAAGGCAC", b"AAGCTCAAC", b"TCCAGTCGA", b"AGCCATCAC",
              b"AACGAGAAG", b"CTACAGAAC", b"AGAGCTATG", b"GAGGATGGA",
              b"TGTACCTTA", b"ACACACAAA", b"TCAGGAGGA", b"GAGGTGCTA",
              b"ACCCTGACC", b"ACAAGGATC", b"ATCCCGGAG", b"TATGTGGCA",
              b"GCTGCCAAT", b"ATCAGAGCT", b"TCGAAGTGA", b"ATAGACGAG",
              b"AGCCCAATC", b"CAGAATCGT", b"ATCTCCACA", b"ACGAAAGGT",
              b"TAGCTTGTA", b"ACACGAGAT", b"AACCGCCTC", b"ATTTAGATG",
              b"CAAGCAAGC", b"CAAAGTGTG", b"GGCAAGCAA", b"GAGCCAATA",
              b"ATGTAATGG", b"CCTGAGCAA", b"GAGTACATT", b"TGCGATCTA"
            ];
        c2s = vec![
              b"TACAGGATA", b"CACCAGGTA", b"TGTGAAGAA", b"GATTCATCA",
              b"CACCCAAAG", b"CACAAAGGC", b"GTGTGTCGA", b"CTAGGTCCT",
              b"ACAGTGGTA", b"TCGTTAGCA", b"AGCGACACC", b"AAGCTACTT",
              b"TGTTCTCCA", b"ACGCGAAGC", b"CAGAAATCG", b"ACCAAAATG",
              b"AGTGTTGTC", b"TAGGGATAC", b"AGGGCTGGT", b"TCATCCTAA",
              b"AATCCTGAA", b"ATCCTAGGA", b"ACGACCACC", b"TTCCATTGA",
              b"TAGTCTTGA", b"ACTGTTAGA", b"ATTCATCGT", b"ACTTCGAGC",
              b"TTGCGTACA", b"CAGTGCCCG", b"GACACTTAA", b"AGGAGGCGC",
              b"GCCTGTTCA", b"GTACATCTA", b"AATCAGTTT", b"ACGATGAAT",
              b"TGACAGACA", b"ATTAGGCAT", b"GGAGTCTAA", b"TAGAACACA",
              b"AAATAAATA", b"CCGACAAGA", b"CACCTACCC", b"AAGAGTAGA",
              b"TCATTGAGA", b"GACCTTAGA", b"CAAGACCTA", b"GGAATGATA",
              b"AAACGTACC", b"ACTATCCTC", b"CCGTATCTA", b"ACACATGTC",
              b"TTGGTATGA", b"GTGCAGTAA", b"AGGATTCAA", b"AGAATGGAG",
              b"CTCTCTCAA", b"GCTAACTCA", b"ATCAACCGA", b"ATGAGTTAC",
              b"ACTTGATGA", b"ACTTTAACT", b"TTGGAGGTA", b"GCCAATGTA",
              b"ATCCAACCG", b"GATGAACTG", b"CCATGCACA", b"TAGTGACTA",
              b"AAACTGCGC", b"ATTACCAAG", b"CACTCGAGA", b"AACTCATTG",
              b"CTTGCTTCA", b"ACCTGAGTC", b"AGGTTCGCT", b"AAGGACTAT",
              b"CGTTCGGTA", b"AGATAGTTC", b"CAATTGATC", b"GCATGGCTA",
              b"ACCAGGTGT", b"AGCTGCCGT", b"TATAGCCCT", b"AGAGGACCA",
              b"ACAATATGG", b"CAGCACTTC", b"CACTTATGT", b"AGTGAAAGG",
              b"AACCCTCGG", b"AGGCAGCTA", b"AACCAAAGT", b"GAGTGCGAA",
              b"CGCTAAGCA", b"AATTATAAC", b"TACTAGTCA", b"CAACAACGG"
            ];
        c3s = vec![
              b"AAGCCTTCT", b"ATCATTCTG", b"CACAAGTAT", b"ACACCTTAG",
              b"GAACGACAA", b"AGTCTGTAC", b"AAATTACAG", b"GGCTACAGA",
              b"AATGTATCG", b"CAAGTAGAA", b"GATCTCTTA", b"AACAACGCG",
              b"GGTGAGTTA", b"CAGGGAGGG", b"TCCGTCTTA", b"TGCATAGTA",
              b"ACTTACGAT", b"TGTATGCGA", b"GCTCCTTGA", b"GGCACAACA",
              b"CTCAAGACA", b"ACGCTGTTG", b"ATATTGTAA", b"AAGTTTACG",
              b"CAGCCTGGC", b"CTATTAGCC", b"CAAACGTGG", b"AAAGTCATT",
              b"GTCTTGGCA", b"GATCAGCGA", b"ACATTCGGC", b"AGTAATTAG",
              b"TGAAGCCAA", b"TCTACGACA", b"CATAACGTT", b"ATGGGACTC",
              b"GATAGAGGA", b"CTACATGCG", b"CAACGATCT", b"GTTAGCCTA",
              b"AGTTGCATC", b"AAGGGAACT", b"ACTACATAT", b"CTAAGCTTC",
              b"ACGAACCAG", b"TACTTCGGA", b"AACATCCAT", b"AGCCTGGTT",
              b"CAAGTTTCC", b"CAGGCATTT", b"ACGTGGGAG", b"TCTCACGGA",
              b"GCAACATTA", b"ATGGTCCGT", b"CTATCATGA", b"CAATACAAG",
              b"AAAGAGGCC", b"GTAGAAGCA", b"GCTATGGAA", b"ACTCCAGGG",
              b"ACAAGTGCA", b"GATGGTCCA", b"TCCTCAATA", b"AATAAACAA",
              b"CTGTACGGA", b"CTAGATAGA", b"AGCTATGTG", b"AAATGGAGG",
              b"AGCCGCAAG", b"ACAGTAAAC", b"AACGTGTGA", b"ACTGAATTC",
              b"AAGGGTCAG", b"TGTCTATCA", b"TCAGATTCA", b"CACGATCCG",
              b"AACAGAAAC", b"CATGAATGA", b"CGTACTACG", b"TTCAGCTCA",
              b"AAGGCCGCA", b"GGTTGGACA", b"CGTCTAGGT", b"AATTCGGCG",
              b"CAACCTCCA", b"CAATAGGGT", b"ACAGGCTCC", b"ACAACTAGT",
              b"AGTTGTTCT", b"AATTACCGG", b"ACAAACTTT", b"TCTCGGTTA",
              b"ACTAGACCG", b"ACTCATACG", b"ATCGAGTCT", b"CATAGGTCA"
            ];

        let mut i: u32 = 0;
        for kmer_u8 in &c1s{
            for kmer in needletail::kmer::Kmers::new( *kmer_u8, 9 ) { // exactly 1
                let km = Kmer::from(kmer).into_u64();
                csl1.insert( km, i );
                i += 1;
            }
        }
        i = 0;
        for kmer_u8 in &c2s{
            for kmer in needletail::kmer::Kmers::new( *kmer_u8, 9 ) { // exactly 1
                let km = Kmer::from(kmer).into_u64();
                csl2.insert( km, i );
                i += 1;
            }
        }
        i = 0;
        for kmer_u8 in &c3s{
            for kmer in needletail::kmer::Kmers::new( *kmer_u8, 9 ) { // exactly 1
                let km = Kmer::from(kmer).into_u64();
                csl3.insert( km, i );
                i += 1;
            }
        }

        Self {
            c1s,
            c2s,
            c3s,
            csl1,
            csl2,
            csl3
        }
    }

    pub fn to_cellid (&self, r1: &[u8], c1: Vec<usize>, c2: Vec<usize>, c3: Vec<usize>  )-> Result< u32, &str>{
        let mut cell_id:u32 = 0;
        let max:u32 = 96;
        //println!("to_cellid should print something! {:?}", &r1[c1[0]..c1[1]]);
        
        //println!("to_cellid the c1 seq: {:?}", std::str::from_utf8( &r1[c1[0]..c1[1]] ) );
        //println!("to_cellid the c2 seq: {:?}", std::str::from_utf8( &r1[c2[0]..c2[1]] ) );
        //println!("to_cellid the c3 seq: {:?}", std::str::from_utf8( &r1[c3[0]..c3[1]] ) );

        for km in [ &r1[c1[0]..c1[1]], &r1[c2[0]..c2[1]], &r1[c3[0]..c3[1]] ]{
            for nuc in km{  
                if *nuc ==b'N'{
                    //let tmp = std::str::from_utf8(km)?;
                    return Err::<u32, &str>( "NnuclError");
                    //Err::<i32, NnuclError<'_>>( NnuclError::new( &tmp ));
                }
            }
        }

        for kmer in needletail::kmer::Kmers::new(&r1[c1[0]..c1[1]], 9 ) { // exactly 1
            let km = Kmer::from(kmer).into_u64();
            
            cell_id += match self.csl1.get( &km ){
                Some(c1) => {
                        println!("to_cellid the c1 {}", c1 );
                        *c1 * max * max
                    },
                None => return Err::<u32, &str>( "no match 1" ),
            };
            //println!("to_cellid 1: {}", cell_id);
        }
        for kmer in needletail::kmer::Kmers::new(&r1[c2[0]..c2[1]], 9 ) {
            
            let km = Kmer::from(kmer).into_u64();
            cell_id +=match self.csl2.get( &km ){
                Some(c2) => {
                        //println!("to_cellid the c2 {}", c2 );
                         *c2 * max 
                    }
                    ,
                None =>  return Err::<u32, &str>( "no match 2" ),

            };
            //println!("to_cellid 2: {}", cell_id);
        }
        for kmer in needletail::kmer::Kmers::new(&r1[c3[0]..c3[1]], 9 ) {
            let km = Kmer::from(kmer).into_u64();
            cell_id +=match self.csl3.get( &km ){
                Some(c3) => {
                    //println!("to_cellid the c3 {}", c3 );
                    *c3
                },
                None =>  return Err::<u32, &str>( "no match 3" ),
            };
            //println!("to_cellid 3: {}", cell_id);
        }

        cell_id += 1;
        //println!("to_cellid final id {}", cell_id);
        return Ok(cell_id)
    }

    pub fn to_sequence(&self, index:u32) -> Vec<&[u8; 9]>{
        let mut idx: u32 = index - 1;
        let code1 = ((idx / (96 * 96)) as f64).floor() as u32;
        idx = idx - code1 * (96 * 96);
        let code2 = ((idx / 96) as f64).floor() as u32;
        idx = idx - code2 * 96;
        println!("index {} -> I translated to the ids {}, {}, {}", index, code1, code2, idx );
        let ret = vec![ self.c1s[code1 as usize], self.c2s[code2 as usize], self.c3s[idx as usize]];
        return ret;
    }

}


// macro_rules! to_cellid {
//    ($r1: expr, $r2: expr) => {
//       //1-9,22-30,44-52
//       $r1.to_cellid($r1, $r2, [0,8], [21,29], [43,51]  )
//    };
// }

#[cfg(test)]
mod tests {
    #[test]
    fn getsamples() {
        let cells = crate::CellIds::new();

        let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
        let mut id:u32 = 1;
        let mut exp= ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        
        let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
        assert_eq!( cells.to_sequence( exp ), exp2 );
        // 3, 3, 3
        primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
        id = 3;
        exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
        exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );

        // and the last one
        primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
        id = 96;
        exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
        exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
        assert_eq!( 884735+1 , exp);
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );        
    }
}



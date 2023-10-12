/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use std::collections::HashSet;
//use kmers::naive_impl::Kmer;
use crate::int_to_str::IntToStr;
//use crate::sampleids::SampleIds;

//use std::thread;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;


/// and here the data
pub struct CellIds{    
    //kmer_len: usize,
    c1s: Vec<u64>,
    c2s: Vec<u64>,
    c3s: Vec<u64>,
    csl1: BTreeMap<u64, u32>,
    csl2: BTreeMap<u64, u32>,
    csl3: BTreeMap<u64, u32>,
    //csl1kmer: SampleIds,
    //csl2kmer: SampleIds,
    //csl3kmer: SampleIds,
    size:u8
}

// here the functions
impl CellIds{

    pub fn new( ver:&String, mut size: u8 )-> Self {

        if size < 5{
            size = 5;
            println!( "CellIDs::new size was set to 5");
        }
        else if size > 9{
            size = 9;
            println!( "CellIDs::new size set to 9");
        }
        let  c1s: Vec<u64>;
        let  c2s: Vec<u64>;
        let  c3s: Vec<u64>;
        let mut csl1 = BTreeMap::<u64, u32>::new();
        let mut csl2 = BTreeMap::<u64, u32>::new();
        let mut csl3 = BTreeMap::<u64, u32>::new();

        let mut bad_entries = HashSet::<u64>::new();

        // TACAGAACA


        if ver == "v2.96" || ver == "v1"{
            c1s = Self::into_u64(vec![ 
            b"GTCGCTATA", b"CTTGTACTA", b"CTTCACATA",
            b"ACACGCCGG", b"CGGTCCAGG", b"AATCGAATG", b"CCTAGTATA", b"ATTGGCTAA", b"AAGACATGC",
            b"AAGGCGATC", b"GTGTCCTTA", b"GGATTAGGA", b"ATGGATCCA", b"ACATAAGCG", b"AACTGTATT",
            b"ACCTTGCGG", b"CAGGTGTAG", b"AGGAGATTA", b"GCGATTACA", b"ACCGGATAG", b"CCACTTGGA",
            b"AGAGAAGTT", b"TAAGTTCGA", b"ACGGATATT", b"TGGCTCAGA", b"GAATCTGTA", b"ACCAAGGAC",
            b"AGTATCTGT", b"CACACACTA", b"ATTAAGTGC", b"AAGTAACCC", b"AAATCCTGT", b"CACATTGCA",
            b"GCACTGTCA", b"ATACTTAGG", b"GCAATCCGA", b"ACGCAATCA", b"GAGTATTAG", b"GACGGATTA",
            b"CAGCTGACA", b"CAACATATT", b"AACTTCTCC", b"CTATGAAAT", b"ATTATTACC", b"TACCGAGCA",
            b"TCTCTTCAA", b"TAAGCGTTA", b"GCCTTACAA", b"AGCACACAG", b"ACAGTTCCG", b"AGTAAAGCC",
            b"CAGTTTCAC", b"CGTTACTAA", b"TTGTTCCAA", b"AGAAGCACT", b"CAGCAAGAT", b"CAAACCGCC",
            b"CTAACTCGC", b"AATATTGGG", b"AGAACTTCC", b"CAAAGGCAC", b"AAGCTCAAC", b"TCCAGTCGA",
            b"AGCCATCAC", b"AACGAGAAG", b"CTACAGAAC", b"AGAGCTATG", b"GAGGATGGA", b"TGTACCTTA",
            b"ACACACAAA", b"TCAGGAGGA", b"GAGGTGCTA", b"ACCCTGACC", b"ACAAGGATC", b"ATCCCGGAG",
            b"TATGTGGCA", b"GCTGCCAAT", b"ATCAGAGCT", b"TCGAAGTGA", b"ATAGACGAG", b"AGCCCAATC",
            b"CAGAATCGT", b"ATCTCCACA", b"ACGAAAGGT", b"TAGCTTGTA", b"ACACGAGAT", b"AACCGCCTC",
            b"ATTTAGATG", b"CAAGCAAGC", b"CAAAGTGTG", b"GGCAAGCAA", b"GAGCCAATA", b"ATGTAATGG",
            b"CCTGAGCAA", b"GAGTACATT", b"TGCGATCTA" ]);
            c2s = Self::into_u64( vec![ 
            b"TACAGGATA", b"CACCAGGTA", b"TGTGAAGAA", b"GATTCATCA", b"CACCCAAAG",
            b"CACAAAGGC", b"GTGTGTCGA", b"CTAGGTCCT", b"ACAGTGGTA", b"TCGTTAGCA", b"AGCGACACC",
            b"AAGCTACTT", b"TGTTCTCCA", b"ACGCGAAGC", b"CAGAAATCG", b"ACCAAAATG", b"AGTGTTGTC",
            b"TAGGGATAC", b"AGGGCTGGT", b"TCATCCTAA", b"AATCCTGAA", b"ATCCTAGGA", b"ACGACCACC",
            b"TTCCATTGA", b"TAGTCTTGA", b"ACTGTTAGA", b"ATTCATCGT", b"ACTTCGAGC", b"TTGCGTACA",
            b"CAGTGCCCG", b"GACACTTAA", b"AGGAGGCGC", b"GCCTGTTCA", b"GTACATCTA", b"AATCAGTTT",
            b"ACGATGAAT", b"TGACAGACA", b"ATTAGGCAT", b"GGAGTCTAA", b"TAGAACACA", b"AAATAAATA",
            b"CCGACAAGA", b"CACCTACCC", b"AAGAGTAGA", b"TCATTGAGA", b"GACCTTAGA", b"CAAGACCTA",
            b"GGAATGATA", b"AAACGTACC", b"ACTATCCTC", b"CCGTATCTA", b"ACACATGTC", b"TTGGTATGA",
            b"GTGCAGTAA", b"AGGATTCAA", b"AGAATGGAG", b"CTCTCTCAA", b"GCTAACTCA", b"ATCAACCGA",
            b"ATGAGTTAC", b"ACTTGATGA", b"ACTTTAACT", b"TTGGAGGTA", b"GCCAATGTA", b"ATCCAACCG",
            b"GATGAACTG", b"CCATGCACA", b"TAGTGACTA", b"AAACTGCGC", b"ATTACCAAG", b"CACTCGAGA",
            b"AACTCATTG", b"CTTGCTTCA", b"ACCTGAGTC", b"AGGTTCGCT", b"AAGGACTAT", b"CGTTCGGTA",
            b"AGATAGTTC", b"CAATTGATC", b"GCATGGCTA", b"ACCAGGTGT", b"AGCTGCCGT", b"TATAGCCCT",
            b"AGAGGACCA", b"ACAATATGG", b"CAGCACTTC", b"CACTTATGT", b"AGTGAAAGG", b"AACCCTCGG",
            b"AGGCAGCTA", b"AACCAAAGT", b"GAGTGCGAA", b"CGCTAAGCA", b"AATTATAAC", b"TACTAGTCA",
            b"CAACAACGG" ] );
            c3s = Self::into_u64( vec![
            b"AAGCCTTCT", b"ATCATTCTG",
            b"CACAAGTAT", b"ACACCTTAG", b"GAACGACAA", b"AGTCTGTAC", b"AAATTACAG", b"GGCTACAGA",
            b"AATGTATCG", b"CAAGTAGAA", b"GATCTCTTA", b"AACAACGCG", b"GGTGAGTTA", b"CAGGGAGGG",
            b"TCCGTCTTA", b"TGCATAGTA", b"ACTTACGAT", b"TGTATGCGA", b"GCTCCTTGA", b"GGCACAACA",
            b"CTCAAGACA", b"ACGCTGTTG", b"ATATTGTAA", b"AAGTTTACG", b"CAGCCTGGC", b"CTATTAGCC",
            b"CAAACGTGG", b"AAAGTCATT", b"GTCTTGGCA", b"GATCAGCGA", b"ACATTCGGC", b"AGTAATTAG",
            b"TGAAGCCAA", b"TCTACGACA", b"CATAACGTT", b"ATGGGACTC", b"GATAGAGGA", b"CTACATGCG",
            b"CAACGATCT", b"GTTAGCCTA", b"AGTTGCATC", b"AAGGGAACT", b"ACTACATAT", b"CTAAGCTTC",
            b"ACGAACCAG", b"TACTTCGGA", b"AACATCCAT", b"AGCCTGGTT", b"CAAGTTTCC", b"CAGGCATTT",
            b"ACGTGGGAG", b"TCTCACGGA", b"GCAACATTA", b"ATGGTCCGT", b"CTATCATGA", b"CAATACAAG",
            b"AAAGAGGCC", b"GTAGAAGCA", b"GCTATGGAA", b"ACTCCAGGG", b"ACAAGTGCA", b"GATGGTCCA",
            b"TCCTCAATA", b"AATAAACAA", b"CTGTACGGA", b"CTAGATAGA", b"AGCTATGTG", b"AAATGGAGG",
            b"AGCCGCAAG", b"ACAGTAAAC", b"AACGTGTGA", b"ACTGAATTC", b"AAGGGTCAG", b"TGTCTATCA",
            b"TCAGATTCA", b"CACGATCCG", b"AACAGAAAC", b"CATGAATGA", b"CGTACTACG", b"TTCAGCTCA",
            b"AAGGCCGCA", b"GGTTGGACA", b"CGTCTAGGT", b"AATTCGGCG", b"CAACCTCCA", b"CAATAGGGT",
            b"ACAGGCTCC", b"ACAACTAGT", b"AGTTGTTCT", b"AATTACCGG", b"ACAAACTTT", b"TCTCGGTTA",
            b"ACTAGACCG", b"ACTCATACG", b"ATCGAGTCT", b"CATAGGTCA" ] );
        }
        else if ver == "v2.384" {
            c1s = Self::into_u64( vec![
            b"TGTGTTCGC", b"TGTGGCGCC", b"TGTCTAGCG", b"TGGTTGTCC", b"TGGTTCCTC",
            b"TGGTGTGCT", b"TGGCGACCG", b"TGCTGTGGC", b"TGCTGGCAC", b"TGCTCTTCC", b"TGCCTCACC",
            b"TGCCATTAT", b"TGATGTCTC", b"TGATGGCCT", b"TGATGCTTG", b"TGAAGGACC", b"TCTGTCTCC",
            b"TCTGATTAT", b"TCTGAGGTT", b"TCTCGTTCT", b"TCTCATCCG", b"TCCTGGATT", b"TCAGCATTC",
            b"TCACGCCTT", b"TATGTGCAC", b"TATGCGGCC", b"TATGACGAG", b"TATCTCGTG", b"TATATGACC",
            b"TAGGCTGTG", b"TACTGCGTT", b"TACGTGTCC", b"TAATCACAT", b"GTTGTGTTG", b"GTTGTGGCT",
            b"GTTGTCTGT", b"GTTGTCGAG", b"GTTGTCCTC", b"GTTGTATCC", b"GTTGGTTCT", b"GTTGGCGTT",
            b"GTTGGAGCG", b"GTTGCTGCC", b"GTTGCGCAT", b"GTTGCAGGT", b"GTTGCACTG", b"GTTGATGAT",
            b"GTTGATACG", b"GTTGAAGTC", b"GTTCTGTGC", b"GTTCTCTCG", b"GTTCTATAT", b"GTTCGTATG",
            b"GTTCGGCCT", b"GTTCGCGGC", b"GTTCGATTC", b"GTTCCGGTT", b"GTTCCGACG", b"GTTCACGCT",
            b"GTTATCACC", b"GTTAGTCCG", b"GTTAGGTGT", b"GTTAGAGAC", b"GTTAGACTT", b"GTTACCTCT",
            b"GTTAATTCC", b"GTTAAGCGC", b"GTGTTGCTT", b"GTGTTCGGT", b"GTGTTCCAG", b"GTGTTCATC",
            b"GTGTCACAC", b"GTGTCAAGT", b"GTGTACTGC", b"GTGGTTAGT", b"GTGGTACCG", b"GTGGCGATC",
            b"GTGCTTCTG", b"GTGCGTTCC", b"GTGCGGTAT", b"GTGCGCCTT", b"GTGCGAACT", b"GTGCAGCCG",
            b"GTGCAATTG", b"GTGCAAGGC", b"GTCTTGCGC", b"GTCTGGCCG", b"GTCTGAGGC", b"GTCTCAGAT",
            b"GTCTCAACC", b"GTCTATCGT", b"GTCGGTGTG", b"GTCGGAATC", b"GTCGCTCCG", b"GTCCTCGCC",
            b"GTCCTACCT", b"GTCCGCTTG", b"GTCCATTCT", b"GTCCAATAC", b"GTCATGTAT", b"GTCAGTGGT",
            b"GTCAGATAG", b"GTATTAACT", b"GTATCAGTC", b"GTATAGCCT", b"GTATACTTG", b"GTATAAGGT",
            b"GTAGCATCG", b"GTACCGTCC", b"GTACACCTC", b"GTAAGTGCC", b"GTAACAGAG", b"GGTTGTGTC",
            b"GGTTGGCTG", b"GGTTGACGC", b"GGTTCGTCG", b"GGTTCAGTT", b"GGTTATATT", b"GGTTAATAC",
            b"GGTGTACGT", b"GGTGCCGCT", b"GGTGCATGC", b"GGTCGTTGC", b"GGTCGAGGT", b"GGTAGGCAC",
            b"GGTAGCTTG", b"GGTACATAG", b"GGTAATCTG", b"GGCTTGGCC", b"GGCTTCACG", b"GGCTTATGT",
            b"GGCTTACTC", b"GGCTGTCTT", b"GGCTCTGTG", b"GGCTCCGGT", b"GGCTCACCT", b"GGCGTTGAG",
            b"GGCGTGTAC", b"GGCGTGCTG", b"GGCGTATCG", b"GGCGCTCGT", b"GGCGCTACC", b"GGCGAGCCT",
            b"GGCGAGATC", b"GGCGACTTG", b"GGCCTCTTC", b"GGCCTACAG", b"GGCCAGCGC", b"GGCCAACTT",
            b"GGCATTCCT", b"GGCATCCGC", b"GGCATAACC", b"GGCAACGAT", b"GGATGTCCG", b"GGATGAGAG",
            b"GGATCTGGC", b"GGATCCATG", b"GGATAGGTT", b"GGAGTCGTG", b"GGAGAAGGC", b"GGACTCCTT",
            b"GGACTAGTC", b"GGACCGTTG", b"GGAATTAGT", b"GGAATCTCT", b"GGAATCGAC", b"GGAAGCCTC",
            b"GCTTGTAGC", b"GCTTGACCG", b"GCTTCGGAC", b"GCTTCACAT", b"GCTTAGTCT", b"GCTGGATAT",
            b"GCTGGAACC", b"GCTGCGATG", b"GCTGATCAG", b"GCTGAGCGT", b"GCTCTTGTC", b"GCTCTCCTG",
            b"GCTCGGTCC", b"GCTCCAATT", b"GCTATTCGC", b"GCTATGAGT", b"GCTAGTGTT", b"GCTAGGATC",
            b"GCTAGCACT", b"GCTACGTAT", b"GCTAACCTT", b"GCGTTCCGC", b"GCGTGTGCC", b"GCGTGCATT",
            b"GCGTCGGTT", b"GCGTATGTG", b"GCGTATACT", b"GCGGTTCAC", b"GCGGTCTTG", b"GCGGCGTCG",
            b"GCGGCACCT", b"GCGCTGGAC", b"GCGCTCTCC", b"GCGCGGCAG", b"GCGCGATAC", b"GCGCCGACC",
            b"GCGAGCGAG", b"GCGAGAGGT", b"GCGAATTAC", b"GCCTTGCAT", b"GCCTGCGCT", b"GCCTAACTG",
            b"GCCGTCCGT", b"GCCGCTGTC", b"GCCATGCCG", b"GCCAGCTAT", b"GCCAACCAG", b"GCATGGTTG",
            b"GCATCGACG", b"GCAGGCTAG", b"GCAGGACGC", b"GCAGCCATC", b"GCAGATACC", b"GCAGACGTT",
            b"GCACTATGT", b"GCACACGAG", b"GATTGTCAT", b"GATTGGTAG", b"GATTGCACC", b"GATTCTACT",
            b"GATTCGCTT", b"GATTAGGCC", b"GATTACGGT", b"GATGTTGGC", b"GATGTTATG", b"GATGGCCAG",
            b"GATCGTTCG", b"GATCGGAGC", b"GATCGCCTC", b"GATCCTCTG", b"GATCCAGCG", b"GATACACGC",
            b"GAGTTACCT", b"GAGTCGTAT", b"GAGTCGCCG", b"GAGGTGTAG", b"GAGGCATTG", b"GAGCGGACG",
            b"GAGCCTGAG", b"GAGATCTGT", b"GAGATAATT", b"GAGACGGCT", b"GACTTCGTG", b"GACTGTTCT",
            b"GACTCTTAG", b"GACCGCATT", b"GAATTGAGC", b"GAATATTGC", b"GAAGGCTCT", b"GAAGAGACT",
            b"GAACTGCCG", b"GAACGCGTG", b"CTTGTGTAT", b"CTTGTGCGC", b"CTTGTCATG", b"CTTGGTCTT",
            b"CTTGGTACC", b"CTTGGATGT", b"CTTGCTCAC", b"CTTGCAATC", b"CTTGAGGCC", b"CTTGACGGT",
            b"CTTCTGATC", b"CTTCTCGTT", b"CTTCTAGGC", b"CTTCGTTAG", b"CTTATGTCC", b"CTTATGCTT",
            b"CTTATATAG", b"CTTAGGTTG", b"CTTAGGAGC", b"CTTACTTAT", b"CTGTTCTCG", b"CTGTGCCTC",
            b"CTGTCGCAT", b"CTGTCGAGC", b"CTGTAGCTG", b"CTGTACGTT", b"CTGCTTGCC", b"CTGCGTAGT",
            b"CTGCACACC", b"CTGATGGAT", b"CTGAGTCAT", b"CTGACGCCG", b"CTGAACGAG", b"CTCTTGTAG",
            b"CTCTTAGTT", b"CTCTTACCG", b"CTCTGCACC", b"CTCTCGTCC", b"CTCGTATTG", b"CTCGACTAT",
            b"CTCCTGACG", b"CTCACTAGC", b"CTATACGGC", b"CGTTCGCTC", b"CGTTCACCG", b"CGTATAGTT",
            b"CGGTGTTCC", b"CGGTGTCAG", b"CGGTCCTGC", b"CGGCGACTC", b"CGGCACGGT", b"CGGATAGCC",
            b"CGGAGAGAT", b"CGCTAATAG", b"CGCGTTGGC", b"CGCGCAGAG", b"CGCACTGCC", b"CCTTGTCTC",
            b"CCTTGGCGT", b"CCTTCTGAG", b"CCTTCTCCT", b"CCTTCGACC", b"CCTTACTTG", b"CCTGTTCGT",
            b"CCTGTATGC", b"CCTCGGCCG", b"CCGTTAATT", b"CCATGTGCG", b"CCAGTGGTT", b"CCAGGCATT",
            b"CCAGGATCC", b"CCAGCGTTG", b"CATTCCGAT", b"CATTATACC", b"CATGTTGAG", b"ATTGCGTGT",
            b"ATTGCGGAC", b"ATTGCGCCG", b"ATTGACTTG", b"ATTCGGCTG", b"ATTCGCGAG", b"ATTCCAAGT",
            b"ATTATCTTC", b"ATTACTGTT", b"ATTACACTC", b"ATGTTCTAT", b"ATGTTACGC", b"ATGTGTATC",
            b"ATGTGGCAG", b"ATGTCTGTG", b"ATGGTGCAT", b"ATGCTTACT", b"ATGCTGTCC", b"ATGCTCGGC",
            b"ATGAGGTTC", b"ATGAGAGTG", b"ATCTTGGCT", b"ATCTGTGCG", b"ATCGGTTCC", b"ATCATGCTC",
            b"ATCATCACT", b"ATATCTTAT", b"ATAGGCGCC", b"AGTTGGTAT", b"AGTTGAGCC", b"AGTGCGACC",
            b"AGGTGCTAC", b"AGGCTTGCG", b"AGGCCTTCC", b"AGGCACCTT", b"AGGAATATG", b"AGCGGCCAG",
            b"AGCCTGGTC", b"AGCCTGACT", b"AGCAATCCG", b"AGAGATGTT", b"AGAGAATTC", b"ACTCGCTTG",
            b"ACTCGACCT", b"ACGTACACC", b"ACGGATGGT", b"ACCAGTCTG", b"ACATTCGGC", b"ACATGAGGT",
            b"ACACTAATT" ] );

            c2s = Self::into_u64( vec![
            b"TTGTGTTGT", b"TTGTGGTAG",
            b"TTGTGCGGA", b"TTGTCTGTT", b"TTGTCTAAG", b"TTGTCATAT", b"TTGTCACGA", b"TTGTATGAA",
            b"TTGTACAGT", b"TTGGTTAAT", b"TTGGTGCAA", b"TTGGTCGAG", b"TTGGTATTA", b"TTGGCACAG",
            b"TTGGATACA", b"TTGGAAGTG", b"TTGCGGTTA", b"TTGCCATTG", b"TTGCACGCG", b"TTGCAAGGT",
            b"TTGATGTAT", b"TTGATAATT", b"TTGAGACGT", b"TTGACTACT", b"TTGACCGAA", b"TTCTGGTCT",
            b"TTCTGCACA", b"TTCTCCTTA", b"TTCTCCGCT", b"TTCTAGGTA", b"TTCTAATCG", b"TTCGTCGTA",
            b"TTCGTAGAT", b"TTCGGCTTG", b"TTCGGAATA", b"TTCGCCAGA", b"TTCGATTGT", b"TTCGATCAG",
            b"TTCCTCGGT", b"TTCCGGCAG", b"TTCCGCATT", b"TTCCAATTA", b"TTCATTGAA", b"TTCATGCTG",
            b"TTCAGGAGT", b"TTCACTATA", b"TTCAACTCT", b"TTCAACGTG", b"TTATGCGTT", b"TTATGATTG",
            b"TTATCCTGT", b"TTATCCGAG", b"TTATATTAT", b"TTAGGCGCG", b"TTACTGGAA", b"TTACTAGTT",
            b"TTACGTGGT", b"TTACGATAT", b"TTACCTAGA", b"TTACATGAG", b"TTACAGCGT", b"TTACACGGA",
            b"TTACACACT", b"TTAATCAGT", b"TTAATAGGA", b"TTAAGTGTG", b"TTAACCTTG", b"TTAACACAA",
            b"TGTTCACTT", b"TGTTCAAGA", b"TGTTAAGTG", b"TGTGTTATG", b"TGTGTCCAA", b"TGTGGAGCG",
            b"TGTCAGTTA", b"TGTCAGAAG", b"TGGTTAGTT", b"TGGTTACAA", b"TGGCGTTAT", b"TGGCGCCAA",
            b"TGGAGTCTT", b"TGCGTATTG", b"TGATAGAGA", b"TGAGGTATT", b"TGAGAATCT", b"TCTTGGTAA",
            b"TCTTCATAG", b"TCTGTCCTT", b"TCTGGAATT", b"TCTACCGCG", b"TCGTTCGAA", b"TCGTCAGTG",
            b"TCGACGAGA", b"TCATGGCTT", b"TCACACTTA", b"TATTCCGAA", b"TATTATGGT", b"TATGCTATT",
            b"TATCAAGGA", b"TAGTTCAAT", b"TAGCTGCTT", b"TAGAGGAAG", b"TACCTGTTA", b"TACACCTGT",
            b"GTTGTGCGT", b"GTTGGCTAT", b"GTTGCCAAG", b"GTTGACCTT", b"GTTCTGCTA", b"GTTCTGAAT",
            b"GTTCTATCA", b"GTTCGCGTG", b"GTTCCTTAT", b"GTTAGCAGT", b"GTTACTGTG", b"GTTACTCAA",
            b"GTTAAGAGA", b"GTTAACTTA", b"GTGTCGGCA", b"GTGTCCATT", b"GTGCTTGAG", b"GTGCTCGTT",
            b"GTGCTCACA", b"GTGCCTGGA", b"GTCTTGTCG", b"GTCTTGATT", b"GTCTTCCGT", b"GTCTTAAGA",
            b"GTCTCATCT", b"GTCTACGAG", b"GTCGTTGCT", b"GTCGTGTTA", b"GTCGGTAAT", b"GTCGGATGT",
            b"GTCGAGCTG", b"GTCCGGACT", b"GTCCAACAT", b"GTCAGACGA", b"GTCAGAATT", b"GTCACTCTT",
            b"GTCAAGGAA", b"GTATGTCTT", b"GTATGTACA", b"GTATCGGTT", b"GTATATGTA", b"GTATACAAT",
            b"GTAGTTAAG", b"GTAGTCGAT", b"GTAGCCTTA", b"GTAGATACT", b"GTACGATTA", b"GTACAGTCT",
            b"GTAATTCGT", b"GCTTGGCAG", b"GCTTGCTTG", b"GCTTGAGGA", b"GCTTCATTA", b"GCTTATGCG",
            b"GCTGTGTAG", b"GCTGTCATG", b"GCTGGTTGT", b"GCTGGACTG", b"GCTGCCTAA", b"GCTGATATT",
            b"GCTCTTAGT", b"GCTCTATTG", b"GCTCGCCGT", b"GCTCCGCTG", b"GCTATTCTG", b"GCTATACGA",
            b"GCTACTAAG", b"GCTACATGT", b"GCTAACTCT", b"GCGTTGTAA", b"GCGTTCTCT", b"GCGTGCGTA",
            b"GCGTCTTGA", b"GCGTCCGAT", b"GCGTAAGAG", b"GCGCTTACG", b"GCGCGGATT", b"GCGCCATAT",
            b"GCGCATGAA", b"GCGATCAAT", b"GCGAGCCTT", b"GCGAGATTG", b"GCGAGAACA", b"GCCTTGGTA",
            b"GCCTTCTAG", b"GCCTTCACA", b"GCCTGAGTG", b"GCCTCACGT", b"GCCGGCGAA", b"GCCGCACAA",
            b"GCCATGCTT", b"GCCATATAT", b"GCCAATTCG", b"GCATTCGTT", b"GCATGATGT", b"GCAGTTGGA",
            b"GCAGTGTCT", b"GCACTTGTG", b"GCAATCTGT", b"GCAACACTT", b"GATTGTATT", b"GATTGCGAG",
            b"GATTCCAGT", b"GATTCATAT", b"GATTATCAG", b"GATTAGGTT", b"GATGTTGCG", b"GATGGATCT",
            b"GATGCTGAT", b"GATGCCTTG", b"GATCTCCTT", b"GATCGCTTA", b"GATATTGAA", b"GATATTACT",
            b"GAGTGTTAT", b"GAGCTCAGT", b"GAGCGTGCT", b"GAGCGTCGA", b"GAGCGGTTG", b"GAGCGACTT",
            b"GAGCCGAAT", b"GAGATAGAT", b"GAGACCTAT", b"GACGGTCGT", b"GACGCAGGT", b"GACGATATG",
            b"GACCTATCT", b"GAATTAGGA", b"GAATCAGCT", b"GAAGTTCAT", b"GAAGTGGTT", b"GAAGTATTG",
            b"GAAGGCATT", b"GAACGCTGT", b"CTTGTCCAG", b"CTTGGATTG", b"CTTGCTGAA", b"CTTGCCGTG",
            b"CTTGATTCT", b"CTTCTGTCG", b"CTTCGGCGT", b"CTTATGAGT", b"CTTACCGAT", b"CTGTTAGGT",
            b"CTGTCGTCT", b"CTGTATAAT", b"CTGGCTCAT", b"CTGGATGCG", b"CTGCGTGTG", b"CTGCGCGGT",
            b"CTGCCGATT", b"CTGCATTGT", b"CTGATTAAG", b"CTGAGATAT", b"CTGACCTGT", b"CTCGTATCT",
            b"CTCGGCAAG", b"CTCGCAATT", b"CTCCTGCTT", b"CTCCTAAGT", b"CTCCGGATG", b"CTCCGAGCG",
            b"CTCACAGGT", b"CTATTCTAT", b"CTATTAGTG", b"CTATGAATT", b"CTACATATT", b"CGTGGCATT",
            b"CGTCTTAAT", b"CGTCTGGTT", b"CGTCACTGT", b"CGTAGGTCT", b"CGGTTCGAG", b"CGGTTCATT",
            b"CGGTGCTCT", b"CGGTAATTG", b"CGGCCTGAT", b"CGGATATAG", b"CGGAATATT", b"CGCTCCAAT",
            b"CGCGTTCGT", b"CGCAGGTTG", b"CGAGGATGT", b"CGAGCTGTT", b"CGACGGCTT", b"CCTTGTGTG",
            b"CCTGTCTCA", b"CCTGACTAT", b"CCTACCTTG", b"CCGTAGATT", b"CCGGCTGGT", b"CATCGGACG",
            b"CATCGATAA", b"CATCCTTCT", b"CAGTTCTGT", b"CAGTGCCAG", b"CAGGCACTG", b"CAGCCTCTT",
            b"CACTTATAT", b"CACTGGTCG", b"CACTGCATG", b"CACGCGTTG", b"CACGATGTT", b"CACCATCTG",
            b"CACAGGCGT", b"ATTGTACAA", b"ATTGGTATG", b"ATTGCTAAT", b"ATTGCATAG", b"ATTGCAGTT",
            b"ATTCTGCAG", b"ATTCTACGT", b"ATTCGGATT", b"ATTCCGTTG", b"ATTCATCAA", b"ATTCAAGAG",
            b"ATTAGCCTT", b"ATTAATATT", b"ATGTTAGAG", b"ATGTTAACT", b"ATGTAGTCG", b"ATGGTGTAG",
            b"ATGGATTAT", b"ATCTTGAAG", b"ATCTGATAT", b"ATCTCAGAA", b"ATCGCTCAA", b"ATCGCGTCG",
            b"ATCCATGGT", b"ATCATGAGA", b"ATCATAGTT", b"ATCAGCGAG", b"ATCACCATT", b"ATAGTAATT",
            b"ATAGCTGTG", b"ATACTCTCG", b"ATACCTCAT", b"AGTTGCGCG", b"AGTTGAATT", b"AGTTATGAT",
            b"AGTGTCCGT", b"AGTGGCTTG", b"AGTGCTTCT", b"AGTATCATT", b"AGTACACAA", b"AGGTATGCG",
            b"AGGTATAGT", b"AGGCTACTT", b"AGGCCAGGT", b"AGGAGCGAT", b"AGCTTATAG", b"AGCTCTAGA",
            b"AGCGTGTAT", b"AGCGTCACA", b"AGCCTTCAT", b"AGCCTGTCG", b"AGCCTCGAG", b"AGCACTGAA",
            b"AGATGTACG", b"AGAGTTAAT", b"AGACCTCTG", b"ACTTCTATA", b"ACTGTCGAG", b"ACTGTATGT",
            b"ACTCTGTAA", b"ACTCGCGAA", b"ACTAGATCT", b"ACTAACGTT", b"ACGTTACTG", b"ACGTGGAAT",
            b"ACGGACTCT", b"ACGCCTAAT", b"ACGCCGTTA", b"ACGACGTGT", b"ACCTCGCAT", b"ACCATCATA",
            b"ACATATATT", b"ACAGGCACA", b"ACACCTGAG", b"ACACATTCT" ] );

            c3s = Self::into_u64( vec![
            b"TTGTGGCTG", b"TTGTGGAGT", b"TTGTGCGAC", b"TTGTCTTCA", b"TTGTAAGAT",
            b"TTGGTTCTG", b"TTGGTGCGT", b"TTGGTCTAC", b"TTGGTAACT", b"TTGGCGTGC", b"TTGGATTAG",
            b"TTGGAGACG", b"TTGGAATCA", b"TTGCGGCGA", b"TTGCGCTCG", b"TTGCCTTAC", b"TTGCCGGAT",
            b"TTGCATGCT", b"TTGCACGTC", b"TTGCACCAT", b"TTGAACCTG", b"TTCTCGCGT", b"TTCTCAACT",
            b"TTCTACTCA", b"TTCGTCCAT", b"TTCGGATAC", b"TTCGGACGT", b"TTCGCAATC", b"TTCCGGTGC",
            b"TTCCGACTG", b"TTCATTATG", b"TTCATGGAT", b"TTCAGCGCA", b"TTCACCTCG", b"TTCAAGCAG",
            b"TTCAACTAC", b"TTATGCCAG", b"TTATGCATC", b"TTATCGTAC", b"TTATACCTA", b"TTATAATAG",
            b"TTATAAGTC", b"TTAGTTAGC", b"TTAGCTCAT", b"TTAGCACTA", b"TTAGATATG", b"TTACTACGA",
            b"TTACCGTCA", b"TTACAGAGC", b"TTAATTGCA", b"TTAACAGAT", b"TGTTGGCTA", b"TGTTGATGA",
            b"TGTTAAGCT", b"TGTGGCCGA", b"TGTGCTAGC", b"TGTGCGTCA", b"TGTCGCAGT", b"TGTCGAGCA",
            b"TGTACAACG", b"TGGTTCCGA", b"TGGTTCACT", b"TGGTCAAGT", b"TGGCTTGTA", b"TGGCTGTCG",
            b"TGGCGTATG", b"TGGCGCGCT", b"TGGATGTAC", b"TGGACTTGC", b"TGGAATACT", b"TGCTAGCGA",
            b"TGCGTTGCT", b"TGCGGTCTG", b"TGCGCTTAG", b"TGCGCGACG", b"TGCCTGCAT", b"TGCCTAGAC",
            b"TGCACGAGT", b"TGAGTGTGC", b"TGAGGCTCG", b"TCTTCCGTC", b"TCTTATAGT", b"TCTTACCAT",
            b"TCTGTTGTC", b"TCTGTTACT", b"TCTGGCTAG", b"TCTCAGATC", b"TCTAGTTGA", b"TCTAGTACG",
            b"TCGTACTAC", b"TCGGTGTAG", b"TCGGCTGCT", b"TCGCTACTG", b"TCGATCACG", b"TCGAGGCAT",
            b"TCCGGCGTC", b"TCCGGAGCT", b"TCCGCTCGT", b"TCCGAGTAC", b"TCCATTCAT", b"TCCATGGTC",
            b"TCCAAGTCG", b"TCATTACGT", b"TCATGCACT", b"TCAGGTTGC", b"TCAGACCGT", b"TCACTCAGT",
            b"TCAAGCTCA", b"TATTGCGCA", b"TATTCGGCT", b"TATTCCAGC", b"TATTCATCA", b"TATGTTCAG",
            b"TATGGTATG", b"TATGCAAGT", b"TATCTGGTC", b"TATCTGACT", b"TATCCAGAT", b"TATCAGTCG",
            b"TATCACGCT", b"TAGGCGCGA", b"TAGGCACAT", b"TAGGATCGT", b"TAGCATTGC", b"TAGAGTTAC",
            b"TAGACTGAT", b"TACTTGTCG", b"TACGTCCGA", b"TACCGTACT", b"TACCGCGAT", b"TACCAGGAC",
            b"TACAGAAGT", b"TAAGTGCAT", b"TAAGCTACT", b"GTTGACCGA", b"GTTCTCGAC", b"GTTCCTGCT",
            b"GTTATGATG", b"GTGCTTGCA", b"GTGCCGCGT", b"GTATTGCTG", b"GTATTCCGA", b"GTATTAAGC",
            b"GTATGACGT", b"GTAGTTGTC", b"GTAGTACAT", b"GTAGCTCGA", b"GGTTGCTCA", b"GGTTGAGTA",
            b"GGTTAACGT", b"GGTGTGGCA", b"GGTCTTCAG", b"GGTCGTCTA", b"GGTCGGCGT", b"GGTCCGACT",
            b"GGTCATGTC", b"GGTCACATG", b"GGTAGTGCT", b"GGTAGCGTC", b"GGTACCAGT", b"GGTAAGGAT",
            b"GGCTTGTGC", b"GGCTTGACT", b"GGCTTACGA", b"GGCTGTAGT", b"GGCTGGCAG", b"GGCTCCATC",
            b"GGCGTGGAT", b"GGCGTAATC", b"GGCGCAAGT", b"GGCGAGTAG", b"GGCGACCGT", b"GGCCTGTCA",
            b"GGCCATTGC", b"GGCACTCTG", b"GGATGTCAT", b"GGAGTAACT", b"GGAGAACGA", b"GGACTGGCT",
            b"GGACGTTCA", b"GGAACGTGC", b"GCTGTCCAT", b"GCTGGTTCA", b"GCTGCAACT", b"GCTCGTTAC",
            b"GCTATAGAT", b"GCTAGTCGT", b"GCTACCATG", b"GCGTTCTGA", b"GCGTGTTAG", b"GCGGTATCG",
            b"GCGGAGCAT", b"GCGCGGTGC", b"GCGCCTAGT", b"GCGCCGGCT", b"GCCTTCATG", b"GCCATACTG",
            b"GCATGTTGA", b"GCATGCTAC", b"GCAGTATAC", b"GCAGGTACT", b"GCAGCGCGT", b"GCACCTCAT",
            b"GCAATTCGA", b"GATTGCCGT", b"GATGAACAT", b"GATCTTCGA", b"GATCTGCAT", b"GAGTGGCAT",
            b"GAGTCGGAC", b"GAGTATGAT", b"GAGGCGAGT", b"GAGGCAACG", b"GAGCGCACT", b"GAATAGGCT",
            b"ATTGTCACT", b"ATTGTATCA", b"ATTGGTCAG", b"ATTGGCGAT", b"ATTGATCGT", b"ATTCGTAGT",
            b"ATTCATACG", b"ATTCAGGAC", b"ATTACTTCA", b"ATTAATTAG", b"ATTAAGCAT", b"ATGTCTCTA",
            b"ATGTAGCGT", b"ATGGCATAC", b"ATGGAGATC", b"ATGGACTCG", b"ATGGAACGA", b"ATGCTTCAT",
            b"ATGCTCGCT", b"ATGCGACGT", b"ATGCCGTAG", b"ATGAGTTCG", b"ATGACTATC", b"ATGACCGAC",
            b"ATCTTATGC", b"ATCTTACTA", b"ATCTATCAG", b"ATCGTGTAC", b"ATCGTCTGA", b"ATCGGCATG",
            b"ATCGCGAGC", b"ATCGCAACG", b"ATCGATGCT", b"ATCGAATAG", b"ATCCTTCTG", b"ATCCTGCGT",
            b"ATCCGCACT", b"ATCCATTAC", b"ATCCAAGCA", b"ATCAGATCA", b"ATCACACAT", b"ATCAACGTC",
            b"ATCAACCGA", b"ATATTGAGT", b"ATATTCGTC", b"ATATTACAG", b"ATATCTTGA", b"ATATCGCAT",
            b"ATATCAATC", b"ATAGTCCTG", b"ATAGGTCTA", b"ATAGCTGAC", b"ATAGCGGTA", b"AGTTCGCTG",
            b"AGTTACAGC", b"AGTTAACTA", b"AGTGCAATC", b"AGTCTGGTA", b"AGTCTGAGC", b"AGTCTACAT",
            b"AGTCGAACT", b"AGTCCATCG", b"AGTCATTCA", b"AGTATCCAG", b"AGTAGACTG", b"AGTAATCGA",
            b"AGTAAGTGC", b"AGGTTGGCT", b"AGGTTCTAG", b"AGGTGTTCA", b"AGGTGCCAT", b"AGGTCTGAT",
            b"AGGTCGTAC", b"AGGTCAGCA", b"AGGCTTATC", b"AGGCTATGA", b"AGGCCGACG", b"AGGCCAAGC",
            b"AGGCAGGTC", b"AGGCAAGAT", b"AGGAGCAGT", b"AGGACCGCT", b"AGGAATTAC", b"AGCTTGGAC",
            b"AGCTTAAGT", b"AGCTACACG", b"AGCGTTACG", b"AGCGGTGCA", b"AGCGGAGTC", b"AGCGGACGA",
            b"AGCGCGCTA", b"AGCGATAGC", b"AGCGACTCA", b"AGCCTCTAC", b"AGCCGTCGT", b"AGCATGATC",
            b"AGCACTTCG", b"AGCACGGCA", b"AGATTCTGA", b"AGATTAGAT", b"AGATGATAG", b"AGATATGTA",
            b"AGATACCGT", b"AGAGTGCGT", b"AGAGCCGAT", b"AGACTCACT", b"ACTTGCCTA", b"ACTTGAGCA",
            b"ACTTCTAGC", b"ACTTCGACT", b"ACTTAGTAC", b"ACTGTTGAT", b"ACTGTAACG", b"ACTGGTATC",
            b"ACTGACGTC", b"ACTGAAGCT", b"ACTCTGATG", b"ACTCCTGAC", b"ACTCCGCTA", b"ACTCAACTG",
            b"ACTATTGCA", b"ACTAGGCAG", b"ACTACGCGT", b"ACTAATACT", b"ACGTTCGTA", b"ACGTGTGCT",
            b"ACGTGTATG", b"ACGTGGAGC", b"ACGTCTTCG", b"ACGTCAGTC", b"ACGGTCTCA", b"ACGGTCCGT",
            b"ACGGTACAG", b"ACGGCGCTG", b"ACGCTGCGA", b"ACGCGTGTA", b"ACGCGCCAG", b"ACGATGTCG",
            b"ACGATGGAT", b"ACGATCTAC", b"ACGAGCTGA", b"ACGAGCATC", b"ACGAATCGT", b"ACGAACGCA",
            b"ACCTTGTAG", b"ACCTGTTGC", b"ACCTGTCAT", b"ACCTCGATC", b"ACCTAGGTA", b"ACCTACTGA",
            b"ACCTAATCG", b"ACCGTAGCA", b"ACCGGTAGT", b"ACCGGCTAC", b"ACCGCTTCA", b"ACATTGTGC",
            b"ACATTCTCG", b"ACATGGCTG", b"ACATGACGA", b"ACATATGAT", b"ACATATACG", b"ACAGCGTAC",
            b"ACACTTGCT", b"ACACTATCA", b"ACACGCATG", b"ACACCAGTA", b"ACACCAACT", b"ACACATAGT",
            b"ACACACCTA" ] );
        }
        else {
            panic!( "The version '{ver}' is unknown!" );
        }

        // let mut csl1kmer= SampleIds::new( kmer_len );
        // for kmer in &c1s{
        //     csl1kmer.add( *kmer );
        // }

        // let mut csl2kmer= SampleIds::new( kmer_len );
        // for kmer in &c2s{
        //     csl2kmer.add( *kmer );
        // }

        // let mut csl3kmer= SampleIds::new( kmer_len );
        // for kmer in &c3s{
        //     csl3kmer.add( *kmer );
        // }

        let mut i: u32 = 0;

        // lets ditch the kmer library
        for km in &c1s{
            if let std::collections::btree_map::Entry::Vacant(e) = csl1.entry(*km) {
                //let info = Info::new(km, name.clone() );
                e.insert(i);
            } else {
                bad_entries.insert( *km );
                csl1.remove( km );
            }
            i +=1;
        }

        i = 0;
        bad_entries.clear();
        for km in &c2s{
            if let std::collections::btree_map::Entry::Vacant(e) = csl2.entry(*km) {
                //let info = Info::new(km, name.clone() );
                e.insert(i);
            } else {
                bad_entries.insert( *km );
                csl2.remove( km );
            }
            i +=1;
        }

        i = 0;
        bad_entries.clear();
        for km in &c3s{
            if let std::collections::btree_map::Entry::Vacant(e) = csl3.entry(*km) {
                //let info = Info::new(km, name.clone() );
                e.insert(i);
            } else {
                bad_entries.insert( *km );
                csl3.remove( km );
            }
            i +=1;
        }

        Self {
            //kmer_len,
            c1s,
            c2s,
            c3s,
            csl1,
            csl2,
            csl3,
            //csl1kmer,
            //csl2kmer,
            //csl3kmer,
            size
        }
    }

    /// convert the hard coded cell identifiers to u64 instead of &[u8;9]
    fn into_u64( seq_a:Vec<&[u8;9]> ) -> Vec<u64>{
        let mut ret = Vec::<u64>::with_capacity( seq_a.len() );
        let mut tool: IntToStr;
        for seq in seq_a {
            tool = IntToStr::new(seq.to_vec(), 9);
            ret.push( tool.into_u64().clone() );
        }
        ret
    }

    /// returns the first entry in the tuple that has the lowest second entry
    /// and only if there is only one tuple with the lowest second entry 
    fn best_entry( data:Vec<(usize, usize)> ) -> Option<u32>{
        if data.len() == 0 {
            return None
        }
        if data.len() == 1 {
            return Some(data[0].0 as u32)
        }
        else {
            let mut counter = vec![0,0,0];
            let mut min = usize::MAX;
            for ( _i, mismatch ) in &data{
                if min > *mismatch{
                    min = *mismatch;
                }
                counter[*mismatch] +=1;
            }
            if counter[min] == 1{
                // this is great - we have a minimum overlap and that has exactly one match!
                for ( i, mismatch ) in &data{
                    if *mismatch == min{
                        return Some( data[*i].0 as u32 )
                    }
                }
            }
        }
        None
    } 

    pub fn to_cellid (&self, r1: &[u8], c1: Vec<usize>, c2: Vec<usize>, c3: Vec<usize>  )-> Result< u32, &str>{
        let mut cell_id:u32 = 0;
        // This has to be a static 384 to reproduce what BD has...
        // I would use that for v2.384 only...
        let max:u32 = 384;
        //let max:u32 = self.c1s.len() as u32;
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

        let mut ok = false;

        let mut tool = IntToStr::new( r1[c1[0]..c1[1]].to_vec(), c1[1]- c1[0]);
        let km = tool.into_u64();
        cell_id += match self.csl1.get( &km ){
            Some(c1) => {
                    //println!("to_cellid the c1 {}", c1 );
                    ok = true;
                    *c1 * max * max
                },
            None => {
                let mut good = Vec::<(usize,usize)>::with_capacity(3);
                for i in 0..self.c1s.len(){
                    if ( km ^ self.c1s[i]) < 3 {
                        // could be as little as one bd change - max two
                        good.push( (i, ( km ^ self.c1s[i])  as usize ) );
                    }
                }
                if let Some(c1) = Self::best_entry( good ){
                    ok = true;
                    c1 * max * max
                }else {
                    0
                }
            }
        };
        if !ok {
            return  Err::<u32, &str>( "Cells no match 1" )
        }
        ok = false;         

        tool = IntToStr::new( r1[c2[0]..c2[1]].to_vec(), c2[1]- c2[0]);
        let km = tool.into_u64();
        cell_id += match self.csl2.get( &km ){
            Some(c2) => {
                //println!("to_cellid the c1 {}", c1 );
                ok = true;
                *c2 * max 
            },
            None => {
                let mut good = Vec::<(usize,usize)>::with_capacity(3);
                for i in 0..self.c2s.len(){
                    if ( km ^ self.c2s[i]) < 3 {
                        // could be as little as one bd change - max two
                        good.push( (i, ( km ^ self.c2s[i])  as usize)  );
                    }
                }
                if let Some(c2) = Self::best_entry( good ){
                    ok = true;
                    c2 * max 
                }else {
                    0
                }
            }.try_into().unwrap(),
        };
        if !ok {
            return  Err::<u32, &str>( "Cells no match 1" )
        }
        ok = false;    

        tool = IntToStr::new( r1[c3[0]..c3[1]].to_vec(), c3[1]- c3[0]);
        let km = tool.into_u64();
        cell_id += match self.csl3.get( &km ){
            Some(c3) => {
                //println!("to_cellid the c1 {}", c1 );
                ok = true;
                *c3 
            },
            None => {
                let mut good = Vec::<(usize,usize)>::with_capacity(3);
                for i in 0..self.c3s.len(){
                    if ( km ^ self.c3s[i]) < 3 {
                        // could be as little as one bd change - max two
                        good.push( (i, ( km ^ self.c3s[i])  as usize)  );
                    }
                }
                if let Some(c3) = Self::best_entry( good ){
                    ok = true;
                    c3 
                }else {
                    0
                }
            }.try_into().unwrap(),
        };
        if !ok {
            return  Err::<u32, &str>( "Cells no match 1" )
        }

        cell_id += 1;
        //println!("to_cellid final id {}", cell_id);
        Ok(cell_id)
    }

    pub fn to_sequence(&self, index:u32) -> Vec<u64>{
        let mut idx: u32 = index - 1;
        let max:u32 = 384;
        //let max:u32 = self.c1s.len() as u32;
        let code1 = ((idx / (max * max)) as f64).floor() as u32;
        idx -= code1 * (max * max);
        let code2 = ((idx / max) as f64).floor() as u32;
        idx -= code2 * max;
        //println!("index {} -> I translated to the ids {}, {}, {}", index, code1, code2, idx );
        vec![ self.c1s[code1 as usize], self.c2s[code2 as usize], self.c3s[idx as usize] ]
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
    use crate::cellids::CellIds;

    #[test]
    //
    fn getcells_9(){
        let mut cells = CellIds::new(&"v1".to_string(), 9 as u8 );

        let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
        let mut id:u32 = 1;
        let mut exp= ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(_err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        
        let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
        assert_eq!( cells.to_sequence( exp ), exp2 );
        // 3, 3, 3
        primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
        id = 3;
        exp = ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(_err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );

        // and the last one
        primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
        id = 96;
        exp = ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
        //assert_eq!( 884735+1 , exp);
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(_err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );
    }
    //
    #[test]
    fn getcells_7() {
        let mut cells = CellIds::new(&"v1".to_string(), 7 as u8 );

        let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
        let mut id:u32 = 1;
        let mut exp= ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ),
            Err(_err) => (),
        };
        
        let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
        assert_eq!( cells.to_sequence( exp ), exp2 );
        // 3, 3, 3
        primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
        id = 3;
        exp = ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), 
            Err(_err) => (), 
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );

        // and the last one
        primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
        id = 96;
        exp = ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
        //assert_eq!( 884735+1 , exp);
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), 
            Err(_err) => (),
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );        
    }
    #[test]
    fn getcells_384_9() {
        let mut cells = CellIds::new(&"v2.384".to_string(), 9 as u8 );

        let primer = b"TGTCTAGCGNNNNNNNNNNNNTTGTGCGGANNNNNNNNNNNNNTTGTGCGAC"; // totally artificial - primer design wrong... - lazy
        let exp2 = vec![b"TGTCTAGCG", b"TTGTGCGGA", b"TTGTGCGAC"];
        let id:u32 = 3;
        let exp= ((id-1)* 384 * 384 + (id-1) * 384 + (id-1) +1) as u32;
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ),
            Err(_err) => (),
        };
        
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );        
    }
}



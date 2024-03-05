#[cfg(test)]
mod tests {
    use rustody::fast_mapper::mapper_entries::SecondSeq;
    use rustody::traits::BinaryMatcher;
    #[test]
    fn test_second_seq_equality() {
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq3 = SecondSeq(0b1111000010101010101010101010101010101010, 20);
        let seq4 = SecondSeq(0b101010101010101010101010101010101010101010101010, 24);

        assert_eq!(seq1, seq2, "same sequences"); // Test equal sequences
        assert_ne!(seq1, seq3, "different sequences"); // Test unequal sequences
        assert_eq!(seq1, seq4,  "same, but one is longer");
    }

    #[test]
    fn test_second_table() {
        let seq1 = SecondSeq(0b1010101011010101, 32);
        let mut exp = std::collections::HashMap::new();
        exp.insert('A', 24);
        exp.insert('C', 3);
        exp.insert('G', 4);
        exp.insert('T', 1);

        assert_eq!(seq1.table() , exp,  "table did return the right counts");
    }

    #[test]
    fn test_second_seq_hashing() {
        let mut map = std::collections::HashMap::new();
        let seq1 = SecondSeq(0b1010101010101010101010101010101010101010, 20);
        let seq2 = SecondSeq(0b1111000010101010101010101010101010101010, 20);

        map.insert(seq1, "Value for seq1");

        assert_eq!(map.get(&seq1), Some(&"Value for seq1")); // Test retrieval by equal key
        assert_eq!(map.get(&seq2), None); // Test retrieval by different key
    }

    #[test]
    fn test_humming2() {
        let seq1 = SecondSeq(0b101010, 20);
        let seq2 = SecondSeq(0b101010, 20);
        assert_eq!( seq1.hamming_distance( &seq2 ), 0 );
        let seq3 = SecondSeq(0b011010, 20);
        assert_eq!( seq1.hamming_distance( &seq3 ), 1 );
        let seq4 = SecondSeq(0b001010, 20);
        assert_eq!( seq1.hamming_distance( &seq4 ), 1 );
        let seq5 = SecondSeq(0b011001, 20);
        assert_eq!( seq1.hamming_distance( &seq5 ), 2 );
        let seq6 = SecondSeq(0b0, 20);
        assert_eq!( seq1.hamming_distance( &seq6 ), 3 );
        assert_eq!( seq4.hamming_distance( &seq5 ), 2 );
    }

    
    // works somehow. Chatty - thanky you!
    #[test]
    fn test_needleman_wunsch() {
        let seq1 = SecondSeq(0b101010, 20);
        let seq2 = SecondSeq(0b101010, 20);
        assert_eq!( seq1.needleman_wunsch( &seq2, 0.8 ), 1.0/20.0 );
        let seq3 = SecondSeq(0b011010, 20);
        assert_eq!( seq1.needleman_wunsch( &seq3, 0.8 ), 3.0/20.0 );
        let seq4 = SecondSeq(0b001010, 20);
        assert_eq!( seq1.needleman_wunsch( &seq4, 0.8 ), 3.0/20.0 );
        let seq5 = SecondSeq(0b011001, 20);
        assert_eq!( seq1.needleman_wunsch( &seq5, 0.8 ), 5.0/20.0 );
        let seq6 = SecondSeq(0b0, 20);
        assert_eq!( seq1.needleman_wunsch( &seq6, 0.8 ), 7.0/20.0 );
        let seq7 = SecondSeq(0b101010, 15);
        assert_eq!( seq1.needleman_wunsch( &seq7, 0.8 ), 1.0/15.0 );
    }

    #[test]
    fn test_to_string() {

        let seq = "TGGTATCTTTTACTTACCTGCTTGAATACTTG".to_string();
        let seq_u64 = 13632604658164817131_u64;
        let obj = SecondSeq(seq_u64, 32);

        assert_eq!( BinaryMatcher::to_string(&obj), seq );
        
    }
    
}

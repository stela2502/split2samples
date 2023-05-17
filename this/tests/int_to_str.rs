// IntToString.rs

// Define your library code here

// Include a module for tests

#[cfg(test)]
mod tests {

   use this::int_to_str::IntToStr;

    #[test]
    fn test_u64_to_str(){

        let tool = IntToStr::new( b"CGATATT".to_vec(), 32);

        let num:u64 = tool.into_u64();
        println!("I have this number for the sting 'CGATATT' {num}");

        //let num:u64 = 15561;

        //println!("This is the number I want to convert to a Sequence: {num:b}");
        // 11110011001001
        // T T A T A C G 
        //let tool = IntToStr::new();
        let mut data:String = "".to_string();
        tool.to_string( 7, &mut data );
        assert_eq!( data, "CGATATT".to_string() )
    } 

     #[test]
    fn check_conversion_4bp() {

     let seq = b"AGGC";
     //         C G G A
     //         01101000
     let tool = IntToStr::new( seq.to_vec() , 32);

     assert_eq!( tool.len(),  1 ); 
     //panic!("{:b}", binary[0] );
     //                          G G C 
     assert_eq!( tool, vec![ 0b1101000 ]);

     let mut decoded:String = "".to_string();

     tool.to_string( 15, &mut decoded );
     assert_eq!( decoded, "AGGC" );      
                            
    }



     #[test]
    fn check_conversion_15bp() {
        //          0000111122223333   
     let seq = b"AGGCTTGATAGCGAG";
     let tool = IntToStr::new(seq.to_vec(),32);

     assert_eq!( tool.len(),  4 );

     //panic!("{:b} {:b} {:b} {:b}", binary[0], binary[1], binary[2], binary[3] );
     //                                                A G C A
     //                                              0b00100100  
     //                          G G C     T T G A     T A G C     G A G

     //panic!("the binaries I get: {:b} {:b} {:b} {:b} ", binary[0], binary[1], binary[2], binary[3]);
     //assert_eq!( binary, vec![ 0b1101000, 0b101111, 0b1100011, 0b100010 ]);
     assert_eq!( tool, vec![ 0b100010, 0b1100011, 0b101111, 0b1101000   ]);
     let mut decoded:String = "".to_string();

     tool.to_string(15, &mut decoded );
     assert_eq!( decoded, "AGGCTTGATAGCGAG" );
    }

     #[test]
    fn check_conversion_1bp() {

     let seq = b"C";
     let tool = IntToStr::new( seq.to_vec(), 10 );

     assert_eq!( tool.len(),  1 ); 
     //                                                A G C A
     //                                              0b00100100  
     //                          G G C     T T G A     T A G C     G A G
     assert_eq!( tool, vec![ 0b1 ]);

     let mut decoded:String = "".to_string();

     tool.to_string(1, &mut decoded );
     assert_eq!( decoded, "C" );
    }

    #[test]
    fn check_conversion_oneA() {

     let seq = b"A";
     let tool = IntToStr::new(seq.to_vec(), 32);

     assert_eq!( tool.len(),  1 ); 
     //                                                A G C A
     //                                              0b00100100  
     //                          G G C     T T G A     T A G C     G A G
     assert_eq!( tool, vec![ 0b0 ]);

     let mut decoded:String = "".to_string();

     tool.to_string(1, &mut decoded );
     assert_eq!( decoded, "A" );
    }


    #[test]
    fn check_conversion_4_from_15bp() {

     let seq = b"AGGCCTGTATGA";
     let tool = IntToStr::new( seq.to_vec(), 10);

     assert_eq!( tool.len(),  3 ); 

     //                          G G C 
     assert_eq!( tool[2], 0b1101000 );

     let mut decoded:String = "".to_string();

        tool.to_string(4, &mut decoded );
        assert_eq!( decoded, "AGGC" );
        decoded.clear();
        tool.to_string(3, &mut decoded );
     assert_eq!( decoded, "AGG" );                
    }

    #[test]
    fn check_longer_string() {

     let seq = b"CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTC";
     let tool = IntToStr::new(seq.to_vec(), 32);

     assert_eq!( tool.len(), 12 ); 

     //                       CT G G
     //                       G G T C
     assert_eq!( tool[11], 0b10101101  );

     let mut decoded:String = "".to_string();

        tool.to_string(4,  &mut decoded );
        assert_eq!( decoded, "CTGG" );      
        decoded.clear();
        tool.to_string(tool.len()*4 -3, &mut decoded );
        assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTC" );       

        decoded.clear();
        tool.to_string(tool.len()*4 -6,  &mut decoded );
        assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGG" );               
    }

    #[test]
    fn check_u64_decode() {

     let seq = b"CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTC";
     let mut tool = IntToStr::new(seq.to_vec(), 8);
     let binary = tool.into_u64(  );


     assert_eq!( binary, 7706141192143397037 ); 

     let mut decoded:String = "".to_string();
     tool.to_string( 32, &mut decoded);
     assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGC".to_string() );

     assert_eq!( tool.next(), Some( (24749_u16,55990_u64) ) );
     assert_eq!( tool.next(), Some( (55990_u16,46741_u64) ) );
     assert_eq!( tool.next(), Some( (46741_u16,27377_u64) ) );
     assert_eq!( tool.next(), Some( (27377_u16,58859_u64) ) );
     assert_eq!( tool.next(), Some( (38955_u16,30381_u64) ) );

     assert_eq!( tool.next(), Some( (30381_u16,28069_u64) ) );
     assert_eq!( tool.next(), Some( (28069_u16,55996_u64) ) );
     assert_eq!( tool.next(), Some( (55996_u16,47482_u64) ) );
     assert_eq!( tool.next(), Some( (26122_u16,23979_u64) ) );
     assert_eq!( tool.next(), Some( (23979_u16,7017_u64) ) );

     assert_eq!( tool.next(), Some( (7017_u16,46767_u64) ) );



     //                       CT G G
     //                       G G T C
     // assert_eq!( binary[11], 0b10101101  );

     // let mut decoded:String = "".to_string();

     // tool.decode_vec(4, binary.clone(), &mut decoded );
     // assert_eq!( decoded, "CTGG" );      
        // decoded.clear();
        // tool.decode_vec(binary.len()*4 -3, binary.clone(), &mut decoded );
     // assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTC" );       

     // decoded.clear();
        // tool.decode_vec(binary.len()*4 -6, binary.clone(), &mut decoded );
     // assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGG" );               
    }

    // #[test]
    // fn u128_to_str(){
    //  let x:u128  = 8587300487894951391; // this fits into a u64
    //  let mut data = "".to_string();
    //  let tool = IntToStr::new(b"".to_vec(), 32);
    //  tool.u128_to_str( 32, &x, &mut data);

    //  assert_eq!( data, "TCTCATGAAGTATGACAGCTACAGCCGCTTCT".to_string() );

    //  data.clear();
    //  tool.u128_to_str( 12, &x, &mut data);
    //  assert_eq!( data, "TCTCATGAAGTA".to_string() );


    //  data.clear();
    //  tool.u128_to_str( 9, &x, &mut data);
    //  assert_eq!( data, "TCTCATGAA".to_string() );

    //  data.clear();
    //  tool.u128_to_str( 40, &x, &mut data);
    //  assert_eq!( data, "TCTCATGAAGTATGACAGCTACAGCCGCTTCTAAAAAAAA".to_string() );

    // }



}
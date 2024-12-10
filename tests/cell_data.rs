#[cfg(test)]
mod tests {
    use rustody::singlecelldata::cell_data::CellData;
    use rustody::singlecelldata::IndexedGenes;
    use rustody::singlecelldata::cell_data::GeneUmiHash;

    #[test]
    fn merge_new_approach() {
        // define objects
        let mut cell1 = CellData::new(1);
        let mut cell2 = CellData::new(1);
        let mut index1 = IndexedGenes::empty(None);
        let mut index2 = IndexedGenes::empty(None);

        // create some genes
        let _ = index1.get_gene_id("Gene1");
        let _ = index2.get_gene_id("Gene1");

        let _ = index1.get_gene_id("GeneX");

        let _ = index1.get_gene_id("Gene2");
        let _ = index2.get_gene_id("Gene3");

        let _ = index2.get_gene_id("GeneX");

        assert_eq!( index1.get_all_gene_names(), vec!["Gene1".to_string(), "GeneX".to_string(), "Gene2".to_string()], "Gene names in cell1" );

        assert_eq!( index2.get_all_gene_names(), vec!["Gene1".to_string(), "Gene3".to_string(), "GeneX".to_string()], "Gene names in cell2" );

        // create one read for "Gene1" with UMI 1
        let counts1 = GeneUmiHash( 0, 1_u64 );
        // add the same read into both cells
        cell1.add( counts1.clone() ); // Gene1, 1
        cell2.add( counts1 );         // Gene1 1 -> total of 1 in merge

        // create one read for "Gene1" with UMI 3
        let counts2 = GeneUmiHash( 0, 3_u64 );
        // add that only to cell1
        cell1.add( counts2 );        // Gene1, 3 -> total of 2 in merge

        // create one read for gene 2 (GeneX in 1 and Gene3 in 2) with UMI 3
        let counts3 = GeneUmiHash( 1, 3_u64 );
        cell1.add( counts3.clone() ); // GeneX 3 -> total of 1 in merge
        cell2.add( counts3 );         // Gene3 3 -> total of 1 in merge

        // create one read for gene 2 (Gene2 in 1 and GeneX in 2) with UMI 3
        let counts4 = GeneUmiHash( 2, 3_u64 );
        cell1.add( counts4.clone() ); // Gene2, 3 -> total of 1 in merge
        cell2.add( counts4 );         // GeneX, 3 -> total of 1 in merge

        assert_eq!( cell1.n_umi() ,4, "total counts cell1" );
        assert_eq!( cell2.n_umi() ,3, "total counts cell2" );

        let translation = index1.merge( &index2 );
        assert_eq!( translation, vec![0,3,1], "cell2 genes mapped to cell1 order" );

        assert_eq!( cell1.n_umi_4_gene_id( &0 ), 2, "Gene1 in cell"); //Gene1
        assert_eq!( cell1.n_umi_4_gene_id( &1 ), 1, "GeneX in cell1" ); //GeneX
        assert_eq!( cell1.n_umi_4_gene_id( &2 ), 1, "Gene2 in cell1" ); //Gene2

        assert_eq!( cell2.n_umi_4_gene_id( &0 ), 1, "Gene1 in cell2"); //Gene1
        assert_eq!( cell2.n_umi_4_gene_id( &1 ), 1, "GeneX in cell2" ); //Gene3
        assert_eq!( cell2.n_umi_4_gene_id( &2 ), 1, "Gene3 in cell2" ); //GeneX

        cell1.merge_re_id_genes( &cell2, &translation );

        assert_eq!( index1.get_all_gene_names(), vec!["Gene1".to_string(), "GeneX".to_string(), "Gene2".to_string(), "Gene3".to_string()], 
            "new gene names for the merged index" );

        assert_eq!( cell1.n_umi_4_gene_id( &0 ), 2, "Gene1 counts in merged cell" ); //Gene1
        assert_eq!( cell1.n_umi_4_gene_id( &1 ), 1, "GeneX counts in merged cell" ); //GeneX
        assert_eq!( cell1.n_umi_4_gene_id( &2 ), 1, "Gene2 counts in merged cell" ); //Gene2
        assert_eq!( cell1.n_umi_4_gene_id( &3 ), 1, "Gene3 counts in merged cell" ); //Gene3

    }

}
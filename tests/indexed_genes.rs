#[cfg(test)]
mod tests {
    use rustody::singlecelldata::IndexedGenes;

    #[test]
    fn test_new_approach () {
        let mut genes = IndexedGenes::empty( Some(0) );
        let id1 = genes.get_gene_id("geneA");
        let id2 = genes.get_gene_id("geneB");

        assert_eq!(id1, 0);
        assert_eq!(id2, 1);
        assert_eq!( genes.get_all_gene_names(), vec!["geneA".to_string(), "geneB".to_string()] );
    }

    #[test]
    fn test_merge() {
        let mut self_genes = IndexedGenes::empty(Some(0));
        let mut other_genes = IndexedGenes::empty( Some(0) );

        self_genes.get_gene_id("geneA");

        other_genes.get_gene_id("geneX");
        other_genes.get_gene_id("geneY");

        let res = self_genes.merge(&other_genes);
        assert_eq!(res, vec![1,2]);

        assert_eq!( other_genes.get_all_gene_names(), vec![ "geneX".to_string(), "geneY".to_string()] );

        assert_eq!( self_genes.get_all_gene_names(), vec!["geneA".to_string(), "geneX".to_string(), "geneY".to_string()] );
    }
}

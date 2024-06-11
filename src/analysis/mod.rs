/// all analysis classes implement an analysis strategy.
/// They are not only data containers, but logics containers - if that makes sense.

pub mod minimal_sam;
pub mod fast_mapper;
pub mod gene_mapper;
pub mod genomic_mapper;
pub mod analysis_te;

pub use crate::analysis::minimal_sam::MinimalSam as MinimalSam;
pub use crate::analysis::fast_mapper::Analysis as Analysis;
pub use crate::analysis::gene_mapper::AnalysisGeneMapper as AnalysisGeneMapper;
pub use crate::analysis::genomic_mapper::AnalysisGenomicMapper as AnalysisGenomicMapper;
pub use crate::analysis::analysis_te::AnalysisTE as AnalysisTE;




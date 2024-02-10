use clap::Parser;

#[derive(Parser, Clone)]
#[clap(version = "1.2.5", author = "Stefan L. <stefan.lang@med.lu.se>")]
pub struct OptsTE {
    /// the input R1 reads file
    #[clap(short, long)]
    pub reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    pub file: String,
    /// the specie of the library [mouse, human]
    #[clap(short, long)]
    pub specie: String,
    /// the outpath
    #[clap(short, long)]
    pub outpath: String,
    /// the index build from the TE transcripts
    #[clap(short, long)]
    pub tfs: Option<String>,
    /// the index containing the 'normal' transcripts
    #[clap(short, long)]
    pub expression: Option<String>,
    /// the minimum (UMI) reads per cell (sample + genes + expression combined)
    #[clap(short, long)]
    pub min_umi: usize,
    /// the version of beads you used v1, v2.96 or v2.384
    #[clap(short, long)]
    pub version: String,
    /// Optional: end the analysis after processing <max_reads> cell fastq entries 
    #[clap(default_value_t=usize::MAX, long)]
    pub max_reads: usize,
    /// minimal sequencing quality 
    #[clap(default_value_t=25.0, long)]
    pub min_quality: f32,
    /// minimal sequencing quality 
    #[clap(default_value_t=32, long)]
    pub gene_kmers: usize,
    /// how many threads to use to analyze this (default max available)
    #[clap(short, long)]
    pub num_threads: Option<usize>,
    /// this is a BD rhapsody (bd) or a 10x expression experiment( default 10x)? 
    #[clap(default_value="10x", long)]
    pub exp: String,
    /// this is a BD rhapsody (bd) or a 10x expression experiment( default 10x)? 
    #[clap(default_value_t=10_000, long)]
    pub chunk_size: usize,
}

impl OptsTE {
    // Create a public wrapper function for parse
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
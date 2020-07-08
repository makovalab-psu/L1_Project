### Pipeline for additional analysis of: 
1. Restriction enzyme MSPI and TaqI sites (sub-folder `Enzyme_site_analysis/`)
2. L1 target motifs: (sub-folder `L1_target_motif_analysis/`)
3. Genome-wide polyA/T distribution: (sub-folder `PolyAT_analysis/`)

#### Enzyme_site_analysis/
- The folder contains 2 subfolders: `Clustering/` and `Enzyme_site_distribution/`
- `Clustering/` contains 2 Rscripts and 1 RData for clustering analysis of restriction enzyme sites (MSPI and TaqI) with other genomic features in the study
- `Enzyme_site_distribution/` contains 1 Rscript and 1 RData for the genome-wide distance distribution analysis of MSPI and TaqI sites

#### L1_target_motif_analysis/
- The folder contains 2 subfolders: `Distance_L1_target/` and `IWT_L1_target/`
- `Distance_L1_target/` contains the scripts and RData for the analysis of distances between estimated de novo L1 insertions (both complete and strigently filtered set) and concensus target motifs 
- `IWT_L1_target/` contains the scripts and processed Rdata object for IWTomics test of L1 target motifs on the strigently filtered de novo L1 dataset

#### PolyAT_analysis/
- The folder contains 1 shell script and 6 interval files for the genome-wide polyA/T sequence analysis 


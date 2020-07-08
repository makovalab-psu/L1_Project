restriction enzyme sites, 
L1 target motifs
polyA/T distribution


### Pipeline for additional analysis of: 
1. Restriction enzyme MSPI and TaqI sites (sub-folder `Enzyme_site_analysis/`)
2. L1 target motifs: (sub-folder `L1_target_motif_analysis/`)
3. Genome-wide polyA/T distribution: (sub-folder `PolyAT_analysis/`)

#### Enzyme_site_analysis/
- The folder contains 2 subfolders: `Clustering/` and `Enzyme_site_distribution/`
- `Clustering/` contains 2 Rscripts and 1 RData for clustering analysis of restriction enzyme sites (MSPI and TaqI) with other genomic features in the study
- `Enzyme_site_distribution/` contains 1 Rscript and 1 RData for the genome-wide distance distribution analysis of MSPI and TaqI sites

#### L1_target_motif_analysis/
- The folder contains 2 Rscripts and 1 RData 
- `IWTomic_newL1_revision.R` contains the pipeline for IWTomics tests based on new L1s from recent studies (Flasch et al 2019)(Sultana et al 2019)
- `L1_complete_new_revision.RData` contains processed datasets and extracted features to be used in the IWTomics test

#### PolyAT_analysis/
- The folder contains 1 shell script and 6 interval files for the genome-wide polyA/T sequence analysis 


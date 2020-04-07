### IWTomics pipeline for: 
1. High-resolution features (sub-folder `high_resolution/`)
2. Low resolution features: (sub-folder `low_resolution/`)
3. Additional de novo L1s from recent studies: (sub-folder `newL1_recent_study/`)

#### high_resolution/
- The folder contains 3 Rscripts, 2 RData, and 2 subfolders 
- The 3 scripts should be run in the following order: `IWTomic_high_resolution.R`, `Median_pvalues.R`, and `Plot_features_reordered.R`
- IWTomics results and plots are generated in the 2 sub-folders(`median_scaled` and `plots_reordered`) accordingly

#### low_resolution/
- The folder contains 2 Rscripts and 1 RData 
- Rscript `IWT_L1_low_resolution.r` is the main pipeline, which loads `L1_low_resolution_complete.RData` and source the `IWTomicsData_low_resolution.r`

#### newL1_recent_study/
- The folder contains 1 Rscripts and 1 RData 
- `IWTomic_newL1_revision.R` contains the pipeline for IWTomics tests based on new L1s from recent studies (Flasch et al 2019)(Sultana et al 2019)
- `L1_complete_new_revision.RData` contains processed datasets and extracted features to be used in the IWTomics test

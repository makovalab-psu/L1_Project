## Multiple Functional Logistic Regression (mFLR) and variable selection

### The folder contains 2 R script files, 2 input txt files and 1 sub-folder:

- `func_log_regression.r`: main R script for mFLR. For each of the six comparisons, performs variable selection and fit the selected mFLR model.

- `funLASSO_model_modified_for_L1.r`: R functions needed to perform variable selection in a mFLR model, based on group LASSO as proposed in *Matsui Hidetoshi (2014) Variable and boundary selection for functional data via multiclass logistic regression modeling. Computational Statistics & Data Analysis, 78, 176-185*, and sourced by the main R script for mFLR. **These functions are based on the code provided by Hidetoshi Matsui**.

- `localization_table.txt`: table indicating which high-resolution features must be included in the mFLR model as *functional* or *scalar* predictors. Input to the main R script for mFLR.

- `localization_table_LowFeatures.txt`: table indicating which low-resolution features must be included in the mFLR model as *scalar* predictors. Input to the main R script for mFLR.

### In the subfolder `random1`, there are 2 RData and 6 pdf files:
1. `L1_transformed_random_1.RData`: transformed variables for random sample 1, post-transformation variables can be found in “result_mean” after loading the data.
(the raw files, plots and transformation pipeline can be find at folder:`L1/transformed_data` on google drive).

2. `L1_lowFeatures_transformed_random1.RData`: transformed variables for 7 low-resolution features random sample 1, post-transformation variables can be found in “regionsFeatures” after loading the data.

3. 6 pdf files are the BIC and classification error plots used to select the tuning parameter in the feature selections, with each comparison written separately: 
- `comp1_feature_selection_type1.pdf`
- `comp2_feature_selection_type1.pdf`
- `comp3_feature_selection_type1.pdf`
- `comp4_feature_selection_type1.pdf`
- `comp5_feature_selection_type1.pdf`
- `comp6_feature_selection_type1.pdf`

## Multiple Functional Logistic Regression (mFLR) and variable selection

### The folder contains 2 R script files, 2 input txt files and 1 sub-folder:

- `func_log_regression.r`: 

- `funLASSO_model_modified_for_L1.r`: R functions needed to perform variable selection in a multiple Functional Logistic Regression model, based on group LASSO as proposed in *Matsui Hidetoshi (2014) Variable and boundary selection for functional data via multiclass logistic regression modeling. Computational Statistics & Data Analysis, 78, 176-185*, and sourced by the main R script for mFLR. **These functions are based on the code provided by Hidetoshi Matsui**.

- `localization_table.txt`: table indicating which high-resolution features must be included in the mFLR model as *functional* or *scalar* predictors. Input to the main R script for mFLR.

- `localization_table_LowFeatures.txt`: table indicating which low-resolution features must be included in the mFLR model as *scalar* predictors. Input to the main R script for mFLR.

Sub-folder:
- `random1`:

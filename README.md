# Folder Organization

## scripts/
Inside this folder you find all the scripts used to generate the Target Variable, representing Chromatin Accessibility during Tumor Onset.
They are organized as follows:

- ./obsolete : scripts I do not use anymore, will probably be eliminated in the near future
- start_jupyter.sh starts a jupyter server on the workstation
- filter_fragments.sh filters the CellRanger output fragments.tsv file, retaining only fragments falling in a defined set of peaks.
- HELPER_FUNCTIONS.R contains all smaller helper functions used along all the scripts
- All files starting with H refer to the main helper function. The number after the H refers to the script they refer to
- All main processing steps are divided in different scripts and numbered according to the order to be followed

## ML/
Inside this folder you find all the Machine Learning Related jupyter notebooks and all the objects saved as an output. 
Each Machine Learning approach (Regression, Binary Classification, Multiclass Classification) is split among 3 different notebooks:

  -  \<approach\>_Pipeline: Contains the model training and saves all needed files in \<approach\>_Output/
  -  \<approach\>_Feature_Interpretation: Containing explorative analysis for model debugging and feature interpretation
  -  \<approach\>_Test_Set_Evaluation: Contain final evaluations on Test Set

    NOTE: All random seeds have been set in order to ensure a fully reproducible analysis


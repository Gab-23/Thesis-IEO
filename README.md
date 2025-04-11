# Folder Organization

## scripts
Inside this folder you find all the scripts used to generate the Target Variable, representing Chromatin Accessibility during Tumor Onset.
They are organized as follows:

- ./obsolete : scripts I do not use anymore, will probably be eliminated in the near future
- start_jupyter.sh starts a jupyter server on the workstation
- filter_fragments.sh filters the CellRanger output fragments.tsv file, retaining only fragments falling in a defined set of peaks.
- HELPER_FUNCTIONS.R contains all smaller helper functions used along all the scripts
- All files starting with H refer to the main helper function. The number after the H refers to the script they refer to
- All main processing steps are divided in different scripts and numbered according to the order to be followed
- 


## ML

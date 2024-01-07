### Scripts to reproduce BANKSY analyses. 

Note: This is not the main BANKSY repo. That can be found here: [R version](https://github.com/prabhakarlab/Banksy), [Python version](https://github.com/prabhakarlab/Banksy_py/tree/main).

The scripts in this repo run with BANKSY version 0.1.5. This is a legacy version that uses a 
custom BANKSY object (i.e., not the SingleCellExperiment class in the latest version).
You can install it using the line: 
remotes::install_github("prabhakarlab/Banksy@main") # to install 0.1.5.

There are the following major folders in this repository: 
- fig3-crc
- fig3-hypothalamus
- fig4-hippocampus
- fig5-dlpfc
- suppfig-codex
- suppfig-cosmx
- artificial-data

The scripts to regenerate the corresponding results are in the src/ folders in 
these top level directories. The data must be downloaded and placed in the data/ 
folders, and the instructions to download the data can be found in either the 
data/README file, or in the beginning of the respective script in the src directory itself. 

The SlideSeq and STARmap data analyses can be found in the Jupyter notebooks at the [banksy_manuscript branch](https://github.com/prabhakarlab/Banksy_py/tree/Banksy_manuscript) of the Python repo. 

This folder contains config files (.cfg) for all original traits, SGCT (ID_SH.cfg) and UGCTs from the set of anthropometric traits. Herein is also a script running DEPICT analysis using config files for all traits (.sh). The analysis was performed with significance threshold 5e-08 using GWAS data without genomic control.

To start DEPICT it is necessary to install required packages and pandas 0.17.0, to do it, please, use command:

pip install pandas==0.17.0

inside conda environment for python 2.7

## 01_antro_heatmap.R
This script creates a heatmap visualizing scaled matrix of overlapping gene sets between anthropometric traits.


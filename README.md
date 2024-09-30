# TARP_motu_comparison Repository

##Re-connecting ecosystems: Integrating coral reefs into the monitoring of island restoration##
###Benkwitt et al. In Review###

This repository holds the data and analyses for our comparative study (Benkwitt et al. In Review) among motu on Tetiaroa, French Polynesia that dispaly varying rat invasion histories.
All analyses were done on data collected from Tetiaroa in November 2021 

###Repository Organisation###

All data can be found in this repository EXCEPT for raw sequence data for microbial communities. 
Coral and seawater microbial 16S rRNA amplicon sequences can be accessed through the NCBI Sequence Read Archive (SRA) under the accession # PRJNA11146751. 

Data streams included in analysis here: 
>1) Algae isotopes
2) Coral and water microbes
3) Benthic cover
4) Fish surveys
5) Bird surveys
6) Algae surveys

Each directory includes the sub-directories: input (or raw_data), output and procedure.

>input = any necessary input data files or raw data files (some directories may not include an input file)
output = any output files, including data tables and figures, from scripts
procedure = any code/script file (most included here are either R, R Markdown or bash scripts)


**For the integrated analyses (e.g., metric vs. seabirds or algal N15) used in the associated manuscript, please see the cross_dataset_comparisons directory**

*Important Note:*
Deprecated scripts and files (i.e., those that include old or inaccurate R code and associated output) are included in either 'deprecated' or 'old_results' directories throughout this repository. 
These are kept only for reference and should NOT be used for replication of results.


**For questions, please contact C. Benkwitt or H. Epstein.**
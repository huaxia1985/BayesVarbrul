# BayesVarbrul

A unified multidimensional analysis of language change in a speaker community

Citation: Hua, X. 2021. BayesVarbrul: A unified multidimensional analysis of language change in a speaker community. BioRxiv.

BayesVarbrul.R implements a Bayesian hierarchical model that adapts the concepts in genome-environment association studies to study language evolution in a speaker community.
instruction.R gives step-by-step instructions to use the method, using the Gurindji Kriol case study in the paper. It also gives the code to generate each figure in the paper.

rawdata.csv and var.csv are the data files for the Gurindji Kriol case study. rawdata.csv gives the social factor and the langauge usage pattern of each speaker. Variants of the same variable i are all named Vi. The second row in the file gives the variant type for each corresponding variant. var.csv gives the name and type of each variable in rawdata.csv.

For corpus data, the data files are in the same format, except that the data is not 0/1, but the number of times a speaker uses a particular variant in the corpus.

To use the code:
1) download all the files to a local folder
2) install the required libraries (MCMCpack, mvtnorm, and stringr).
3) to use the code for user's own dataset, prepare the data files in the same format as rawdata.csv and var.csv
4) follow instruction.R to run each step in R

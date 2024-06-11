# Ccap_Phylogeography_and_Climate_adaptation
This repo contains R code to run analyses performed in the submitted manuscript titled "Extensive admixture on a global scale and local climatic adaptation of the Mediterranean fruit fly (Diptera, Tephritidae: Ceratitis capitata) revealed by whole genome sequencing" by Deschepper et al. 2024

findgraph.R: A script to start admixture graph searching based on f2 statistics. This script can be called with findgraph.slurm, in which we specify the number of admixture events in the graph M and the number of graph searching iterations I. The output are files containing the best model for each model searching iteration. With testmodels.R, you can combine the outputted models into a dataframe and test whether a specific model has significantly better fit than another one.

Deschepper_etal_2024_Phylogeography_ClimateAdaptation_Ccap.R: This R script contains the code to run analyses concerning gene-climate association using PCadapt, RDA and LFMM. The script is annotated and organisized in different sections, including a section on how to prepore the genotypes dataframe.

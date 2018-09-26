# Drosophila pseudoobscura SR Analyses
------------------------------------------
### This repository contains data and Python and R scripts used in "Extensive recombination suppression and chromosome-wide differentiation of a segregation distorter in Drosophila" 
Fuller Z.L., Koury S.A., Leonard C.J., Young R.E , Richards S., Schaeffer S.W. & Phadnis N. (2018)

Descriptions for the scripts contained in the repository:
 - dxy_freq.py - Python script to calculate Dxy and Fst in genomic windows and estimate bootstrapped confidence intervals. 
 - SR_analyses.R - R scripts and commands used to estimate all statistics and generate plots that appear in the manuscript for all analyses, but the tests for differential gene expression
 - SR_gene_expression.R - R script used to perform the test for differential gene expression and generate Figure 5 in the mansucript
 
Data dervided from the raw Fst files and that is used as input into all of the above scripts are included in the data/ directory. These files contain genotype counts for each strain at each variable site, and are seperated by chromosome.

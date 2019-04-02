# LG_to_Chromosome_SyntenyPlots
This repository contains a custom R script to plot Circos plots using Circlize to identify synteny between linkage groups and chromosomes in a reference genome

If you use/adapt any of this script please cite this repository: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2622015.svg)](https://doi.org/10.5281/zenodo.2622015)

The script is: LinkageMap_Chromosome_Syntey_circosplot.R

Input files are:
1. FalconChromosomeLengths.txt
- the lengths of fasta headings and 40 chromosome sequences

2. LM_stats.csv
- summarised from De-Kayne and Feulner 2018 including linkage group names and lengths (cM)

3. SA_linkagemap_FalconPhaseMapped_Filtered.csv
- a modified sam file from mapping linkage map markers to the reference genome

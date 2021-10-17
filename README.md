#This sum-based method for GWAS signal enrichment analyses 

1.Introduction

The sum-based method uses signals of all markers within a pre-defined candidate feature. Briefly, we calculated the following summary statistics for candidate region: 
<img width="322" alt="image" src="https://user-images.githubusercontent.com/36602011/137618373-c5fc8cf2-7e6e-4a70-aea7-55f90805a6d5.png">

In which, T_sum is the summary statistics for a tested feature group. m_g is the number of SNPs located in candidate feature, and β is the estimate of marker effect obtained from GWAS summary statistics. Using formula 1, we calculated the T_sum for candidate feature. 

Usage:

perl GWAS_enrichment.pl -a [genome_region.bed] -b [GWAS_summaries.txt] -g [specific.regions] -e [genome_region_extention] -t [permutation_times]

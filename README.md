#This sum-based method for GWAS signal enrichment analyses 

1.Introduction

The sum-based method uses signals of all markers within a pre-defined candidate feature. Briefly, we calculated the following summary statistics for candidate regions: 

       <img width="160" alt="image" src="https://user-images.githubusercontent.com/36602011/137618373-c5fc8cf2-7e6e-4a70-aea7-55f90805a6d5.png">

In which, <img width="40" alt="image" src="https://user-images.githubusercontent.com/36602011/137618450-738015d4-7dce-4b08-98ec-86aeb4154e06.png">is the summary statistics for a tested feature group. <img width="30" alt="image" src="https://user-images.githubusercontent.com/36602011/137618470-8dba4886-6880-4adb-b97f-5d53b40b35f9.png">is the number of SNPs located in candidate feature, and β is the estimate of marker effect obtained from GWAS summary statistics. Using formula 1, we calculated the T_sum for candidate regions. 

Usage:

perl GWAS_enrichment.pl -a [genome_region.bed] -b [GWAS_summaries.txt] -g [specific.regions] -e [genome_region_extention] -t [permutation_times]

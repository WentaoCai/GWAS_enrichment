# This Sum-based method for GWAS Signal Enrichment analysis (SumGSE)

## 1. Introduction

SumGSE is a tool to integrate genomic information of biological mechanisms with GWAS summary statistics for complex traits. Almost all of the software here is command-line based.

The sum-based method uses signals of all markers within a pre-defined candidate feature. Briefly, we calculated the following summary statistics for candidate regions: 

   <img width="160" alt="image" src="https://user-images.githubusercontent.com/36602011/137618373-c5fc8cf2-7e6e-4a70-aea7-55f90805a6d5.png">

In which, <img width="40" alt="image" src="https://user-images.githubusercontent.com/36602011/137618450-738015d4-7dce-4b08-98ec-86aeb4154e06.png">is the summary statistics for a tested feature group. <img width="30" alt="image" src="https://user-images.githubusercontent.com/36602011/137618470-8dba4886-6880-4adb-b97f-5d53b40b35f9.png">is the number of SNPs located in candidate feature, and β is the estimate of marker effect obtained from GWAS summary statistics. Using this formula, we calculated the <img width="40" alt="image" src="https://user-images.githubusercontent.com/36602011/137618450-738015d4-7dce-4b08-98ec-86aeb4154e06.png"> for candidate regions. 

## 2. Getting Started

In order to download SumGSE, you should clone this repository via the commands

   ```shell
   git clone https://github.com/WentaoCai/GWAS_enrichment.git 
   cd GWAS_enrichment
   ```   
Once the above has completed, you can run:
   `SumGSE.pl -h`

Usage:

   `SumGSE.pl -a [genome_region.bed] -b [GWAS_summaries.txt] -g [specific.regions] -e [genome_region_extention] -t [permutation_times]`

## 3. Citation

If you use the software, please cite:

[Integrated Small RNA Sequencing, Transcriptome and GWAS Data Reveal microRNA Regulation in Response to Milk Protein Traits in Chinese Holstein Cattle. Frontiers in Genetics, 2021.](https://www.frontiersin.org/articles/10.3389/fgene.2021.726706/full "Citation Paper")


## Author

Wentao Cai, Institute of Animal Science of CAAS
Issues with SumGSE? Email: wtaocai@gmail.com

# This Sum-based method for GWAS Signal Enrichment analysis (SumGSE)

## 1. Introduction

SumGSE is a tool to integrate genomic information of biological mechanisms with GWAS summary statistics for complex traits. All of the software here is command-line based.

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



### Usage 1: The enrichment of GWAS signals for your chosed regions.

If you want to check if your chosed regions were more enriched with GWAS signals. you can used:

   `perl SumGSE.pl -i [genome_region.bed] -g [GWAS_summaries.txt]`
 
 The usage 1 is useful for the feature regions, such as lncRNAs, ChIP/ATAC peak et al.  
 
 Example: 
 
 `perl SumGSE.pl -i lncRNA.test.bed -g GWAS_statistics.txt -e 50 -n 1000`
 
 
   
### Usage 2: The enrichment of GWAS signals of your chosed regions limited in specific region.

If you want to check if your chosed regions were more enriched with GWAS signals in the specific regions. you can used:

   `perl SumGSE.pl -i [genome_region.bed] -b [GWAS_summaries.txt] -s [specific_regions.bed]`
   
The usage 2 may be useful to check the enrichment of differentailly expressed (genes/methylation/...）compare to that of all (genes/methylation...) in genome.

Example: 
 
 `perl SumGSE.pl -i gene.test.bed -g GWAS_statistics.txt -s swine.gene.bed -e 50 -n 1000`


### Options:

        -i    input file in bed format (Required). The input file should be genome feature regions(such as DEGs，miRNAs targets, Peaks from ChIP-seq, ATAC, et al.)  The first three columns should be chromosome, start, end. Example: 1   567821 573421  EEF1D

        -g    GWAS summary statistics (Required). The first two columns should be chromosome and position, the last column shoud be effect values,such as t value or beta value. Example: 1  123089 rs0011345  0.00045  -1.4625

        -e    extended range (KB) for genome feature regions (Optional). For example, the -e 100 means genome feature regions should also include their unstream/downstream 100kb region. Default -e is 0.

        -n    repeat n times for the permutation test (Optional). Default -n is 1000.

        -o    output file (Optional). Default the output file name is "SumGSE_permutation.out".

        -s    specific regions (Optional). If assuming -s, the permutation SNPs will be limited in these specific regions.

## 3. Citation

If you use the software, please cite:

[Integrated Small RNA Sequencing, Transcriptome and GWAS Data Reveal microRNA Regulation in Response to Milk Protein Traits in Chinese Holstein Cattle. Frontiers in Genetics, 2021.](https://www.frontiersin.org/articles/10.3389/fgene.2021.726706/full "Citation Paper")


## Author

Wentao Cai, Institute of Animal Science of CAAS

Issues with SumGSE? Email: wtaocai@gmail.com

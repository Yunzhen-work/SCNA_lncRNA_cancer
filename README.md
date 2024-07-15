# lncRNAs occur frequent somatic copy number alternation (SCNA)

## Project Overview

This project is to identify lung-adnocarcinoma-related SCNA lncRNAs and explore the expression characteristics and biological functions of the SCNA lncRNAs. Several R scripts are used for specific tasks. Check the following list of functions and select the one needed.

## Functions

- 450k_rawdata.r
    Process raw data of methylation merge it into matrix
  
- Diff_gene.r
    Process differential expression analysis using SAM package
  
- Dosage_sensitivity_score_calculate.r
    Calculate dosage sensitivity score for each lncRNA
  
- Expr_Normalized.r
    Normalization of TCGA expression data
  
- Filtered_ceRNA.r
    Filter lncRNA-miRNA-mRNA ceRNA relationship
  
- Hallmark_analysis.r
    LncRNA-related Hallmark analysis
  
- LncRNA_methylation.r
    Identify methylation in lncRNA
  
- LncRNA_nearby_gene.r
    Search for nearby genes of lncRNA and calculate correlation of each lncRNA and its neighboring genes
  
- LncRNA_promoter_indentified.r
    Identify promoter region of lncRNA
    
- Multiple_survival_analysis.r
    Multiple regression analysis
  
- Process_NA_methylation.r
    Process missing methylation values (Using mean value)
  
- TCGA_expr_process.r
    Process gene expression data downloaded from TCGA
  
- Univariate_survival_analysis.r
    Univariate cox analysis for lncRNA

## Usage

1. Clone the repository
git clone https://github.com/Yunzhen-work/SCNA_lncRNA_cancer.git

2. Enter R environment

3. Install the following dependencies
    impute matrixStats samr ggplot2

4. Configure the working directory in the script and run the script as needed

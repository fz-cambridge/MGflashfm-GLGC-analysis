# MGflashfm-GLGC-analysis
The MGflashfm-GLGC analysis scripts and instructions

## Replication

This repository holds the scripts and instructions used to generate the results for the MGflashfm fine-mapping paper (cited below). 

For the repository containing the source code of MGflashfm R package, see: https://github.com/...

Paper citation:

> *Leveraging information between multiple population groups and traits improves fine-mapping resolution* <br />
> F Zhou, O Soremekun, T Chikowore, S Fatumo, I Barroso, AP Morris, JL Asimit <br />
> bioRxiv ...; doi: ...

For more details, see the folder (GLGC_data_scripts/) and the paper (Methods and Supplementary Materials)

## Background
GLGC performed a multi-group genome-wide meta-analysis of lipid levels in 1.65 million people[^1]. Their multi-group meta-analysis included five genetically similar groups that are labelled by continent. For consistency, we have kept the same group labels as those published: admixed African or African (AFR, N=99,432, 6.0% of the sample); East Asian (EAS, N=146,492, 8.9%); European (EUR, N=1,320,016, 79.8%); Hispanic (HIS, N=48,057, 2.9%); and South Asian (SAS, N=40,963, 2.5%). We consider four of their five blood lipids traits: low-density lipoprotein cholesterol (LDL), high-density lipoprotein cholesterol (HDL), triglycerides (TG), and total cholesterol (TC); non-high-density lipoprotein cholesterol (nonHDL-C) is excluded due to its higher number of missing variants in any given region, compared to the other four traits. 

We used MGflashfm and MGfm to fine-map signals among HDL, LDL, TG, and TC in the five groups, and provide the following scripts and data:

- *GLGC_analysis.R:* the main analysis script that prepares the data for fine-mapping analysis and runs the analyses
- *QCfunctions.R:* a dependency in the GLGC_analysis.R script, here we provide functions that we used to perform quality control on the GWAS summary data from each population group
- *GLGC_input_50regions_chr_pos_index_100k.csv:* Details of fine-mapped regions

## Reference
[^1]: Graham, S. E. et al. (2021) *The power of genetic diversity in genome-wide association studies of lipids.*. Nature 600, 675–679. [https://www.nature.com/articles/s41467-020-14791-2](https://www.nature.com/articles/s41586-021-04064-3).

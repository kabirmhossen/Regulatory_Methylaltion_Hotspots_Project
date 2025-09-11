# Identification of Regulatory Methylation Hotspots in HIV Infection through Integrated Omics Approach and Machine Learning Framework

## Overview
This study explores the interplay between transcriptomic and epigenetic alterations in chronic HIV infection, aiming to uncover regulatory methylation hotspots (RMHs) that influence immune dysfunction. Chronic HIV infection is known to disrupt immune homeostasis, and this project hypothesizes that such disruption is partly due to coordinated transcriptomic and DNA methylation changes.
Using an integrated omics approach combining microarray expression and methylation data, the study identified 85 genes under coordinated regulation. Machine learning was applied to this integrated dataset to prioritize seven key RMHs: **_LTBP4, STAT1, SULF2, SLC2A14, SNX2, DUOX1, and TUBGCP2_**. These genes were found to be functionally linked to pathways involving TGF-β signaling, immune suppression, and redox imbalance.

## Poster Presented in the Final Year Project Demonstration
![Poster Presentation](figures/High_Res_Poster-01.jpg)

## Dependencies
# R Packages (Illumina 450k/EPIC arrays)
- minfi v1.50.0
- ChAMP v2.30.0
- DMRcate v2.16.0
- limma v3.60.6
- missMethyl v1.36.0 (gometh for GO/KEGG from CpGs)
- sva v3.50.0 (ComBat batch correction)
- IlluminaHumanMethylation450kanno.ilmn12.hg19 v0.6.1
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19 v0.6.0
- VennDiagram v1.7.3
- ComplexHeatmap v2.20.0
- ggplot2 v3.5.1
- corrplot v0.92
- GenomicRanges v1.54.1
- annotatr v1.24.0
- BiocParallel v1.38.0

# Python Libraries (optional / workflow & downstream)
- pandas v2.1.4
- numpy v1.26.4
- scikit-learn v1.4.2
- scipy v1.11.4
- xgboost v1.7.6
- lightgbm v4.3.0
- methylprep v1.10.0 (EPIC/450k IDAT processing in Python)
- methylcheck v1.6.0
- matplotlib v3.8.2
- seaborn v0.13.2

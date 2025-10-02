# MCMV_transcriptome

**Supports the publication:**  
**_Global transcriptional reprogramming by cytomegalovirus infection suppresses antigen presentation while enhancing migration machinery in murine dendritic cells._**

This repository provides processed transcriptomic data and network files generated from the RNA-seq analysis of murine dendritic cells infected with wild-type MCMV or a signaling-deficient M33 mutant (M33NQY). The data support the figures and supplementary tables described in the manuscript.

## 📁 Repository Contents

### Analysis Script
- `analysis.R`  added 3/10/25
  R script used to generate the processed expression data and differential gene expression tables included here (`Normlist_mouse_genes.csv`, `Normlist_MCMVgenes.csv`, `S2_Mockvs_IC2.csv`, `S3_Mockvs_NQY.csv`, and `S4_IC2vsNQY_DEGs.csv`).  
  Gene set enrichment (Tables S5–S7) and the Cytoscape network file (`Interactome_Map.cys`) were created using the methods described in the manuscript and are not produced by this script.

### Differential Gene Expression Tables
- `S2_Mockvs_IC2.csv`  
  **Table S2**: Differentially expressed genes in wild-type (IC2) MCMV-infected cells vs. mock-infected controls.

- `S3_Mockvs_NQY.csv`  
  **Table S3**: Differentially expressed genes in M33NQY-infected cells vs. mock-infected controls.

- `S4_IC2vsNQY_DEGs.csv`  
  **Table S4**: Differentially expressed genes in wild-type (IC2) vs. M33NQY-infected cells.

### Gene Set Enrichment Analysis (easyGSEA) Output
- `S5_Mockvs_IC2_KEGG-RA-BP-CC-MF_IC2_vs_Mock.csv`  
  **Table S5**: easyGSEA results (KEGG, Reactome, GO:BP/CC/MF) comparing IC2 vs. mock.

- `S6_Mockvs_NQY_KEGG-RA-BP-CC-MF_NQY_vs_Mock.csv`  
  **Table S6**: easyGSEA results comparing M33NQY vs. mock.

- `S7_IC2vsNQY_DEGs_KEGG-RA-BP-CC-MF_IC2_verses_NQY.csv`  
  **Table S7**: easyGSEA results comparing IC2 vs. M33NQY.

### Expression Data
- `Normlist_mouse_genes.csv`  
  TMM-normalized count data for mouse genes (mm10), generated via edgeR and used for expression analysis in Figures 1–3.

- `Normlist_MCMVgenes.csv`  
  TMM-normalized count data for MCMV genes generated via edgeR and used for expression analysis in Figure 4.

### Network Files
- `Interactome_Map.cys`  
  Cytoscape network file for Figure 3C. This includes merged antigen presentation and leukocyte migration networks, annotated with log₂ fold change and betweenness centrality.

## Computational Environment
Analyses were run in:
- R 4.4.2 (2024-10-31)
- RStudio 2024.12.0+467
- Major packages: tidyverse 2.0.0 (ggplot2 3.5.2, dplyr 1.1.4, tidyr 1.3.1), lubridate 1.9.4, pheatmap 1.0.13  
Full session information including all package versions and dependencies is provided in `sessionInfo.txt`.

---

## 🧬 Citation  
If using this dataset, please cite:  
Parry, R.H., McMillan, C.L.D., Bruce, K., Stevenson, P.G., Farrell, H.E. (2025). *Global transcriptional reprogramming by cytomegalovirus infection suppresses antigen presentation while enhancing migration machinery in murine dendritic cells.*

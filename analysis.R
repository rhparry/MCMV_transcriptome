############################################################
# analysis.R
# Rhys Parry r.parry@uq.edu.au
# Global transcriptional reprogramming by MCMV in DCs
# edgeR analysis to produce README-listed CSV outputs
#
# Inputs:
#   - Mouse_CMV_counts_fixed.csv   (host gene counts; 1st col = gene symbol)
#   - MCMV_counts.csv              (viral gene counts; 1st col = gene id)
#
# Outputs (exact filenames to match README):
#   - Normlist_mouse_genes.csv
#   - Normlist_MCMVgenes.csv
#   - S2_Mockvs_IC2.csv
#   - S3_Mockvs_NQY.csv
#   - S4_IC2vsNQY_DEGs.csv
#
# Notes:
# - Sections marked "EXPLORATORY – NOT USED IN FINAL PAPER" were used for QC/
#   visualisation and do not drive reported results.
# - GSEA (S5–S7) and Cytoscape session (Interactome_Map.cys) are generated
#   outside R as described in the manuscript Methods.
############################################################

# Reproducibility & defaults -----------------------------------------------------
set.seed(42)
options(stringsAsFactors = FALSE)

# Packages ----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  # library(EnhancedVolcano)
})

# I/O paths ---------------------------------------------------------------------
host_counts_file  <- "Mouse_CMV_counts_fixed.csv"
virus_counts_file <- "MCMV_counts.csv"

stopifnot(file.exists(host_counts_file),
          file.exists(virus_counts_file))

# Load counts (assumes first column is gene identifier) -------------------------
rawcounts      <- read.csv(host_counts_file,  check.names = FALSE)
rawcounts_MCMV <- read.csv(virus_counts_file, check.names = FALSE)

# Sample metadata ---------------------------------------------------------------
# IMPORTANT: Order these to match the column order of your CSVs.
# Host (9 samples): Mock x3, IC2 x3, NQY x3  (as used in the manuscript text)
group <- factor(c("Mock","Mock","Mock",
                  "IC2","IC2","IC2",
                  "NQY","NQY","NQY"),
                levels = c("IC2","Mock","NQY"))  # explicit order for contrasts

# Virus (6 samples): IC2 x3, NQY x3 (no Mock)
group_virus <- factor(c("IC2","IC2","IC2","NQY","NQY","NQY"),
                      levels = c("IC2","NQY"))

# Basic input sanity checks -----------------------------------------------------
stopifnot(ncol(rawcounts)      == 1 + length(group))       # 1 gene col + 9 samples
stopifnot(ncol(rawcounts_MCMV) == 1 + length(group_virus)) # 1 gene col + 6 samples

# DGEList objects ---------------------------------------------------------------
CMV_RNA <- DGEList(counts = as.matrix(rawcounts[, 2:ncol(rawcounts)]),
                   group  = group,
                   genes  = data.frame(gene = rawcounts[, 1], check.names = FALSE))

CMV_virusgenome <- DGEList(counts = as.matrix(rawcounts_MCMV[, 2:ncol(rawcounts_MCMV)]),
                           group  = group_virus,
                           genes  = data.frame(gene = rawcounts_MCMV[, 1], check.names = FALSE))

# Filtering (matches manuscript thresholds) -------------------------------------
keep_host  <- filterByExpr(CMV_RNA,         min.total.count = 10)
keep_virus <- filterByExpr(CMV_virusgenome, min.count       = 2)

CMV_RNA         <- CMV_RNA[keep_host, ,        keep.lib.sizes = FALSE]
CMV_virusgenome <- CMV_virusgenome[keep_virus, , keep.lib.sizes = FALSE]

# Normalisation (TMM) -----------------------------------------------------------
CMV_RNA         <- calcNormFactors(CMV_RNA)
CMV_virusgenome <- calcNormFactors(CMV_virusgenome)

# EXPLORATORY – NOT USED IN FINAL PAPER: quick library size plot ----------------
# par(mar = c(10, 6, 2, 2))
# col_group <- c("#4e5f95","#4e5f95","#4e5f95", "#5ac0c2","#5ac0c2","#5ac0c2", "#9aeeba","#9aeeba","#9aeeba")
# barplot(CMV_RNA$samples$lib.size * 1e-6,
#         names.arg = c("Uninfected_1","Uninfected_2","Uninfected_3",
#                       "IC2_1","IC2_2","IC2_3","NQY_1","NQY_2","NQY_3"),
#         ylab = "Library size (millions)", col = col_group, las = 2)

# Design matrices ---------------------------------------------------------------
# Host (3-level): columns = c("IC2","Mock","NQY")
design_host <- model.matrix(~0 + group)
colnames(design_host) <- levels(group)

# Virus (2-level): columns = c("IC2","NQY")
design_virus <- model.matrix(~0 + group_virus)
colnames(design_virus) <- levels(group_virus)

# Dispersion & model fitting ----------------------------------------------------
CMV_RNA         <- estimateDisp(CMV_RNA,         design_host,  robust = TRUE)
fit_host        <- glmQLFit(CMV_RNA,             design_host,  robust = TRUE)

CMV_virusgenome <- estimateDisp(CMV_virusgenome, design_virus, robust = TRUE)
fit_virus       <- glmQLFit(CMV_virusgenome,     design_virus, robust = TRUE)

# Normalised CPM outputs (README filenames) -------------------------------------
norm_host  <- cpm(CMV_RNA, normalized.lib.sizes = TRUE)
rownames(norm_host) <- CMV_RNA$genes$gene
write.csv(norm_host, file = "Normlist_mouse_genes.csv", row.names = TRUE)

norm_virus <- cpm(CMV_virusgenome, normalized.lib.sizes = TRUE)
rownames(norm_virus) <- CMV_virusgenome$genes$gene
write.csv(norm_virus, file = "Normlist_MCMVgenes.csv", row.names = TRUE)

# Differential expression (README filenames) ------------------------------------
# Column order is c("IC2","Mock","NQY"), so contrasts are numeric triplets in that order.
# S2: IC2 vs Mock
S2_Mockvs_IC2 <- topTags(
  glmQLFTest(fit_host, contrast = c(1, -1, 0)),
  n = Inf, sort.by = "PValue", adjust.method = "BH"
)
write.csv(S2_Mockvs_IC2, file = "S2_Mockvs_IC2.csv", row.names = FALSE)

# S3: NQY vs Mock
S3_Mockvs_NQY <- topTags(
  glmQLFTest(fit_host, contrast = c(0, -1, 1)),
  n = Inf, sort.by = "PValue", adjust.method = "BH"
)
write.csv(S3_Mockvs_NQY, file = "S3_Mockvs_NQY.csv", row.names = FALSE)

# S4: IC2 vs NQY
S4_IC2vsNQY <- topTags(
  glmQLFTest(fit_host, contrast = c(1, 0, -1)),
  n = Inf, sort.by = "PValue", adjust.method = "BH"
)
write.csv(S4_IC2vsNQY, file = "S4_IC2vsNQY_DEGs.csv", row.names = FALSE)


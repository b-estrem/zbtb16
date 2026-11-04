library(tidyverse)
library(Rsamtools)
library(MEDIPS)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9.masked)

out_path <- "~/szabo/brandon_estrem/03.zbtb16_mouse/analysis/output/"

# Constant Variables
genome="BSgenome.Mmusculus.UCSC.mm9.masked"
uniq=1e-3 
extend=300 # Assume single end reads are from 300 bp fragments
shift=0
ws=500 # Window size to divide genome - may need to change
chr.select="chr19"

##############################################
##### Load MeDIP samples
RS1 <- paste0(out_path, "bam/mCap_RS_r1.bam")
RS2 <- paste0(out_path, "bam/mCap_RS_r2.bam")
KIT1 <- paste0(out_path, "bam/mCap_KIT_r1.bam")
KIT2 <- paste0(out_path, "bam/mCap_KIT_r2.bam")
THY1 <- paste0(out_path, "bam/mCap_THY_r1.bam")
THY2 <- paste0(out_path, "bam/mCap_THY_r2.bam")

##############################################
##### Saturation analysis - QC (all look good)
# files <- list(
#   RS1  = RS1,
#   RS2  = RS2,
#   KIT1 = KIT1,
#   KIT2 = KIT2,
#   THY1 = THY1,
#   THY2 = THY2
# )
# for (name in names(files)) {
#     file <- files[[name]]
#     sr = MEDIPS.saturation(file = file, BSgenome = genome, uniq = uniq, extend = extend, shift = shift, window_size = ws, chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)
#     pdf(paste0("../plots/medips/", name, "_saturation.pdf"), width = 7, height = 7)
#     MEDIPS.plotSaturation(sr)
#     dev.off()
# }

##############################################
##### Create MEDIPS sets
RS_MeDIP <- MEDIPS.createSet(file = RS1, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)
RS_MeDIP <- c(RS_MeDIP, MEDIPS.createSet(file = RS2, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select))

KIT_MeDIP <- MEDIPS.createSet(file = KIT1, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)
KIT_MeDIP <- c(KIT_MeDIP, MEDIPS.createSet(file = KIT2, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select))

THY_MeDIP <- MEDIPS.createSet(file = THY1, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)
THY_MeDIP <- c(THY_MeDIP, MEDIPS.createSet(file = THY2, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select))

pdf("../plots/medips/correlation.pdf", width = 7, height = 7)
correlation <- MEDIPS.correlation(MSets=c(RS_MeDIP, KIT_MeDIP, THY_MeDIP), plot = T, method="pearson")
dev.off()

# ## CpG normalization vector
# CS <- MEDIPS.couplingVector(pattern = "CG", refObj = THY_MeDIP[[1]])

# run_medips <- function(name, df1, df2, cs) {
#   ## Do the comparison
#     diff_cov <- MEDIPS.meth(MSet1 = df1, MSet2 = df2, CSet = cs, p.adj = "bonferroni", diff.method = "edgeR", MeDIP = T, minRowSum = 10)

#     diff_cov.s <- MEDIPS.selectSig(results = diff_cov, p.value = 0.05, adj = T, ratio = NULL, bg.counts = NULL, CNV = F)
#     ## Subset significant regions into gain/loss
#     diff_cov.s.gain <- diff_cov.s[which(diff_cov.s[, grep("logFC", colnames(diff_cov.s))] > 0), ]
#     diff_cov.s.loss <- diff_cov.s[which(diff_cov.s[, grep("logFC", colnames(diff_cov.s))] < 0), ]
#     ## Merge neigbhoring signficant regions
#     diff_cov.s.gain.m <- MEDIPS.mergeFrames(frames = diff_cov.s.gain, distance = 1)
#     diff_cov.s.loss.m <- MEDIPS.mergeFrames(frames = diff_cov.s.loss, distance = 1)

#     ## Get data for ROIs
#     columns <- names(diff_cov)[grep("counts", names(diff_cov))]
#     gain.s <- MEDIPS.selectROIs(results = diff_cov, rois = diff_cov.s.gain.m, columns = columns, summarize = "median")
#     loss.s <- MEDIPS.selectROIs(results = diff_cov, rois = diff_cov.s.loss.m, columns = columns, summarize = "median")
#     # gq 
#     # zbtb16
#     # shuffled?

#     return(list(gain = gain.s, loss = loss.s))
# }

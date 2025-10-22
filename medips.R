library(tidyverse)
library(Rsamtools)
library(MEDIPS)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10.masked)

out_path <- "~/szabo/brandon_estrem/02.zbtb16/analysis/output/"

# Constant Variables
genome="BSgenome.Mmusculus.UCSC.mm10.masked"
uniq=1e-3 
extend=300 # Assume singl end reads are from 300 bp fragments
shift=0
ws=500 # Window size to divide genome - may need to change

RS1 <- BamFile(paste0(out_path, "bam/mCap_RS_r1.bam"))
RS2 <- BamFile(paste0(out_path, "bam/mCap_RS_r2.bam"))

RS_MeDIP <- MEDIPS.createSet(file = RS1, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws)
RS_MeDIP <- c(RS_MeDIP, MEDIPS.createSet(file = RS2, BSgenome = genome, extend = extend, shift = shift, uniq = uniq, window_size = ws))
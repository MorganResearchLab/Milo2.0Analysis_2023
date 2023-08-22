#! /usr/bin/env Rscript

library(SingleCellExperiment)

big.sce <- readRDS("Newcastle/Cambridge_Sanger_Newcastle_Tcell_reducedMNN.RDS")

tcell.annot <- read.table("Subclusters/Tcell_annotations_ext.tsv",
	                  sep="\t", header=TRUE, stringsAsFactors=FALSE)

all.meta <- read.table("Newcastle/COVID19_scMeta-data.tsv",
	               sep="\t", header=TRUE, stringsAsFactors=FALSE)

meta.merge <- merge(tcell.annot, all.meta, by='CellID')
rownames(meta.merge) <- meta.merge$CellID

# UMAPs
tcell.umap <- read.table("Newcastle/Cambridge_Sanger_Newcastle_Tcell_UMAPs.tsv",
	                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(tcell.umap) <- tcell.umap$CellID

gender.umap <- read.table("Newcastle/Cambridge_Sanger_Newcastle_Tcell-Gender_UMAPs.tsv",
	                  sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(gender.umap) <- gender.umap$CellID


# add annotations and other metadata
big.sce <- big.sce[, intersect(colnames(big.sce), rownames(meta.merge))]

colData(big.sce)$Annotation <- meta.merge[colnames(big.sce), ]$Sub.Annotation
colData(big.sce)$D0_status_summary <- meta.merge[colnames(big.sce), ]$D0_status_summary
colData(big.sce)$Age <- meta.merge[colnames(big.sce), ]$Age
colData(big.sce)$Sex <- meta.merge[colnames(big.sce), ]$Sex
colData(big.sce)$sample_id <- meta.merge[colnames(big.sce), ]$sample_id
colData(big.sce)$patient_id <- meta.merge[colnames(big.sce), ]$patient_id
colData(big.sce)$Collection_Day <- meta.merge[colnames(big.sce), ]$Collection_Day
colData(big.sce)$D0_status <- meta.merge[colnames(big.sce), ]$D0_status
colData(big.sce)$Days_from_onset <- meta.merge[colnames(big.sce), ]$Days_from_onset
colData(big.sce)$time_after_LPS <- meta.merge[colnames(big.sce), ]$time_after_LPS
colData(big.sce)$Smoker <- meta.merge[colnames(big.sce), ]$Smoker

### add UMAP co-ordinates
reducedDim(big.sce, "Regressed.UMAP") <- gender.umap[colnames(big.sce), c(1, 2)]
reducedDim(big.sce, "Tcell.UMAP") <- tcell.umap[colnames(big.sce), c(1, 2)]


saveRDS(big.sce, "Subclusters/Tcell_SCE.RDS")


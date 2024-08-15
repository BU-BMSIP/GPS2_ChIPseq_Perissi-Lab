if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ChIPseeker", "TxDb.Mmusculus.UCSC.mm39.knownGene", "org.Mm.eg.db"))

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)

peaks <- readPeakFile("/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/Treatment/ATF4/results/filtered/ATF4_WT_bf.narrowPeak")

txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")


output_dir <- "/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/WT/ATF3/results/annotation"
output_file <- file.path(output_dir, paste0(basename("ATF3"), "_annotated.txt"))
output_file
write.table(as.data.frame(peakAnno), file=output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Assuming you have your peakAnno object ready
plot <- plotAnnoPie(peakAnno)

# Specify the directory
dir <- "/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/WT/ATF3/results/plots"

# Create the directory if it does not exist
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

# Save the plot in the specified directory
filepath <- file.path(dir, "ATF3_annopiechart.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotAnnoPie(peakAnno)  # Call the plotting function again inside the device
dev.off()

# Bar plot of genomic annotation
filepath <- file.path(dir, "ATF3_annobar.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotAnnoBar(peakAnno)# Call the plotting function again inside the device
dev.off()

# Distribution of peaks relative to TSS
filepath <- file.path(dir, "ATF3_DistToTss.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotDistToTSS(peakAnno)# Call the plotting function again inside the device
dev.off()

filepath <- file.path(dir, "ATF3_covplot.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
covplot(peaks)# Call the plotting function again inside the device
dev.off()

filepath <- file.path(dir, "ATF3_peak_profile_heatmap.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
peak_Profile_Heatmap(peak = peaks,
                     upstream = 3000,
                     downstream = 3000,
                     by = "gene",
                     type = "start_site",
                     TxDb = txdb,
                     nbin = 800)
dev.off()


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peaks, windows=promoter)

filepath <- file.path(dir, "ATF3_tagHeatMap.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
tagHeatmap(tagMatrix)
dev.off()

filepath <- file.path(dir, "ATF3_plotAvgProf.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()
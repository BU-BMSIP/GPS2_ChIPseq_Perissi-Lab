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


output_dir <- "/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/Treatment/ATF4/results/annotation"
output_file <- file.path(output_dir, paste0(basename("ATF4_WT"), "_annotated.txt"))
output_file
write.table(as.data.frame(peakAnno), file=output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Assuming you have your peakAnno object ready
plot <- plotAnnoPie(peakAnno)

# Specify the directory
dir <- "/projectnb/perissilab/Jawahar/GPS2_ChIPseq_Perissi-Lab/Treatment/GPS2/results/plots"

# Create the directory if it does not exist
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

# Save the plot in the specified directory
filepath <- file.path(dir, "gps2_withFCCP_annopiechart.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotAnnoPie(peakAnno)  # Call the plotting function again inside the device
dev.off()

# Bar plot of genomic annotation
filepath <- file.path(dir, "gps2_withFCCP_annobar.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotAnnoBar(peakAnno)# Call the plotting function again inside the device
dev.off()

# Distribution of peaks relative to TSS
filepath <- file.path(dir, "gps2_withFCCP_DistToTss.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotDistToTSS(peakAnno)# Call the plotting function again inside the device
dev.off()

filepath <- file.path(dir, "gps2_withFCCP_covplot.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
covplot(peaks)# Call the plotting function again inside the device
dev.off()

filepath <- file.path(dir, "gps2_withFCCP_peak_profile_heatmap.png")
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

filepath <- file.path(dir, "gps2_withFCCP_tagHeatMap.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
tagHeatmap(tagMatrix)
dev.off()

filepath <- file.path(dir, "gps2_withFCCP_plotAvgProf.png")
png(filename = filepath, width = 8, height = 6, units = "in", res = 300)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()

# Additional plots
# Venn Diagram of overlapping peaks (example with two peak files)
if (length(bed_files) > 1) {
  peaks1 <- readPeakFile(bed_files[1])
  peaks2 <- readPeakFile(bed_files[2])
  vennPeaks(peaks1, peaks2, by="region")
  dev.copy(png, file.path(output_dir, paste0(basename(bed_file), "_vennDiagram.png")))
  dev.off()
}

# Correlation heatmap (example with dummy correlation matrix)
# Replace this with actual correlation data if available
cor_matrix <- matrix(runif(100), nrow=10, ncol=10)  # Dummy correlation matrix
rownames(cor_matrix) <- paste0("Peak_", 1:10)
colnames(cor_matrix) <- paste0("Sample_", 1:10)
corrplot(cor_matrix, method="color", type="full")
dev.copy(png, file.path(output_dir, paste0(basename(bed_file), "_correlationHeatmap.png")))
dev.off()

# Motif enrichment analysis
results <- motifbreakR(sites=peaks, pwmList=MotifDb, threshold=0.85)
motif_output_file <- file.path(output_dir, paste0(basename(bed_file), "_motifEnrichment.txt"))
write.table(results, file=motif_output_file, sep="\t", quote=FALSE, row.names=FALSE)
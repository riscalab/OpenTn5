20230822_OpenTn5.R
#Plotting sample correlation - counts in bins and counts in peaks 
setwd("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/MCF7+K562")
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(readr)
library(tidyverse)

#Counting accross tiled Bins in Genome
#make vector of bamfile paths
Libs = list.files("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/MCF7+K562")
Libs = Libs[grep(".bam$", Libs)]
bamfiles <- Libs
bamfiles

#[1] "K562Lot1pgTn5_5ul.atac.merged.bam"            "K562RiscaLabIllumina.atac.merged.bam"        
#[3] "MCF7Illumina_RiscaLab.atac.sorted.merged.bam" "MCF7Lot1pgTn5.atac.sorted.merged.bam"        
#[5] "MCF7Lot2pgTn5.atac.sorted.merged.bam"   


hg38 <- GRangesForBSGenome("hg38") # get GRanges for hg38 genome
#1000bp bins
hg38_1000bptiles <- tile(hg38, width=1000) # tile hg38 genome into 1000bp windows
hg38_1000bptiles <- hg38_1000bptiles[seqnames(hg38_1000bptiles) %in% paste0("chr", c(1:22, "X", "Y")), ] # keep only rows with chr1-22,chrX & chrY
hg38_1000bptiles <- unlist(hg38_1000bptiles) # unlist into GRanges from GRangesList
# count fragments under each 1000bp bin from bam files
Counts_1000bptiles <- chromVAR::getCounts(paste0(bamfiles), hg38_1000bptiles, paired=TRUE, by_rg=FALSE, format="bam")
Counts_1000bptiles <- assay(Counts_1000bptiles) # convert into counts matrix, matrix object
# add row names based on chromosome coordinates
row.names(Counts_1000bptiles) <- paste0(seqnames(hg38_1000bptiles), ":", start(hg38_1000bptiles), "-", end(hg38_1000bptiles))
save(Counts_1000bptiles, file = "Counts_1000bptiles_ChromVAR.RData")

#10000 bp bins (10kb)
hg38_10000bptiles <- tile(hg38, width=10000) # tile hg38 genome into 10000bp windows
hg38_10000bptiles <- hg38_10000bptiles[seqnames(hg38_10000bptiles) %in% paste0("chr", c(1:22, "X", "Y")), ] # keep only rows with chr1-22,chrX & chrY
hg38_10000bptiles <- unlist(hg38_10000bptiles) # unlist into GRanges from GRangesList
# count fragments under each 1000bp bin from bam files
Counts_10000bptiles <- chromVAR::getCounts(paste0(bamfiles), hg38_10000bptiles, paired=TRUE, by_rg=FALSE, format="bam")
Counts_10000bptiles <- assay(Counts_10000bptiles) # convert into counts matrix, matrix object
# add row names based on chromosome coordinates
row.names(Counts_10000bptiles) <- paste0(seqnames(hg38_10000bptiles), ":", start(hg38_10000bptiles), "-", end(hg38_10000bptiles))
save(Counts_10000bptiles, file = "Counts_10kbptiles_ChromVAR.RData")

###Counting reads under master peakset for k562+mcf7 samples together
peaks <-getPeaks("/lustre/fs4/risc_lab/scratch/landerson01/20220915_K562_OpenTn5/chromvar_mcf7+k562/opentn5master.merged.bed",sort_peaks=TRUE)
resize(peaks, width = 500, fix = "center")

opentn5fragment_counts <- getCounts(bamfiles, 
                                    peaks, 
                                    paired =  TRUE, 
                                    by_rg = FALSE, 
                                    format = "bam")

opentn5fragment_counts = addGCBias(opentn5fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
opentn5fragment_counts <- assay(opentn5fragment_counts)# convert into counts matrix, matrix object
save(opentn5fragment_counts, file = "opentn5_MCF7+K562_ChromVARfragmentcounts.RData")

#add metadata to matrix with coldata dataframe
#setup metadat for experiment for input into deseq
colData <- matrix(data = 0,nrow = 5) 
colData <- as.data.frame(colData)
rownames(colData) <- colnames(Counts_10000bptiles)

# add info about experiment into colData metadata table
colData$celltype <- c("K562","K562","MCF7","MCF7","MCF7")
colData$tn5 <-c("Lot1","Illumina","Illumina","Lot1","Lot2")
colData$sample <- c("K562 Lot1 pgTn5","K562 Illumina","MCF7 Illumina","MCF7 Lot1 pgTn5", "MCF7 Lot2 pg Tn5")


#import counts matrix into DEseq object
#10kb
tiles10kb_dds <- DESeq2::DESeqDataSetFromMatrix(countData=Counts_10000bptiles, design= celltype~tn5, colData= colData)
# rlog normalize DESeq object
tiles10kb_ddsrlog <- rlog(tiles10kb_dds, blind=TRUE)
# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(tiles10kb_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData$sample 
colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(tiles10kb_dds)[,c("celltype","tn5")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap

#1000bp
tiles1000bp_dds <- DESeq2::DESeqDataSetFromMatrix(countData=Counts_1000bptiles, design= celltype~tn5, colData= colData)
# rlog normalize DESeq object
tiles1000bp_ddsrlog <- rlog(tiles1000bp_dds, blind=TRUE)
# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(tiles1000bp_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData$sample 
colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(tiles1000bp_dds)[,c("celltype","tn5")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap


#Plot Spearman correlation of counts under peaks 
#import counts matrix into DEseq object, use colData from above
#peakcounts
peakcounts_dds <- DESeq2::DESeqDataSetFromMatrix(countData=opentn5fragment_counts, design= celltype~tn5, colData= colData)
peakcounts_ddsrlog <- rlog(peakcounts_dds, blind=TRUE) #rlog transform
# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(peakcounts_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData$sample 
colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(peakcounts_dds)[,c("celltype","tn5")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap

###Counting reads under CTCF sites 
peaks <-getPeaks("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/CTCF_bed_Encode_ENCFF278FNP/ENCFF278FNP.bed",sort_peaks=TRUE)
resize(peaks, width = 500, fix = "center")

CTCF_counts <- getCounts(bamfiles, 
                                    peaks, 
                                    paired =  TRUE, 
                                    by_rg = FALSE, 
                                    format = "bam")

CTCFfragment_counts <- assay(CTCF_counts)# convert into counts matrix, matrix object
save(CTCFfragment_counts, file = "opentn5_MCF7+K562_CTCF_ChromVARfragmentcounts.RData")

###September 5 2023 Plotting EColi Reads
setwd("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/AlignLogs")
samplenames <- c("K562Illumina_33","K562Illumina_34","K562Illumina_35","K562Lot1-1ul_36","K562Lot1-1ul_37","K562Lot1-2.5ul_38","K562Lot1-2.5ul_39","K562Lot1-2.5ul_40","K562Lot1-5ul_41","K562Lot1-5ul_42","K562Lot1-5ul_43","K562Lot1-7.5ul_44","K562Lot1-7.5ul_45","K562Lot1-7.5ul_46","MCF7Illumina-001","MCF7Illumina-002","MCF7Lot1-007","MCF7Lot1-008","MCF7Lot2-009","MCF7Lot2-010")
EColi_Reads_R1 <- list() # read 1
EColi_Reads_R2 <- list() # read 2
for(i in 1:20){
  EColi_Reads_R1[[i]] <- as.data.frame(read_delim(paste0(samplenames[[i]],"_R1_001.trim_screen.txt"),delim="\t",skip=1))
  EColi_Reads_R1[[i]]$samples <- samplenames[i] # include sample name 
  EColi_Reads_R1[[i]] <- EColi_Reads_R1[[i]][EColi_Reads_R1[[i]]$Genome %in% "Ecoli", ] # keep only row containing EColi reads info
  EColi_Reads_R2[[i]] <- as.data.frame(read_delim(paste0(samplenames[[i]],"_R2_001.trim_screen.txt"),delim="\t",skip=1))
  EColi_Reads_R2[[i]]$samples <- samplenames[i]
  EColi_Reads_R2[[i]] <- EColi_Reads_R2[[i]][EColi_Reads_R2[[i]]$Genome %in% "Ecoli", ] }

# convert list to dataframe
EColi_Reads_R1 <- do.call(rbind, EColi_Reads_R1)
EColi_Reads_R2 <- do.call(rbind, EColi_Reads_R2)

# remove "#" & "%" characters from column names
colnames(EColi_Reads_R1) <- gsub("#", "", colnames(EColi_Reads_R1))
colnames(EColi_Reads_R2) <- gsub("#", "", colnames(EColi_Reads_R2))
colnames(EColi_Reads_R1) <- gsub("%", "Percent", colnames(EColi_Reads_R1))
colnames(EColi_Reads_R2) <- gsub("%", "Percent", colnames(EColi_Reads_R2))

#make sample metadata info 
colData <- matrix(data = 0,nrow = 20) 
colData <- as.data.frame(colData)
rownames(colData) <- samplenames
colData$Enzyme <- c(rep("Illumina",3),rep("Lot1",11),rep("Illumina",2),rep("Lot1",2),rep("Lot2",2))
colData$CellType <- c(rep("K562",14),rep("MCF7",6))
colData$Replicate <- c("rep1,","rep2","rep3","rep1","rep2","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep1","rep2","rep1","rep2")
colData$Amount <- c(rep("2.5ul",3),rep("1ul",2),rep("2.5ul",3),rep("5ul",3),rep("7.5ul",3),rep("2.5ul",2),rep("5ul",4))
colData$samples <- samplenames

# merge EColi Reads dataframe with additional sample metaData information
EColi_Reads_R1 <- merge(EColi_Reads_R1, colData, by.y = "samples")
EColi_Reads_R2 <- merge(EColi_Reads_R2, colData, by.y = "samples")

# export raw data to csv for Jan
write.csv(EColi_Reads_R1, "ATAC_EColi_Reads_R1.csv")
write.csv(EColi_Reads_R2, "ATAC_EColi_Reads_R2.csv")

# plot mapped EColi reads comparison between Jan's pG-Tn5 vs' Epicypher's pAG-Tn5
library(ggplot2)
library(viridis)
ggplot(EColi_Reads_R1, aes(x=Enzyme, y=PercentOne_hit_one_genome, fill= Enzyme)) +
  geom_boxplot() +
  geom_point(aes(color = CellType)) +
  facet_wrap(~CellType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% EColi Reads") +
  xlab("")

ggplot(EColi_Reads_R1, aes(x=Enzyme, y=PercentOne_hit_one_genome, fill= Enzyme)) +
  geom_boxplot() +
  geom_point(aes(color = Amount)) +
  facet_wrap(~CellType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% EColi Reads") +
  xlab("")



#20240124
#plotting tech rep bams correlation
#Counting accross tiled Bins in Genome
#make vector of bamfile paths
setwd("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/techrepbams")
Libs = list.files("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/techrepbams")
Libs = Libs[grep(".bam$", Libs)]
bamfiles <- Libs
bamfiles

#[1] "K562_Illumina_33_S20_001.sorted.bam"        "K562_Illumina_34_S21_001.sorted.bam"       
#[3] "K562_Illumina_35_S22_001.sorted.bam"        "K562_Lot15ul_41_S39_001.sorted.bam"        
#[5] "K562_Lot15ul_42_S40_001.sorted.bam"         "K562_Lot15ul_43_S41_001.sorted.bam"        
#[7] "MCF7_001Cycling_Illumina.sorted.merged.bam" "MCF7_002Cycling_Illumina.sorted.merged.bam"
#[9] "MCF7_007Lot1Cycling.sorted.merged.bam"      "MCF7_008Lot1Cycling.sorted.merged.bam"     
#[11] "MCF7_009Lot2Cycling.sorted.merged.bam"      "MCF7_010Lot2Cycling.sorted.merged.bam"  

# count fragments under each 1000bp bin from bam files
Counts_1000bptiles <- chromVAR::getCounts(paste0(bamfiles), hg38_1000bptiles, paired=TRUE, by_rg=FALSE, format="bam")
Counts_1000bptiles <- assay(Counts_1000bptiles) # convert into counts matrix, matrix object
# add row names based on chromosome coordinates
row.names(Counts_1000bptiles) <- paste0(seqnames(hg38_1000bptiles), ":", start(hg38_1000bptiles), "-", end(hg38_1000bptiles))
save(Counts_1000bptiles, file = "Counts_1kbptiles__TECHREPS_ChromVAR.RData")

###Counting reads under master peakset for k562+mcf7 samples together
peaks <-getPeaks("/lustre/fs4/risc_lab/scratch/landerson01/20220915_K562_OpenTn5/chromvar_mcf7+k562/opentn5master.merged.bed",sort_peaks=TRUE)
resize(peaks, width = 500, fix = "center")

opentn5fragment_counts <- getCounts(bamfiles, 
                                    peaks, 
                                    paired =  TRUE, 
                                    by_rg = FALSE, 
                                    format = "bam")

opentn5fragment_counts = addGCBias(opentn5fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
opentn5fragment_counts <- assay(opentn5fragment_counts)# convert into counts matrix, matrix object
save(opentn5fragment_counts, file = "Counts_underPeaks_opentn5_techreps_ChromVARfragmentcounts.RData")

#add metadata to matrix with coldata dataframe
#setup metadat for experiment for input into deseq
colData <- matrix(data = 0,nrow = 12) 
colData <- as.data.frame(colData)
rownames(colData) <- colnames(opentn5fragment_counts)

# add info about experiment into colData metadata table
colData$celltype <- c(rep("K562",6),rep("MCF7",6))
colData$tn5 <-c(rep("Illumina",3),rep("pG-Tn5",3),rep("Illumina",2),rep("pG-Tn5",4))
colData$Lot <-c(rep("Illumina",3),rep("pG-Tn5 Lot 1",3),rep("Illumina",2),rep("pG-Tn5 Lot 1",2),rep("pg-Tn5 Lot 2",2))
#colData$sample <- c("K562 Lot1 pgTn5","K562 Illumina","MCF7 Illumina","MCF7 Lot1 pgTn5", "MCF7 Lot2 pg Tn5")

#import count object into DESeq2
#1000bp
tiles1000bp_dds <- DESeq2::DESeqDataSetFromMatrix(countData=Counts_1000bptiles, design= celltype~tn5, colData= colData)
# rlog normalize DESeq object
tiles1000bp_ddsrlog <- rlog(tiles1000bp_dds, blind=TRUE)
# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(tiles1000bp_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- colData$sample 
#colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(tiles1000bp_dds)[,c("celltype","tn5")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap spearman

#Plot Spearman correlation of counts under peaks 
#import counts matrix into DEseq object, use colData from above
#peakcounts
peakcounts_dds <- DESeq2::DESeqDataSetFromMatrix(countData=opentn5fragment_counts, design= celltype~tn5, colData= colData)
peakcounts_ddsrlog <- rlog(peakcounts_dds, blind=TRUE) #rlog transform
save(peakcounts_ddsrlog, file = "Counts_underPeaks_RLOG_techreps_ChromVARfragmentcounts.RData")

# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(peakcounts_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- colData$sample 
#colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(peakcounts_dds)[,c("celltype","tn5","Lot")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
write.csv(cor.matrix, file = "Spearman_CorMatrix_underPeaks_RLOG_techreps.csv")
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap

##### NEW PLOT 
#plot K562 sample correlation with transposome amount varied
setwd("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/K562_transposomeamountvaried")
#make vector of bamfile paths
Libs = list.files("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/K562_transposomeamountvaried")
Libs = Libs[grep(".bam$", Libs)]
bamfiles <- Libs
bamfiles

# count fragments under each 1000bp bin from bam files
Counts_1000bptiles <- chromVAR::getCounts(paste0(bamfiles), hg38_1000bptiles, paired=TRUE, by_rg=FALSE, format="bam")
Counts_1000bptiles <- assay(Counts_1000bptiles) # convert into counts matrix, matrix object
# add row names based on chromosome coordinates
row.names(Counts_1000bptiles) <- paste0(seqnames(hg38_1000bptiles), ":", start(hg38_1000bptiles), "-", end(hg38_1000bptiles))
save(Counts_1000bptiles, file = "Counts_1kbptiles_K562all_ChromVAR.RData")
#add metadata to matrix with coldata dataframe
#setup metadat for experiment for input into deseq
colData <- matrix(data = 0,nrow = 14) 
colData <- as.data.frame(colData)
rownames(colData) <- colnames(Counts_1000bptiles)

# add info about experiment into colData metadata table
colData$tn5 <-c(rep("Illumina",3),rep("pG-Tn5",11))
colData$Amount <-c(rep("2.5uL",3),rep("1uL",2),rep("5uL",3),rep("7.5uL",3),rep("1uL",3))

#import count object into DESeq2
#1000bp
tiles1000bp_dds <- DESeq2::DESeqDataSetFromMatrix(countData=Counts_1000bptiles, design= ~tn5, colData= colData)
# rlog normalize DESeq object
tiles1000bp_ddsrlog <- rlog(tiles1000bp_dds, blind=TRUE)
save(tiles1000bp_ddsrlog, file = "RLOGCounts_1kbptiles_K562all_ChromVAR.RData")
# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(tiles1000bp_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- colData$sample 
#colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(tiles1000bp_dds)[,c("Amount","tn5")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
write.csv(cor.matrix, file = "Spearman_CorMatrix_1kbtiles_K562all.csv")
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap spearman


#2024-07-02 double checking things before submission

peaks <-getPeaks("/rugpfs/fs0/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/GEO_submission/ATAC_opentn5masterpeakset_MCF7andK562.merged.bed",sort_peaks=TRUE)
resize(peaks, width = 500, fix = "center")




#make vector of bamfile paths
setwd("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/techrepbams")
Libs = list.files("/lustre/fs4/risc_lab/scratch/landerson01/20230822_OpenTn5_MCF7_K562/techrepbams")
Libs = Libs[grep(".bam$", Libs)]
bamfiles <- Libs
bamfiles

opentn5fragment_counts <- getCounts(bamfiles, 
                                    peaks, 
                                    paired =  TRUE, 
                                    by_rg = FALSE, 
                                    format = "bam")

opentn5fragment_counts = addGCBias(opentn5fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
opentn5fragment_counts <- assay(opentn5fragment_counts)# convert into counts matrix, matrix object
save(opentn5fragment_counts, file = "opentn5_MCF7+K562_ChromVARfragmentcounts_20220702.RData")

#add metadata to matrix with coldata dataframe
#setup metadat for experiment for input into deseq
colData <- matrix(data = 0,nrow = 12) 
colData <- as.data.frame(colData)
rownames(colData) <- colnames(opentn5fragment_counts)

# add info about experiment into colData metadata table
colData$celltype <- c(rep("K562",6),rep("MCF7",6))
colData$tn5 <-c(rep("Illumina",3),rep("pG-Tn5",3),rep("Illumina",2),rep("pG-Tn5",4))
colData$Lot <-c(rep("Illumina",3),rep("pG-Tn5 Lot 1",3),rep("Illumina",2),rep("pG-Tn5 Lot 1",2),rep("pg-Tn5 Lot 2",2))

#Plot Spearman correlation of counts under peaks 
#import counts matrix into DEseq object, use colData from above
#peakcounts
peakcounts_dds <- DESeq2::DESeqDataSetFromMatrix(countData=opentn5fragment_counts, design= celltype~tn5, colData= colData)
peakcounts_ddsrlog <- rlog(peakcounts_dds, blind=TRUE) #rlog transform
save(peakcounts_ddsrlog, file = "Counts_underPeaks_RLOG_techreps_ChromVARfragmentcounts_20220702.RData")

# get spearman correlation coefficients for pairwise comparisons between samples
sampleDists <- dist(t(assay(peakcounts_ddsrlog))) # get rlog normalized count matrix & convert to distance measure for calculating correlations. 
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- colData$sample 
#colnames(sampleDistMatrix) <- colData$sample # specify col names as sample names 
heatmap_anno <- as.data.frame(colData(peakcounts_dds)[,c("celltype","tn5","Lot")]) #set annotations to include on heatmap
rownames(heatmap_anno) <- rownames(sampleDistMatrix)
cor.matrix <-cor(sampleDistMatrix, method="spearman") # measure correlation 
write.csv(cor.matrix, file = "Spearman_CorMatrix_underPeaks_RLOG_techreps_20220702.csv")
pheatmap(cor.matrix, annotation_col = heatmap_anno) # plot correlation heatmap


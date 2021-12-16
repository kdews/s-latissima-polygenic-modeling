# SETUP
# Remove everything from memory
rm(list=ls())
# Load packages
library(BiocManager)
# To install or update packages
# BiocManager::install("gdsfmt", force = TRUE)
# BiocManager::install("SeqArray", force = TRUE)
# BiocManager::install("SNPRelate", force = TRUE)
# install.packages("RColorBrewer")
library(SNPRelate)
library(SeqArray)
#library(RColorBrewer)
# Read script input
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
# Define VCF filename
vcf_file <- args[1]

# CONVERT INPUT
# Convert VCF to GDS file (if needed)
vcf_base <- sub('\\.vcf.*', '', (basename(vcf_file)))
gds_file <- paste(vcf_base, ".gds", sep="")
if(file.exists(gds_file)) {
  print(paste("Detected GDS file:", gds_file))
} else {
  seqVCF2GDS(vcf_file, gds_file, verbose=FALSE)
}

# RUN PCA ANALYSIS
# Open GDS file
genofile <- seqOpen(gds_file)
samp.id <- seqGetData(genofile, "sample.id")
# Filter out "S. angutissima" subspecies
# seqSetFilter(genofile, sample.id=samp.id[grep("SA",samp.id, invert=TRUE)])
pca <- snpgdsPCA(genofile, num.thread=12,
#                 sample.id=samp.id[grep("SA",samp.id, invert=TRUE)], 
                 autosome.only=FALSE)
tab <- data.frame(sample.id=pca$sample.id,
                  EV1=pca$eigenvect[,1],    # the first eigenvector
                  EV2=pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
# Close GDS file when done using to avoid bugs
seqClose(genofile)

# ANNOTATION
# Annotate populations
pop.group <- unlist(strsplit(tab$sample.id, "-"))
pop.group <- pop.group[seq(2, length(pop.group), 5)]
# Annotate sex
sex <- 1:length(tab$sample.id)
male_subset <- grep("-MG-",tab$sample.id)
female_subset <- grep("-FG-",tab$sample.id)
sex[male_subset] <- rep("male",length(male_subset))
sex[female_subset] <- rep("female",length(female_subset))
# # Annotate subspecies
# subspecies <- sub("SL.*", "latissima", tab$sample.id)
# subspecies <- sub("SA.*", "angutissima", subspecies)
# Add annotations to PCA dataframe
tab_annot <- data.frame(sample.id = pca$sample.id,
                        sex = factor(sex)[match(pca$sample.id, tab$sample.id)],
                        # subspecies = factor(subspecies)[match(pca$sample.id, 
                        #                                       tab$sample.id)],
                        pop = factor(pop.group)[match(pca$sample.id, tab$sample.id)],
                        EV1 = pca$eigenvect[,1],    # the first eigenvector
                        EV2 = pca$eigenvect[,2],    # the second eigenvector
                        stringsAsFactors = FALSE)
# Percent variance explained
perc_1 <- round(pca$varprop[1]*100, digits=2)
perc_2 <- round(pca$varprop[2]*100, digits=2)

# PLOTS
# Set color palettes
#brewer.pal(length(levels(tab_annot$pop)), "Set3")
pop_col <- data.frame(pop=tab_annot$pop,
                      col=rainbow(length(unique(tab_annot$pop)))[tab_annot$pop])
# PCA color by pop
png(file=paste(vcf_base, "_pca_pop.png", sep=""), 
    width = 1000, height = 1000, res = 110)
plot(tab$EV2, tab$EV1, col=pop_col$col, pch=16,
     xlab=paste("PC2 = ", perc_2, "%", sep=""), 
     ylab=paste("PC1 = ", perc_1, "%", sep=""))
legend("topright", legend=unique(pop_col)$pop, pch=16, 
       col=unique(pop_col)$col)
dev.off()
# PCA color by sex
png(file=paste(vcf_base, "_pca_sex.png", sep=""),
    width = 1000, height = 1000, res = 110)
plot(tab$EV2, tab$EV1, col=as.integer(tab_annot$sex), pch=16,
     xlab=paste("PC2 = ", perc_2, "%", sep=""), 
     ylab=paste("PC1 = ", perc_1, "%", sep=""))
legend("topright", legend=levels(tab_annot$sex), pch=16, 
       col=1:nlevels(tab_annot$sex))
dev.off()
# PCA color by subspecies 
# png(file="s_latissima_popgen_all_variants_no_checks_subspecies.png", 
#     width = 1000, height = 1000, res = 110)
# plot(tab$EV2, tab$EV1, col=as.integer(tab_annot$subspecies), pch=16,
#      xlab=paste("PC2 = ", perc_2, "%", sep=""), 
#      ylab=paste("PC1 = ", perc_1, "%", sep=""))
# legend("topright", legend=levels(tab_annot$subspecies), pch=16, 
#        col=1:nlevels(tab_annot$subspecies))
# dev.off()


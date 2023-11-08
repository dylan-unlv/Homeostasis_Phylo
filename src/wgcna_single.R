library(WGCNA)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)


#####
#read in data, separate into data and metadata
####

tdata <- read_delim('data/TPM_5-species.tsv', delim='\t')
genes <- tdata$Gene
tdata <- data.frame(tdata %>% select(-Gene))
rownames(tdata) <- genes
meta <- t(tdata) %>% rownames() %>% str_split_fixed(string = ., pattern='\\.', 4) %>% as_tibble()
colnames(meta) <- c('species', 'temp', 'rep', 'licate')
meta$temp <- as.integer(meta$temp)
meta$id <- colnames(tdata)
data <- t(tdata)

#####
#subset data -- single network for each species (probably need to look for types of camel)
####

species <- 'rat'
data <- data[grep(species, rownames(data)), ]

#make sure meta and data are in the same order
meta <- meta %>% filter(id %in% (data %>% rownames())) %>% arrange(match(id, data %>% rownames()))

#####
# Normalize counts (TPM needs log2(x+1))
####

#filter out genes with 0 counts in 90%+ of samples (1783/9825 ~ 20% in camel)
data <- data[,colSums(data>0)/nrow(data)>0.9]
#normalize
data <- log2(data +1)


#####
# Pick Soft Threshold (0.9 is optimal to keep scale-free topology, 0.8 admissible)
####

#good to check no NA, but unlikely at this step
#goodSamplesGenes(mat, verbose=3)

#check no outlier sampoles
sampleTree <- hclust(dist(data), method = "average")

#save the sample clusters
png(filename = paste0('figs/',species,'_sample_clusters.png'), width=600, height=900, units = 'px')

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h=80, col='red')
dev.off()

#remove outliers if any 
# Determine cluster under the line

#clust <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples <- (clust==1)
#data <- data[keepSamples, ]

#create df in same order as mat, save
#datTraits <- meta[match(rownames(mat), meta$sample),] %>% select(-c(sample))
#nums <- as.numeric(str_sub(datTraits$temp, end=-2))
#nums[is.na(nums)] <- 37
#datTraits$temp <- nums 


#set soft-thresholding power for network generation
powers <- c(c(1:10), seq(from = 12, to=66, by=2))
sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5)

#save thresholding results 
png(filename = paste0('figs/',species,'_sft_soft_threshold.png'), width=700, height=700, units = 'px')

cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
# this line corresponds to using an R^2 cut-off of .9
abline(h=0.90,col="red")
dev.off()



#softPower for each species:
#camel=20
#human=14
#rat=48
#rhino=10
#squirrel=58

#####
# Generate networks
#####

nSamples <- nrow(meta)
softPower <- 48
adj <- adjacency(data, power=softPower)
TOM <- TOMsimilarity(adj)
dissTOM <- 1-TOM

#module detection using dynamic tree cutting
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 30 
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

#this function outputs a ton of modules, we have to merge coexpressed modules
dynamicColors <- labels2colors(dynamicMods)
MEList <- moduleEigengenes(data, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

#based on the threshold below, merge modules by eigengenes
MEDissThres <- 0.05 #0.25 camel human rhino; 0.05 rat squirrel
merged <- mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merged$colors
mergedMEs <- merged$newMEs


#plot merge
pdf(file = paste0('figs/',species,'_merged_clusters.pdf'), wi=9, he=12)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#final merge
moduleColors <- mergedColors
colorOrder<-c("grey", standardColors(50));
moduleLabels<-match(moduleColors, colorOrder)-1;
MEs <- mergedMEs

#recalculate eigengenes with merged modules
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(data, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, meta$temp, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#look at correlations within the module eigengenes to temp
pdf(file = paste0('figs/',species,'_eigengene_temp_corr.pdf'), wi=5, he=9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot (if many meta values, change xlabels)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = 'temp',
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.y = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships "))
dev.off()


#####
#save progress
####
saveRDS(list(TOM=TOM, data=data, MEs=MEs, meta=meta, moduleColors=moduleColors),
        file = paste0('data/wgcna/',species,'_wgcna_objs.rds'))


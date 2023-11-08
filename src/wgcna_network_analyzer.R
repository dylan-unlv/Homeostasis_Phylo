library(tidyverse)
library(biomaRt)
library(WGCNA)
library(ComplexHeatmap)


cons_data <- read_delim('data/cytoscape/pent_consensus.nodes', delim='\t')

#####
#First, KEGG pathway analysis of modules
#


#need to convert to entrez ids for this
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
altNames <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id"),
                  values=cons_data$nodeName,mart= mart,uniqueRows = F)


library(clusterProfiler) #needs to be import after the previous line... seriously.

#iter modules
mods <- cons_data$module %>% unique()
for (imodule in mods){
  
  #subset data, add entrez ids
  data <- cons_data %>% filter(module==imodule) %>% 
    rename('ensembl_gene_id'=nodeName) %>% 
    dplyr::left_join(., altNames, by='ensembl_gene_id')


  #KEGG / GO Terms
  result <- enrichKEGG(data$entrezgene_id, keyType = 'kegg',
            organism = 'hsa', pAdjustMethod = 'hochberg')

  #write out
  result@result %>% write_csv(paste0('data/kegg/',imodule,'_KEGG.csv'))
}



#####
#Correlate consensus modules with metadata
#
cons_R <- readRDS('data/wgcna/pent_consensus_wgcna_objs.rds')
setLabels <- c('camel','human','rat','rhino','squirrel')
nSets <- 5 
tdata <- read_delim('data/TPM_5-species.tsv', delim='\t')
meta <- t(tdata %>% select(-Gene)) %>% rownames() %>% str_split_fixed(string = ., pattern='\\.', 4) %>% as_tibble()
colnames(meta) <- c('species', 'temp', 'rep', 'licate')
meta$temp <- as.integer(meta$temp)
meta$id <- colnames(tdata %>% select(-Gene))


#get module eigengenes

for (setn in 1:nSets){
  tmeta <- meta %>% filter(species==setLabels[setn])
  MEs0 <- moduleEigengenes(cons_R$multiExpr[[setn]]$data, cons_R$moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, tmeta$temp, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(tmeta))

  if (setn==1) {
    gmat <- moduleTraitCor
    pmat <- moduleTraitPvalue
  }
  else{
    gmat <- cbind(gmat, moduleTraitCor[match(rownames(gmat), rownames(moduleTraitCor))] )
    pmat <- cbind(pmat, moduleTraitPvalue[match(rownames(pmat), rownames(moduleTraitPvalue))])
  }
  
}

setLabelsordered<- c('camel','rhino','human','squirrel','rat')

gmat <- data.frame(gmat)
pmat <- data.frame(pmat)
colnames(gmat) <- setLabels
colnames(pmat) <- setLabels
gmat <- as.matrix(gmat[setLabelsordered])
pmat <- as.matrix(pmat[setLabelsordered])

#look at correlations within the module eigengenes to temp
pdf(file = paste0('figs/consensus_5_eigengene_temp_corr.pdf'), wi=9, he=9)
# Will display correlations and their p-values
textMatrix = paste(signif(gmat, 2), "\n(",
                   signif(pmat, 1), ")", sep = "")
dim(textMatrix) = dim(gmat)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot (if many meta values, change xlabels)
labeledHeatmap(Matrix = gmat,
               xLabels = setLabelsordered,
               yLabels = rownames(gmat),
               ySymbols = rownames(gmat),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.y = 0.9,
               zlim = c(-1,1),
               main = paste("Consensus Module-Temp Relationships "))
dev.off()
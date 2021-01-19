library(DESeq2)
library(readr)
library(dplyr)

gene_descriptions <- read.delim("~/Google Drive/Microbiome/FISSEQ/data/Arabidopsis/gene_descriptions.txt")
mRNAcounts <- data.frame(read_delim("~/Google Drive/Microbiome/FISSEQ/data/Arabidopsis/Arabidopsis_RNAseq.txt", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE))

mRNAcounts           <- mRNAcounts[-1,]
rownames(mRNAcounts) <- mRNAcounts$X1
mRNAcounts           <- mRNAcounts[,-1]

#focus on 2 exp
df.f           <- mRNAcounts[,c(4:8,14:18)]
df.f           <- df.f %>% mutate_if(is.character,as.numeric) 
rownames(df.f) <- rownames(mRNAcounts)
colData        <- DataFrame(condition=factor(c("CL28","CL28","CL28","CL28","CL28",
                                        "both", "both","both", "both", "both")))
# create DESeq input matrix
dds <- DESeqDataSetFromMatrix(df.f, colData,
                              formula(~ condition))
# run DEseq
dds <- DESeq(dds)

plotMA(dds)
# get differentially expressed genes
res             <- results(dds)
res$description <- sapply(rownames(res), function(x) gene_descriptions$Gene.Model.Description[which(gene_descriptions$Locus.Identifier == x)])
# order by BH adjusted p-value
resOrdered      <- res[order(res$padj),]
# top of ordered matrix
head(resOrdered)

# get differentially expressed gene matrix
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    abs(resOrdered$log2FoldChange)>=1,]


write.csv(data.frame(sig),file = 'DEG_CL28vsboth.csv')


write.csv(data.frame(sig[VennCommon$...1,]),file = 'DEG_CL28vsboth_filtered.csv')

#Plot heatmap
library(pheatmap)
library(RColorBrewer)
dfx <- mRNAcounts[rownames(sig),]
gene.names <- rownames(dfx)
dfx <- dfx %>% mutate_if(is.character,as.numeric) 
rownames(dfx) <- gene.names
colnames(dfx) <- c('NB','NB','NB','C28','C28','C28','C28','C28',
                   'C14','C14','C14','C14','C14',
                   'both','both','both','both','both')
pheatmap(dfx,
         color = rev(brewer.pal(n = 9,name = 'PuOr')),border_color = 'white',
         fontsize_row = 3,
         scale = 'row')

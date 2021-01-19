#Analysis of Arabidopsis RNA seq from https://doi.org/10.1038/s41586-020-2778-7

#### 0. Load required packages ####
library(DESeq2)
library(readr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#### 1. Load data ####
gene_descriptions <-
  read.delim("data/Arabidopsis/gene_descriptions.txt")
mRNAcounts <-
  data.frame(
    read_delim(
      "data/Arabidopsis/Arabidopsis_RNAseq.txt",
      "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    )
  )

mRNAcounts           <- mRNAcounts[-1, ]
rownames(mRNAcounts) <- mRNAcounts$X1
mRNAcounts           <- mRNAcounts[, -1]

#### 2. Analysis focus on comparing CL28 vs. [CL28+CL14] ####
#select right columns
df.f           <- mRNAcounts[, c(4:8, 14:18)]
df.f           <- df.f %>% mutate_if(is.character, as.numeric)
rownames(df.f) <- rownames(mRNAcounts)
#rename columns
colData        <-
  DataFrame(condition = factor(
    c(
      "CL28",
      "CL28",
      "CL28",
      "CL28",
      "CL28",
      "both",
      "both",
      "both",
      "both",
      "both"
    )
  ))
# create DESeq input matrix
dds <- DESeqDataSetFromMatrix(df.f, colData,
                              formula( ~ condition))
# run DEseq
dds <- DESeq(dds)

# get diferentially expressed genes
res             <- results(dds)
res$description <-
  sapply(rownames(res), function(x)
    gene_descriptions$Gene.Model.Description[which(gene_descriptions$Locus.Identifier == x)])
# order by BH adjusted p-value
resOrdered      <- res[order(res$padj), ]
# top of ordered matrix
head(resOrdered)

# get differentially expressed gene matrix and filter for significant genes
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj < 0.10 &
                    abs(resOrdered$log2FoldChange) >= 1, ]

#Plot heatmap
dfx <- mRNAcounts[rownames(sig), ]
gene.names <- rownames(dfx)
dfx <- dfx %>% mutate_if(is.character, as.numeric)
rownames(dfx) <- gene.names
colnames(dfx) <- c(
  'NB',
  'NB',
  'NB',
  'C28',
  'C28',
  'C28',
  'C28',
  'C28',
  'C14',
  'C14',
  'C14',
  'C14',
  'C14',
  'both',
  'both',
  'both',
  'both',
  'both'
)
plot0 <- pheatmap(
  dfx,
  color = rev(brewer.pal(n = 9, name = 'PuOr')),
  border_color = 'white',
  fontsize_row = 3,
  fontsize_col = 4,
  scale = 'row'
)


#Save output and heatmap
write.csv(data.frame(sig),file = 'res/Arabidopsis/DEG_CL28vsboth.csv')
ggsave(plot0,filename = 'res/Arabidopsis/heatmap.pdf')


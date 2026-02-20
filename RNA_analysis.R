# These codes are taken and edited from FGT_T7_Advanced_Analysis_RNASeq_Using R

library(tidyverse)
library(DESeq2)
library(scatterplot3d)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)

# Load in the feature table
load(file = "./mytable_features")

count_table <- mytable_features$counts
annotation_table <- mytable_features$annotation
stat_table <- mytable_features$stat
targets <- mytable_features$targets

colnames(count_table) <- c("Control_1", "Control_2", "Control_3", "KO_1", "KO_2", "KO_3")

sample_info <- data.frame(row.names = c("Control_1", "Control_2", "Control_3", "KO_1", "KO_2", "KO_3"),
                          condition = c("control", "control", "control", "knockout", "knockout", "knockout"),
                          colour = c("red", "red", "red", "blue", "blue", "blue"))

# DESeq2 object using count matrix and metadata
dds_rnaseq <- DESeqDataSetFromMatrix(countData = count_table,
                                     colData = sample_info,
                                     design = ~ condition)

dim(dds_rnaseq)

slotNames(dds_rnaseq)

colnames(colData(dds_rnaseq))

colnames(assay(dds_rnaseq))


# Quality filter to filter out 3 sample with 10 coutns
keep <- rowSums(counts(dds_rnaseq) >= 10) >= 3

dds_rnaseq <- dds_rnaseq[keep,]

nrow(dds_rnaseq)


# Data normalisation
# Calculate size factors
dds_rnaseq <- estimateSizeFactors(dds_rnaseq)
sizeFactors(dds_rnaseq)

# Normalised
rld_rnaseq <- rlog(dds_rnaseq, blind = T)

boxplot(assay(rld_rnaseq), main = "Data normalisation of sample count",
        xlab = "Sample", ylab="Normalised expression level")


# Exploratory data plots through PCA

# 3d PCA
pca <- prcomp(t(na.omit(assay(rld_rnaseq))), scale=T)


s3d<-scatterplot3d(pca$x[,1:3], pch=19,
                   color= dds_rnaseq$colour,
                   main = "PCA analysis of samples")

s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(rld_rnaseq),
     pos = 3, offset = 0.5, cex=0.5)


# 2d PCA
pca_data <- as.data.frame(pca$x)
pca_data$colour = c("red", "red", "red", "blue", "blue", "blue")
pca_data$Sample <- rownames(pca_data)

ggplot(pca_data, aes(x = PC1, y = PC2, color = colour)) +
  geom_point(size = 2, shape = 19) +
  geom_text(aes(label = Sample), vjust = 1.5, hjust = 0.5, size = 3) +
  theme_minimal() +
  labs(x = "PC1", y = "PC2", title = "PCA analysis") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

#We can also output a distance matrix…
sampleDists <- dist(t(assay(rld_rnaseq)))
sampleDists

# Heatmap
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix)


#easier PCA plot
plotPCA(rld_rnaseq, intgroup = "condition")


# Differential gene experession
dds_rnaseq <- DESeq(dds_rnaseq)


result_condition <- results(dds_rnaseq,contrast=c("condition","control","knockout"))

plotMA(result_condition, main="MA plot of Control vs KnockOut", ylim=c(-2,2))
table(result_condition$padj < 0.01)
table(result_condition$padj < 0.05)

result_condition_selected <- subset(result_condition, padj < 0.05)

result_condition_selected <- result_condition_selected[order(abs(result_condition_selected$log2FoldChange),
                                                 decreasing = TRUE), ]

head(result_condition_selected, 10)


# Choose the top 10
top10 <- rownames(result_condition_selected)[1:10]

pheatmap(assay(rld_rnaseq)[top10,], scale="row", show_rownames=T)


# Biomart for common gene name

ensembl_host <- "https://www.ensembl.org"

head(biomaRt::listMarts(host = ensembl_host), 15)

head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",
                                                 mart = useMart("ENSEMBL_MART_ENSEMBL",
                                                host = ensembl_host))), 40)

mart <- biomaRt::useDataset(dataset = "mmusculus_gene_ensembl",
                            mart = useMart("ENSEMBL_MART_ENSEMBL",
                            host = ensembl_host))

resultAnnot <- biomaRt::getBM(values=rownames(dds_rnaseq),
                              attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","description","strand"),
                              filters="ensembl_gene_id", mart=mart)

names <- resultAnnot[,1]

resultAnnot <- as.data.frame(resultAnnot)

rownames(resultAnnot) = names

idx<-match(rownames(dds_rnaseq),rownames(resultAnnot))

all(rownames(dds_rnaseq) == rownames(resultAnnot))
grr<-resultAnnot[match(rownames(dds_rnaseq), resultAnnot$ensembl_gene_id),]

all(rownames(dds_rnaseq) == rownames(grr))
resultAnnot <- grr
all(rownames(dds_rnaseq) == rownames(resultAnnot))

nice_names<- paste(resultAnnot$ensembl_gene_id,resultAnnot$external_gene_name, sep = '_')
resultAnnot$nice_names <- nice_names
head(resultAnnot)
all(rownames(dds_rnaseq) == rownames(resultAnnot))

rld_rnaseq <- rlog(dds_rnaseq, blind = TRUE)
idx2 <- match(rownames(result_condition_selected)[1:10],rownames(dds_rnaseq))
plotme <-(rld_rnaseq)[rownames(result_condition_selected)[1:10],]
rownames(plotme)<-resultAnnot$nice_names[idx2]

#png(filename="heatmap_top10.png",width=500, height=500)
pheatmap(assay(plotme),scale="row",fontsize_row = 10,cellheight =12,
         cellwidth =12,treeheight_row = 40, treeheight_col = 40)
#dev.off()


# Functional enrichment analysis
result <- as.data.frame(result_condition_selected)
result_10 = result[1:10,]
result_10 <- rownames_to_column(result_10, var = "X")

mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
ens2entrez <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), mart=mart)
head(ens2entrez)

results_10_2 <- inner_join(result_10, ens2entrez, by=join_by("X"=="ensembl_gene_id"))

head(results_10_2)

#remove any that are NA
results_10_2 <-results_10_2[is.na(results_10_2$entrezgene_id)==FALSE,]
dim(results_10_2)

head(results_10_2)

ego_BH <- enrichGO(gene = results_10_2$entrezgene_id,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",       
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)       

dotplot(ego_BH) + ggtitle("GO Enrichment Analysis (Biological Process)")


ego_MF <- enrichGO(gene = results_10_2$entrezgene_id,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF", 
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)       

dotplot(ego_MF) + ggtitle("GO Enrichment Analysis (Molecular Function)")


ego_CC <- enrichGO(gene = results_10_2$entrezgene_id,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)          

dotplot(ego_CC) + ggtitle("GO Enrichment Analysis (Cellular Component)")
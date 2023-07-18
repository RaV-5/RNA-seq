#install.packages("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
library('tidyverse')

#count_data <- read.csv('TCGA-count-data.csv', header = TRUE, row.names = 1)

Count_data = read.table(file = "TCGA-count-data.csv", header = T, sep = ",",row.names=1,check.names = FALSE)

dim(Count_data)

Col_data = read.table(file = "TCGA-column-data.csv", header = T, sep = ",",row.names = 1)

rownames(Col_data)

colnames(Count_data)

all(rownames(Col_data)==colnames(Count_data))

#boxplot(Count_data)
hist(Count_data[,1]) # Plotting only the first sample (column 1)

#count no of NA values in matrix
class(Count_data)
(is.na(Count_data))

which(is.na(Count_data),arr.ind=TRUE)

sum(is.na(Count_data))
#replace missing values in matrix with rowsums
#install.packages("zoo")
library(zoo)
Count_data[]<-t(na.aggregate(t(Count_data)))
Count_data
#replace NA values in matrix with zero
##Count_data[is.na(Count_data)] <- 0
#removing genes with all zero values
##df1[rowSums(df1[])>0,]
#?DESeq
dds = DESeqDataSetFromMatrix(countData = round(Count_data),
                             colData = Col_data,
                             design = ~ condition) # we're testing for the different condidtions
dds$condition <- relevel(dds$condition, ref = "normal")
dds
dds <- DESeq(dds)
res1 <- results(dds)
summary(res1)

###keep only sig results, padj<0.05 and log2FoldChange >1
resSigUp <- subset(res1, padj < 0.05 & log2FoldChange >1)
write.csv(resSigUp, "Upregulatedd.csv")

###keep only sig results, padj<0.05 and log2FoldChange < -1
resSigDown <- subset(res1, padj < 0.05 & log2FoldChange < -1)
write.csv(resSigDown, "Downregulated.csv")

###keep UP and Down in one file with padj<0.05
resSig <- subset(res1, log2FoldChange >1 & padj <0.05 | log2FoldChange < -1 & padj < 0.05)
write.csv(resSig, "DE.csv")

##Write for volcano plot
##write.table(res1, file="Volcano", row.names=F, sep="\t")

test <- res1 %>%
  filter(padj<0.05)

#Enhanced Volcano
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

##rpl <- read.table("Volcano", header = TRUE, sep = "\t")
#x<- res$log2FoldChange
#y<- res$pvalue

EnhancedVolcano(res1,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL,
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )

#Enhanced Volcano, how to save to pdf or png file
png(file="volcano.png")
#png(file="volcano.png",res=72)
pdf(file="volcano.pdf")
EnhancedVolcano(res1,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL,
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                transcriptPointSize = 0.5,
                transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )
dev.off()

#Enhanced Volcano with labels of important genes
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('ENSG00000287048.1','ENSG00000280032.1'),
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                border = 'full' )


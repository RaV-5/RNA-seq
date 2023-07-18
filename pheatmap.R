#https://rdrr.io/bioc/COMPASS/man/pheatmap.html
# Simple heatmaps using a matrix data
install.packages("pheatmap")
library(pheatmap)
getwd()
setwd("/Users/apple/Downloads/RNA-SEQ-Workshop/Heatmap/")
##write.csv(test, "heatmap_data2.csv")
pheat <- read.csv("heatmap_data2.csv", header = T, row.names = 1, check.names = FALSE)
pheat <- as.matrix(pheat)
pheatmap(pheat)
pheatmap (pheat, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap (pheat, cluster_cols = FALSE)
pheatmap(pheat, cluster_row = FALSE)
pheatmap(pheat, kmeans_k = 4)
pheatmap(pheat, scale = "row", clustering_distance_rows = "correlation")
pheatmap(pheat, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(pheat, legend = FALSE)
pheatmap(pheat, display_numbers = TRUE)
pheatmap(pheat, display_numbers = TRUE, number_format = "%.1e")
pheatmap(pheat, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
                                                                            "1e-4", "1e-3", "1e-2", "1e-1", "1"))
pheatmap(pheat, cellwidth = 15, cellheight = 12, main = "Example heatmap")
pheatmap(pheat, cellwidth = 15, cellheight = 12, fontsize = 8)

###Advanced heatmaps with annotations
###https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq")
#library("DESeq")
library("pheatmap")
install.packages("dendextend")
library("dendextend")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#example_file <- system.file ("extra/TagSeqExample.tab", package="DESeq")
#data <- read.delim(example_file, header=T, row.names="gene")
#data_subset <- as.matrix(data[rowSums(data)>50000,])
# create heatmap using pheatmap
pheatmap(pheat)
#The heatmap.2() function has a parameter for scaling the rows; this can be easily implemented.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(pheat, 1, cal_z_score))
##perform hierarchical clustering in the same manner as performed by pheatmap to obtain gene clusters.
my_hclust_gene <- hclust(dist(pheat), method = "complete")
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
#We can form two clusters of genes by cutting the tree with the cutree() function; we can either specific the height to cut the tree or the number of clusters we want.
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
##The 1's and 2's indicate the cluster that the genes belong to
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
##We can add multiple row annotations to a heatmap and below I'll add some random annotations.
set.seed(1984)
my_random <- as.factor(sample(x = 1:2, size = nrow(my_gene_col), replace = TRUE))
my_gene_col$random <- my_random
##I'll add some column annotations and create the heatmap.
my_sample_col <- data.frame(sample = rep(c("tumour", "normal"), c(7,3)))
row.names(my_sample_col) <- colnames(pheat)
pheatmap(pheat, annotation_row = my_gene_col, annotation_col = my_sample_col)
##introduce breaks in the heatmap. I'll break up the heatmap by specifying how many clusters I want from the dendrograms
pheatmap(pheat,
         annotation_row = my_gene_col,
         annotation_col = my_sample_col,
         cutree_rows = 2,
         cutree_cols = 2)
my_heatmap <- pheatmap(pheat,
                       annotation_row = my_gene_col,
                       annotation_col = my_sample_col,
                       cutree_rows = 2,
                       cutree_cols = 2)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "my_heatmap.png")
png(my_heatmap,"heatmap.png")

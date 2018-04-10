#!/usr/bin/env Rscript

# This script takes a loom object, builds a model to predict separability, and
# outputs the separability predictions for 1k-50k cells (by 1k)

################################ Arguments #####################################
# args[1]    Path to loom file
# args[2]    Name of cluster 1 
# args[3]    Fraction of cluster 1 in the full dataset
# args[4]    Name of cluster 2 
# args[5]    Fraction of cluster 2 in the full dataset
# args[6]    Path for output tsv file
################################################################################

suppressMessages(library(Seurat))
suppressMessages(library(loomR))

################################################################################

args <- commandArgs(trailingOnly = TRUE)
lf <- connect(filename = args[1])

# Get UMI matrix and cluster information from the loom object
umi.matrix <- t(lf[["matrix"]][,])
rownames(umi.matrix) <- lf[["row_attrs/gene_names"]][]
colnames(umi.matrix) <- lf[["col_attrs/cell_names"]][]
clusters <- lf[["col_attrs/cluster"]][]

seurat.object <- CreateSeuratObject(raw.data = umi.matrix)
seurat.object <- SetIdent(object = seurat.object,
                          cells.use = seurat.object@cell.names,
                          ident.use = clusters)

cluster.avgs <- AverageExpression(object = seurat.object, display.progress = FALSE)

cluster_1_name = args[3]
cluster_2_name = args[5]
# Find DE genes as those with > log(2) average difference between the two clusters
num.de.genes <- length(which(abs(cluster.avgs[, cluster_1_name] - cluster.avgs[, cluster_2_name]) > log(2)))

# Simple prediction based soley on DE genes
num.cells <- seq(1000, 100000, 1000)
predicted.separability <- 1 / (1 + exp(-num.cells * num.de.genes / 1e7))
output <- data.frame(num.cells, predicted.separability)
colnames(output) <- c("Number of Cells", "Predicted Separability")
write.table(x = output, file = args[6], sep = "\t", quote = FALSE, row.names = FALSE)

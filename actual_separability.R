#!/usr/bin/env Rscript

# This script takes a loom object and two cluster names and computes the actual
# separability

################################ Arguments #####################################
# args[1]    Path to loom file
# args[2]    Name of cluster 1
# args[3]    Name of cluster 2
################################################################################

suppressMessages(library(Seurat))
suppressMessages(library(RANN))
suppressMessages(library(loomR))

################################ Fxn Definitions ###############################
ProcessData <- function(object) {
  object <- FindVariableGenes(object, do.plot = FALSE, display.progress = FALSE)
  object@var.genes <- rownames(head(object@hvg.info, 1000))
  object <- ScaleData(object, genes.use = object@var.genes,
                      display.progress = FALSE, check.for.norm = FALSE)
  object <- RunPCA(object, pcs.compute = 30, do.print = FALSE)
  return(GetCellEmbeddings(object, reduction.type = "pca", dims.use = 1:30))
}

ComputeSeparability <- function(input.data, cells.1, cells.2, k = 20) {
  if (length(cells.1) < 3) {
    stop("Fewer than 3 cells in the first group (cells.1)")
  }
  if (length(cells.2) < 3) {
    stop("Fewer than 3 cells in the second group (cells.2)")
  }
  k <- min(c(k, length(cells.1), length(cells.2)))
  tnn <- nn2(data = input.data[c(cells.1, cells.2), ], k = k + 1)
  idx <- tnn$nn.idx[, -1]
  rownames(idx) <- c(cells.1, cells.2)
  correct_neighbors_c1 <- sapply(cells.1, function(x)
    {
      length(which(rownames(idx)[idx[x,]] %in% cells.1))
    }
  )
  correct_neighbors_c2 <- sapply(cells.2, function(x)
    {
      length(which(rownames(idx)[idx[x,]] %in% cells.2))
    }
  )
  return(mean(c(correct_neighbors_c1, correct_neighbors_c2)) / k)
}
################################################################################

args <- commandArgs(trailingOnly = TRUE)
lf <- connect(filename = args[1])

# Get UMI matrix and cluster information from the loom object
umi.matrix <- t(lf[["matrix"]][,])
rownames(umi.matrix) <- lf[["row_attrs/gene_names"]][]
colnames(umi.matrix) <- lf[["col_attrs/cell_names"]][]
clusters <- lf[["col_attrs/cluster"]][]

seurat.object <- CreateSeuratObject(raw.data = umi.matrix)
seurat.object <- NormalizeData(object = seurat.object, display.progress = FALSE)
processed.mat <- ProcessData(object = seurat.object)
seurat.object <- SetIdent(object = seurat.object,
                          cells.use = seurat.object@cell.names,
                          ident.use = clusters)
cells.1 <- WhichCells(object = seurat.object, ident = args[2])
cells.2 <- WhichCells(object = seurat.object, ident = args[3])
separability <- ComputeSeparability(input.data = processed.mat,
                    cells.1 = cells.1,
                    cells.2 = cells.2)
print(paste0("Separability: ", separability))


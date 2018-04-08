suppressMessages(library(Seurat))
suppressMessages(library(loomR))

# This function will take in a loom object and downsample it to a specified number
# or fraction of cells and UMIs
#
# @param loom       Loom object to downsample
# @param filename   Filename for new loom object to create
# @param umis       Number/percent of UMIs to keep per cell
# @param cells      Randomly downsample to this many/percent of cells
# @param seed       Random seed
# @param verbose    Print messages/progress bars
# @param overwrite  Overwrite existing loom file if it exists

DownsampleLoom <- function(loom, filename, umis = NULL, cells = NULL, seed = 1,
                           verbose = TRUE, overwrite = FALSE) {
  set.seed(seed)
  if (!is.null(cells)) {
    if (length(cells) > 1) {
      stop("cells should be a single number specifying the number or percentage of cells to keep")
    }
    if (cells < 1) {
      cells <- cells * loom$shape[2]
    }
  }
  umi.matrix <- as(loom[["matrix"]][,], "dgCMatrix")
  cells.to.keep <- 1:nrow(umi.matrix)
  if (!is.null(cells)) {
    # randomly select cells to keep
    if (verbose) {
      message("Downsampling cells")
    }
    cells.to.keep <- sample(x = cells.to.keep, size = cells)
    umi.matrix <- umi.matrix[cells.to.keep, ]
  }
  if (!is.null(umis)) {
    if (length(umis) > 1) {
      stop("umis should be a single number specifying the number or percentage of UMIs per cell to keep")
    }
    if (umis < 1) {
      nUMI <- Matrix::rowSums(umi.matrix)
      umis <- round(umis * nUMI)
    }
    # downsample UMIs
    if (verbose) {
      message("Downsampling UMIs")
    }
    umi.matrix <- SampleUMI(data = t(umi.matrix), max.umi = umis, progress.bar = verbose)
  }
  if (verbose) {
    message("Writing new loom file")
  }
  ds.loom <- loomR::create(filename = filename,
                           data = umi.matrix,
                           cell.attrs = list(cell_names = loom[["col_attrs/cell_names"]][cells.to.keep],
                                             cluster = loom[["col_attrs/cluster"]][cells.to.keep]),
                           gene.attrs = list(gene_names = loom[["row_attrs/gene_names"]][]),
                           overwrite = overwrite)
  return(ds.loom)
}

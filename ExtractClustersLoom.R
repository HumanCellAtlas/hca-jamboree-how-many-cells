suppressMessages(library(loomR))

# This function will take in a loom object and downsample it to a specified number
# or fraction of cells and UMIs
#
# @param loom       Loom object to downsample
# @param filename   Filename for new loom object to create
# @param clusters   Vector of clusters to keep
# @param overwrite  Overwrites filename if it exists

ExtractClustersLoom <- function(loom, filename, clusters, overwrite = FALSE) {
  cells.to.keep <- which(loom[["col_attrs/cluster"]][] %in% clusters)
  new.loom <- subset.loom(x = loom, m = cells.to.keep, filename = filename, overwrite = overwrite)
  return(new.loom)
}

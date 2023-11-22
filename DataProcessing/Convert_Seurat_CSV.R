# This script serves two functions:

# 1.) Given two csv files (count matrix and metadata),
# combine them both into a single Seurat/rds object.

# 2.) Extract all the data (count matrix and metadata)
# from an existing Seurat/rds object to two csv files,
# which can then be imported into the script
# "Convert_scanpy_csv.ipynb". In that Jupyter notebook,
# the two csv files can be combined to form a
# single Scanpy/Adata object.


# Create a function that takes in two csv files and 
# makes a Seurat object
library(Seurat)

csv_to_seurat <- function(mtx_path, meta_path, transpose = FALSE,
                          umap = FALSE, cols_to_factorize = NULL,
                          save_name = NULL){
  
  # mtx = file path to count matrix
  # meta = file path to metadata
  # transpose = logical value (TRUE/FALSE) asking whether or not to
  #      transpose the count matrix after loading.
  # umap = logical value (TRUE/FALSE) asking whether or not
  #      to add UMAP coordiantes from the metadata to the Seurat object
  # cols_to_factorize = a vector of column names to turn into factor data
  #      Default value is NULL, meaning no vectors will be factorized
  # save_name = file path to save rds object, if specified
  
  # Load the count matrix, with genes as rows and cells as columns
  mtx <- read.csv(file = mtx_path, header = TRUE, sep = ',', row.names = 1)
  
  # If rows are actually cells, then transpose the matrix
  if(transpose == TRUE){
    mtx <- as.matrix(mtx)        # Turn into a matrix to be able to transpose
    mtx <- t(mtx)                # Transpose
    mtx <- as.data.frame(mtx)    # Turn into a data frame
  }
  
  # Load the metadata
  meta <- read.csv(file = meta_path, header = TRUE, sep = ',', row.names = 1)
  
  # Make sure there are the same number of cells in the count matrix
  # and in the metadata
  if( length(colnames(mtx)) != length(rownames(meta)) ){
    print("Number of cells in matrix is not equal to number of cells in metadata.")
  }
  
  # Make sure all the cells in the matrix (column names) are the same as
  # the cells in the metadata (row names)
  if( all(colnames(mtx) != rownames(meta)) ){
    stop("Cells in count matrix are not the same as the cells in metadata.")
  }
  
  # If both checks are passed, proceed to creating a Seurat object
  rds <- CreateSeuratObject(mtx, meta.data = meta)
  
  # Add UMAP coordinates, if specified
  if(umap == TRUE){
    
    # Get the UMAP coordinates, add row names
    umap <- as.matrix(rds@meta.data[, c("UMAP1", "UMAP2")])
    rownames(umap) <- rownames(rds@meta.data)
    
    # Add the coordinates
    rds[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = "RNA")
    rownames(rds@reductions$umap@cell.embeddings) <- rownames(rds@meta.data)
  }
  
  # Turn specified columns into factor data
  if(length(cols_to_factorize) > 0){
    for(col in cols_to_factorize){
      idx = which(col == colnames(rds@meta.data))           # Find column index
      rds@meta.data[, idx] <- factor(rds@meta.data[, idx])  # Factorize
    }
  }
  
  # Make sure the cell identities correspond to Seurat clustering
  #Idents(object = rds) <- "seurat_clusters"
  
  # Save rds object, if specified
  if(!is.null(save_name)){
    saveRDS(rds, save_name)
  }
  
  # Return rds object
  return(rds)
}


# Example
# RDSobject <- csv_to_seurat(mtx_path  = 'path/to/mtx.csv',
#                            meta_path = 'path/to/meta.csv',
#                            umap = TRUE,
#                            cols_to_factorize = c("Matching", "Seurat_clusters") )

csv_to_seurat(mtx_path = 'data/CSVfiles/MergedDataset_mtx.csv',
              meta_path = 'data/CSVfiles/MergedDataset_meta.csv',
              transpose = TRUE,
              umap = FALSE,
              save_name = 'data/RDSobjects/MergedData.rds')



# Create a function that takes a Seurat object and 
# outputs two csv files

rds_to_csv <- function(rds, mtx_path, meta_path, transpose = FALSE){
  
  # rds = Seurat object
  # mtx_path = file location/name to store count matrix
  # meta_path = file location/name to store metadata
  # transpose = a logical value, setting to TRUE will transpose the count matrix
  
  counts <- rds@assays$RNA@counts
  
  
  if(transpose == TRUE){
    counts <- as.matrix(counts)       # Turn into a matrix to be able to transpose
    counts <- t(counts)               # Transpose
    counts <- as.data.frame(counts)   # Turn into a data frame
  }
  
  # Write out the data
  write.csv(counts, file = mtx_path, quote = FALSE, row.names = TRUE)
  write.csv(rds@meta.data, file = meta_path, quote = FALSE, row.names = TRUE)
  
  # Nothing to return
}


# Example
# rds_to_csv(rds = RDSobject,
#            mtx_path  = 'path/to/mtx.csv',
#            meta_path = 'path/to/meta.csv',
#            transpose = TRUE)





rds_to_csv(integrated,
           mtx_path = 'data/CSVfiles/Integrated_mtx.csv',
           meta_path = 'data/CSVfiles/Integrated_meta.csv',
           transpose = TRUE)




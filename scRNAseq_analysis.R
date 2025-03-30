# Single-Cell RNA-Seq Analysis: Parkinson's Disease vs. Control
# This script performs preprocessing, clustering, and differential expression analysis

# Install and load required packages
install.packages("Seurat")
install.packages("Matrix")
install.packages("EnhancedVolcano")

library(Seurat)
library(Matrix)
library(EnhancedVolcano)

# Set directory where extracted files are saved
extracted_dir <- "filepath/"  # Update this with your actual path

# List all relevant files
matrix_files <- list.files(extracted_dir, pattern = "matrix.mtx$", full.names = TRUE)
barcode_files <- list.files(extracted_dir, pattern = "barcodes.tsv$", full.names = TRUE)
gene_files <- list.files(extracted_dir, pattern = "genes.tsv$", full.names = TRUE)

# Print to verify file paths
print("Matrix files:")
print(matrix_files)
print("Barcode files:")
print(barcode_files)
print("Gene files:")
print(gene_files)

# Function to load a single dataset
load_sample <- function(matrix_path, barcode_path, gene_path) {
  # Read count matrix
  count_matrix <- readMM(matrix_path)
  
  # Read genes and barcodes
  genes <- read.table(gene_path, header = FALSE, sep = "\t")
  barcodes <- read.table(barcode_path, header = FALSE, sep = "\t")
  
  # Assign row and column names
  rownames(count_matrix) <- genes$V1  # Assuming gene IDs are in column 1
  colnames(count_matrix) <- barcodes$V1
  
  return(count_matrix)
}

# Print filenames to verify sample assignments
print("Control 1 (HC) File:")
print(matrix_files[1])
print("Control 2 (HC) File:")
print(matrix_files[2])
print("Disease 1 (PD) File:")
print(matrix_files[3])
print("Disease 2 (PD) File:")
print(matrix_files[4])

# Load all four datasets
control_1 <- load_sample(matrix_files[1], barcode_files[1], gene_files[1]) # First HC sample
control_2 <- load_sample(matrix_files[2], barcode_files[2], gene_files[2]) # Second HC sample
disease_1 <- load_sample(matrix_files[3], barcode_files[3], gene_files[3]) # First PD sample
disease_2 <- load_sample(matrix_files[4], barcode_files[4], gene_files[4]) # Second PD sample

# Check dimensions
print("Sample dimensions:")
print(paste("Control 1:", paste(dim(control_1), collapse = " x ")))
print(paste("Control 2:", paste(dim(control_2), collapse = " x ")))
print(paste("Disease 1:", paste(dim(disease_1), collapse = " x ")))
print(paste("Disease 2:", paste(dim(disease_2), collapse = " x ")))

# Create Seurat objects
control_1_obj <- CreateSeuratObject(counts = control_1, project = "Control1")
control_2_obj <- CreateSeuratObject(counts = control_2, project = "Control2")
disease_1_obj <- CreateSeuratObject(counts = disease_1, project = "Disease1")
disease_2_obj <- CreateSeuratObject(counts = disease_2, project = "Disease2")

# Merge all objects
combined_obj <- merge(control_1_obj, y = c(control_2_obj, disease_1_obj, disease_2_obj), 
                      add.cell.ids = c("C1", "C2", "D1", "D2"))

# Check the final matrix size
print(paste("Combined object dimensions:", paste(dim(combined_obj), collapse = " x ")))

# Quality control and preprocessing
combined_obj <- NormalizeData(combined_obj, normalization.method = "LogNormalize", scale.factor = 10000)
combined_obj <- FindVariableFeatures(combined_obj, selection.method = "vst", nfeatures = 2000)
combined_obj <- ScaleData(combined_obj, features = VariableFeatures(combined_obj))

# Dimensionality reduction
combined_obj <- RunPCA(combined_obj, features = VariableFeatures(combined_obj))
ElbowPlot(combined_obj) # Visualize PCA components

# Clustering and visualization
combined_obj <- FindNeighbors(combined_obj, dims = 1:20)
combined_obj <- FindClusters(combined_obj, resolution = 0.5)
combined_obj <- RunUMAP(combined_obj, dims = 1:20)

# Visualization
DimPlot(combined_obj, reduction = "umap", group.by = "ident")

# Find marker genes for each cluster
markers <- FindAllMarkers(combined_obj, only.pos = TRUE)
head(markers)

# Visualize expression of key genes
FeaturePlot(combined_obj, features = "ENSG00000228696")

# Differential expression between conditions
combined_obj$condition <- ifelse(combined_obj$orig.ident %in% c("Control1", "Control2"), 
                                "Control", "Disease")
Idents(combined_obj) <- combined_obj$condition

deg_results <- FindMarkers(combined_obj, 
                          ident.1 = "Disease", 
                          ident.2 = "Control", 
                          logfc.threshold = 0.25, 
                          min.pct = 0.1)

# Save results
write.csv(deg_results, "DEGs_disease_vs_control.csv")

# Visualization of DEGs
EnhancedVolcano(deg_results,
                lab = NA,  # Removes gene labels
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Differential Expression: Disease vs. Control",
                subtitle = "Log2 Fold Change vs. Adjusted p-value",
                pCutoff = 0.05,
                FCcutoff = 0.5)

# Heatmap of top DEGs
top_genes <- rownames(deg_results)[1:20]  # Select top 20 genes
DoHeatmap(combined_obj, features = top_genes) + scale_fill_viridis_c()

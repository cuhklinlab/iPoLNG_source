library(Seurat)
library(Matrix)

inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")
RNA <- inputdata.10x$`Gene Expression`
ATAC <- inputdata.10x$Peaks
RNAsubset <- RNA[,colnames(ATAC)]
RNA_seurat <- CreateSeuratObject(counts = RNAsubset, project = "RNA", min.cells = 10, min.features = 500)
RNA_seurat <- NormalizeData(RNA_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
RNA_seurat <- FindVariableFeatures(RNA_seurat, selection.method = "vst", nfeatures = 5000)
top5k <- VariableFeatures(RNA_seurat)
ATAC_seurat <- CreateSeuratObject(counts = ATAC, project = "ATAC", min.cells = 5, min.features = 200)
ATAC_seurat <- NormalizeData(ATAC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
ATAC_seurat <- FindVariableFeatures(ATAC_seurat, selection.method = "vst", nfeatures = 20000)
top20k <- VariableFeatures(ATAC_seurat)
RNA_subset <- RNA_seurat@assays[[1]]@counts[top5k,]
ATAC_subset <- ATAC_seurat@assays[[1]]@counts[top20k,]
cellnames <- intersect(colnames(RNA_subset),colnames(ATAC_subset))
RNA_subset <- RNA_subset[,cellnames]
ATAC_subset <- ATAC_subset[,cellnames]
RNA_subset <- RNA_subset[rowSums(RNA_subset)>0,]
ATAC_subset <- ATAC_subset[rowSums(ATAC_subset)>0,]

# save data
writeMM(RNA_subset,file="10xPBMC3k_RNA5k.mtx")
write.csv(rownames(RNA_subset),file="10xPBMC3k_RNA5kgenes.csv",quote=F,row.names=F)
writeMM(ATAC_subset,file="10xPBMC3k_DNA20k.mtx")
write.csv(rownames(ATAC_subset),file="10xPBMC3k_DNA20kbins.csv",quote=F,row.names=F)
write.csv(colnames(RNA_subset),file=file.path("10xPBMC3k_barcodes.csv"),quote=F,row.names=F)

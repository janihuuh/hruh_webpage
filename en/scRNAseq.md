# An example of scRNAseq analysis done in the lab

28th March 2023; added by Jani Huuhtanen (jani.huuhtanen@helsinki.fi); based on the BCL-XL manuscript

## Intro to data

Scripts to reproduce figures and analyses in the manuscript "Erythroid/megakaryocytic differentiation confers BCL-XL dependency and venetoclax resistance in acute myeloid leukemia" Kuusanmäki, Dufva et al. Blood 2022

### To reproduce the results, obtain source data from Synapse:

- Get synapse credentials https://www.synapse.org
- Access synapse project syn24200411 (https://www.synapse.org/bclxl_aml)
- Download project data:
	- Input files individually (see scripts for filenames and download from SYNAPSE) (Recommended) 
	- Programmatic access (synapse, check synapse IDs from synapse):
		```
		pip install synapseclient
		synapse get synapseID
		```
	- Synapse bulk (SIZE):
		```
		pip install synapseclient
		synapse get syn24200411 -r
		```
		
		
### Processed scRNA-seq data is also available at ArrayExpress
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12607


## The index patient (AML1) samples

```{R}
folders        <- list.dirs("data/scRNAseq/AML1/", recursive = T)[-c(1,17)]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
bcl_seurat     <- scrnaseq_files[[1]]
```

## Basic QC
```{R}
getQC <- function(seurat_object){

  ###################

  min_mito     <- 0
  max_mito     <- 15

  min_ribo     <- 5
  max_ribo     <- 50

  min_features <- 300
  max_features <- 5e3

  min_counts   <- 1e3
  max_counts   <- 30e3


  ###################

  seurat_object@meta.data$barcode <- colnames(seurat_object)

  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()

  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()

  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)

  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))

  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))

  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)

}
```

```{R}
bcl_seurat <- bcl_seurat %>% getQC()
```

## Get SingleR predictions; omit predictions from cell types rare than 10 cells

```{R}
getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){

  hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
  blueprint <- SingleR::BlueprintEncodeData()
  ## @ params
  ## cluster = possible cluster vec, if not provided, tries to find in meta.data$cluster
  ## method = if "cluster", then performs preds based on clusters, not cells
  ## sample = to subsample or not

  if(!is.null(sample)){

    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])

  }

  sce       <- as.SingleCellExperiment(seurat_object)

  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)

    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.blu$first.labels
      return(seurat_object)
    }

    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }

  }


  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}
```

```{R}
bcl_seurat             <- bcl_seurat %>% getSingler()
relevant_hpca_clusters <- bcl_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- bcl_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_blueprint_pred)

bcl_seurat$singler_hpca_pred      <- ifelse(bcl_seurat$singler_hpca_pred %in% relevant_hpca_clusters, bcl_seurat$singler_hpca_pred, "rare")
bcl_seurat$singler_blueprint_pred <- ifelse(bcl_seurat$singler_blueprint_pred %in% relevant_blue_clusters, bcl_seurat$singler_blueprint_pred, "rare")
```

## Get doublets

```{R}
getDoublets <- function(seurat_object){

  require(scds)

  # Annotate doublet using co-expression based doublet scoring:
  sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = seurat_object@assays$RNA@counts), colData = seurat_object@meta.data)
  sce_object <- cxds(sce_object)
  sce_object <- bcds(sce_object)
  sce_object <- cxds_bcds_hybrid(sce_object)

  ## Add into Seurat
  seurat_object$cxds_doublet_score   <- SingleCellExperiment::colData(sce_object)$cxds_score
  seurat_object$bcds_doublet_score   <- SingleCellExperiment::colData(sce_object)$bcds_score

  seurat_object$hybrid_doublet_score <- SingleCellExperiment::colData(sce_object)$hybrid_score
  seurat_object$cxds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$cxds_score - min(SingleCellExperiment::colData(sce_object)$cxds_score)) / max(SingleCellExperiment::colData(sce_object)$cxds_score)
  seurat_object$bcds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$bcds_score - min(SingleCellExperiment::colData(sce_object)$bcds_score)) / max(SingleCellExperiment::colData(sce_object)$bcds_score)
  return(seurat_object)

}

```

```{R}
bcl_seurat <- bcl_seurat %>% getDoublets()
bcl_seurat <- subset(bcl_seurat, hybrid_doublet_score < 1.8)
```

## Get Seurat

```{R}
preprocessSeurat <- function(orig_object, cells.to.use){

  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)

  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta

  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)

  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]

  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))

  ## Scale data
  object <- ScaleData(object, features = hvg)

  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))

  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)

  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df

  return(object)

}


getClonalityGenes <- function(object){

  clonality_genes <- c(grep("^TRAV", rownames(object), value = T), grep("^TRBV", rownames(object), value = T),
                       grep("^TRGV", rownames(object), value = T), grep("^TRDV", rownames(object), value = T),
                       grep("^IGLV", rownames(object), value = T), grep("^IGLC", rownames(object), value = T),
                       grep("^IGLL", rownames(object), value = T), grep("^IGKV", rownames(object), value = T),
                       grep("^IGHV", rownames(object), value = T), grep("^IGKC", rownames(object), value = T),
                       grep("^IGH", rownames(object), value = T),  grep("^IGK", rownames(object), value = T))

}



```{R}
getDoublets <- function(seurat_object){

  require(scds)

  # Annotate doublet using co-expression based doublet scoring:
  sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = seurat_object@assays$RNA@counts), colData = seurat_object@meta.data)
  sce_object <- cxds(sce_object)
  sce_object <- bcds(sce_object)
  sce_object <- cxds_bcds_hybrid(sce_object)

  ## Add into Seurat
  seurat_object$cxds_doublet_score   <- SingleCellExperiment::colData(sce_object)$cxds_score
  seurat_object$bcds_doublet_score   <- SingleCellExperiment::colData(sce_object)$bcds_score

  seurat_object$hybrid_doublet_score <- SingleCellExperiment::colData(sce_object)$hybrid_score
  seurat_object$cxds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$cxds_score - min(SingleCellExperiment::colData(sce_object)$cxds_score)) / max(SingleCellExperiment::colData(sce_object)$cxds_score)
  seurat_object$bcds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$bcds_score - min(SingleCellExperiment::colData(sce_object)$bcds_score)) / max(SingleCellExperiment::colData(sce_object)$bcds_score)
  return(seurat_object)

}

getClustering <- function(seurat_object){

  ## Clustering
  res           <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = c(1:ncol(seurat_object@reductions$pca@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}


```

```{R}
clonality_genes <- getClonalityGenes(bcl_seurat)

bcl_seurat <- bcl_seurat %>% preprocessSeurat(cells.to.use = colnames(bcl_seurat))
bcl_seurat <- bcl_seurat %>% getClustering()
```

## Make scVI input

```
getLatentClustering <- function(seurat_object){

  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}

putLatentsSeurat <- function(seurat_object, latent){

  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)

  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)

  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))

  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}

seuratToScvi <- function(seurat_object, file = “results/seurat_object.h5Seurat”){
  
  idents.to.keep <- seurat_object@meta.data %>% group_by(orig.ident) %>% summarise(n=n()) %>% filter(n>3) %>% pull(orig.ident)
  cells.to.keep  <- seurat_object@meta.data %>% filter(orig.ident %in% idents.to.keep) %>% pull(barcode)
  seurat_object  <- subset(seurat_object, cells = cells.to.keep)
  seurat_object_diet <- DietSeurat(seurat_object)
  seurat_object_diet@assays$RNA@data <- seurat_object_diet@assays$RNA@counts
  SeuratDisk::SaveH5Seurat(seurat_object_diet, filename = file)
  SeuratDisk::Convert(file, dest = “h5ad”)
  
}
```

```{R}
dir.create("results/scvi/", showWarnings = F)
dir.create("results/scvi/input_files/", showWarnings = F)
bcl_seurat$orig.ident <- gsub("\\/", "\\_", bcl_seurat$orig.ident)
bcl_seurat %>% getScviInput(folder = "results/scvi/input_files/")
```

## Get scVI results
```{R}
latents    <- fread("results/scvi/results/bcl_latent.csv")
bcl_seurat <- bcl_seurat %>% putLatentsSeurat(latent = latents)
bcl_seurat <- bcl_seurat %>% getLatentClustering() %>% fixSeurat()
```

## Decide on clustering


```{R}
plotClustering <- function(seurat_object){

  res <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]

  q <- NULL; i <- 1

  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }

  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()

}
```


```{R}
bcl_seurat %>% plotClustering()
Idents(bcl_seurat) <- bcl_seurat$RNA_snn_res.0.7 %>% as.character() %>% as.numeric() %>% as.factor()
```

## Get DEGs
```{R}
all_markers <- FindAllMarkers(bcl_seurat, test = "t")
fwrite(all_markers, "results/de.txt", sep = "\t", quote = F, row.names = F)
```

## Save object
saveRDS(bcl_seurat, "results/bcl_seurat.rds")


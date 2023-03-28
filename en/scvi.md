## scVI analysis

28th March 2023; added by Jani Huuhtanen (jani.huuhtanen@helsinki.fi)

## scVI intro

scVI is a Bayesian swiss knife for doing single-cell analysis. In general, I mostly use scVI to gain a _latent_ _representation_ of the scRNAseq count matrices, where different nuisance factors (i.e., batch) can be considered. Essentially the latent is then used in a "basic" Seurat workflow instead of the PC:s, i.e. I count the UMAP from latents and clustering from latents.

Articles that can be of help:

* [Original publication](https://www.nature.com/articles/s41592-018-0229-2)
* [Original editorial](https://www.nature.com/articles/s41592-018-0230-9)
* [Newer publication](https://www.nature.com/articles/s41587-021-01206-w)
* [scVI Github](https://github.com/scverse/scvi-tools)

For (computational) biologists, the best article is probrably [this one](https://www.embopress.org/doi/full/10.15252/msb.20199198). With this article, try to answer the following questions:

1) What is a generative model
2) What is a latent representation
3) What are the pros and cons of PCA; what are the pros and cons of latents

After this assignment, you are ready to move on 

## scVI pragmatically 

The way I use scVI, is in Aalto Universitys Triton supercomputing platform, but also for smaller data sets (tested at ~200 000 cells) local computation is also feasible. There are options for Python (optimal) and for R (suboptimal, but workable), with GPU as well. 

Prgmatically, I use scVI as follows:

1) from Seurat to scVI
2) calculate latents with scVI, where each sample is a batch (note, if this doesn't work, you can add other nuiscance factors to the model or regress e.g., cell cycle)
3) put the latents to Seurat-object
4) Do UMAP from the latents, 
5) Do clustering from the latents

## 1) from Seurat to scVI:

```{python}
seuratToScvi <- function(seurat_object, file = "results/seurat_object.h5Seurat"){
  
  idents.to.keep <- seurat_object@meta.data %>% group_by(orig.ident) %>% summarise(n=n()) %>% filter(n>3) %>% pull(orig.ident)
  cells.to.keep  <- seurat_object@meta.data %>% filter(orig.ident %in% idents.to.keep) %>% pull(barcode)
  seurat_object  <- subset(seurat_object, cells = cells.to.keep)
  seurat_object_diet <- DietSeurat(seurat_object)
  seurat_object_diet@assays$RNA@data <- seurat_object_diet@assays$RNA@counts
  SeuratDisk::SaveH5Seurat(seurat_object_diet, filename = file)
  SeuratDisk::Convert(file, dest = "h5ad")
  
}

## Then use the function
blastpre_seurat %>% seuratToScvi(file = "results/blastpre_seurat.h5Seurat")

```

## 2) calculate latents with scVI, where each sample is a batch

Note that this is done in Python in the Aalto cluster Triton; setting up the scVI is a different topic

```{python}
import scvi; import scanpy as sc; import os; import numpy as np
os.chdir("/scratch/cs/csb/projects/bclxl/")

## Load and preprocess
seurat_data = scvi.data.read_h5ad("results/scvi/input_files/hemap_blast_pre/blastpre_seurat.h5ad")
seurat_data.layers["counts"] = seurat_data.X.copy() # preserve counts
sc.pp.normalize_total(seurat_data, target_sum=1e4)
sc.pp.log1p(seurat_data)
seurat_data.raw = seurat_data # freeze the state in `.raw`

## Train
scvi.data.setup_anndata(seurat_data, layer = "counts", batch_key = "orig.ident")
seurat_data_model=scvi.model.SCVI(seurat_data)
seurat_data_model.train()

## Save it
seurat_data_model.save("results/scvi/output/blast_seurat_diet/")
latent = seurat_data_model.get_latent_representation()
np.savetxt("results/scvi/output/hemap_blast_pre_latent.csv", latent, delimiter=",")
```

## 3) put the latents to Seurat-object

This is back in R.

```{R}
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

latents <- fread("results/scvi/output/hemap_blast_pre_latent.csv")
blastpre_seurat <- blastpre_seurat %>% putLatentsSeurat(latents)

```


## 4) Get latent UMAP

```{R}
getLatentUMAP <- function(seurat_object){

  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df

  return(seurat_object)

}

blastpre_seurat <- blastpre_seurat %>% getLatentUMAP()


```



5) Get latent clustering

```{R}
getLatentClustering <- function(seurat_object){

  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}

blastpre_seurat <- blastpre_seurat %>% getLatentClustering()


```

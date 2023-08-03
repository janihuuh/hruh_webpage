## Running scVI in Atlas

Added 3 August 2023 by Aino Häkkinen (aino-elina.hakkinen@helsinki.fi)

Refer to [scVI intro](scvi.md) by Jani for scVI overview. This is the workflow I've used for running scVI in Atlas, which is somewhat faster than running it locally. Heavily copied from Jason Theodoropoulos. 

### Set up conda environment 
Clone Jason's conda environment with scvi
```python
conda create --clone /csc/mustjoki/anaconda3/envs/scvi --name my-scvi
```

### Convert Seurat object to h5ad locally
Run the following in R to export your Seurat object in h5ad format
```R
library(Seurat)
library(SeuratDisk)

# copied from https://zqfang.github.io/2020-04-28-seurat2scanpy/

top2000 <- head(VariableFeatures(seurat_object), 2000)
top2000 <- seurat_object[top2000]

# slim down a Seurat object to get raw counts and lognorm counts
seu = DietSeurat(
  top2000,
  counts = TRUE, 
  data = TRUE,   
  scale.data = FALSE, 
  features = rownames(top2000),
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), 
  misc = TRUE
)

i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

SaveH5Seurat(seu, filename = "my_data.h5seurat", overwrite = TRUE)
Convert("my_data.h5seurat", "my_data.h5ad", assay="RNA", overwrite = TRUE)

```

### Upload h5ad file to cluster
```
scp /path/to/my_data.h5ad username@atlas.genome.helsinki.fi:/projects/fimm_ngs_mustjoki/path/to/scvi/input/
```

### In cluster, start a screen session
Always use a screen to make sure you don't lose work if connection is lost; read more on screens in [Atlas basics](atlas.md)
```
screen -S scvi
```
### Connect to a compute node
```
qstat -f | grep seq_hugemem2 # list all available nodes
ssh seq_compute*-*
```

### Activate conda environment
```python
source /csc/mustjoki/anaconda3/bin/activate my-scvi
```

### Start python and import libraries
```python
python
import scvi
import scanpy as sc
import numpy as np
```

### Read in h5ad file, set up model and train
```python
adata_obj = scvi.data.read_h5ad("/projects/fimm_ngs_mustjoki/path/to/scvi/input/my_data.h5ad")
adata_obj.layers['counts'] = adata_obj.raw.X.copy()
scvi.data.setup_anndata(adata_obj, batch_key="orig.ident", layer="counts")
model = scvi.model.SCVI(adata_obj, n_latent=20) # modify latents here
scvi.data.view_anndata_setup(model.adata) # optional, view model
model.train()
```

### Save latents
```python
latent = model.get_latent_representation()
np.savetxt("/projects/fimm_ngs_mustjoki/path/to/scvi/output/my_data_latents.csv", latent, delimiter=",")
model.save("/projects/fimm_ngs_mustjoki/path/to/scvi/output/my_data_model")
```

### Exit python, deactivate conda environment and kill screen session
```python
exit()
conda deactivate
```
Kill screen session with ctrl a + k

### Download latents from cluster
```
scp username@atlas.genome.helsinki.fi:/projects/fimm_ngs_mustjoki/path/to/scvi/output/my_data_latents.csv /path/to/local/scvi/output
```

### Put latents back into Seurat object
```R
latent <- read.table("/path/to/local/scvi/output/my_data_latents.csv", sep=",")
latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_object)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_object))

seurat_object <- RunUMAP(seurat_object, dims = 1:15, reduction = "scvi", reduction.name = "latentumap")
```

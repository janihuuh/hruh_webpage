## Running streamlined data integration using Seurat v5 in Atlas

Added 5 October 2023 by Johannes Smolander (johannes.smolander@helsinki.fi)

This example shows you how to run the streamlined data integration vignette of Seurat v5 in Atlas.
https://satijalab.org/seurat/articles/seurat5_integration . This is useful for testing different data integration methods.
The vignette covers 5 different integration methods (CCA, RPCA, scVI, FastMNN, Harmony)

### Prepare the input data for data integration 
Perform this step on your personal computer.

**If you have an older Seurat object (v4)**
```R
library(Seurat)

obj <- readRDS("my_seurat_v4_object.rds")

# Get the raw count data and the meta data. 
# You must have the batch label included in the metadata as a column.
# Specify the batch label in the R script (read the steps below) or define a new column.

obj$batch <- obj$patient

raw_data <- obj@assays$RNA@counts
meta_data <- obj@meta.data

save(raw_data,meta_data,file = "raw_data_and_meta_data.RData")
```

**If you have a Seurat v5 object**
```R
library(Seurat)

obj <- readRDS("my_seurat_v5_object.rds")

# Get the raw count data and the meta data. 
# You must have the batch label included in the metadata as a column.
# Specify the batch label in the R script (read the steps below) or define a new column.

obj$batch <- obj$patient

obj <- JoinLayers(obj)

raw_data <- obj[["RNA"]]$counts
meta_data <- obj@meta.data

save(raw_data,meta_data,file = "raw_data_and_meta_data.RData")
```

### Start a screen session 
It's best to start a screen session so that the user can log out from Atlas at any time during the analysis.
```bash
screen -S seurat_v5
```
Execute Ctrl+A+D to leave the session.

To return to the session, execute:
```bash
screen -r seurat_v5
```

### Connect to a computing node 

Connect to a computing node interactively. This is required for using Singularity.
```bash
qlogin -q interactive.q
```

### Go to your project folder

```bash
cd /csc/mustjoki/my_project
```

### Prepare the R script that you want to run in Seurat v5

For demonstration, we copy the example R script from /csc/mustjoki/singularity
```bash
cp /csc/mustjoki/singularity/run_streamlined_integration_seurat_v5.R .
```


### Run the singularity container

This starts the singularity container with Seurat v5 installed. 
The command binds the project folder ${PWD} to /wrkdir.
Then you can access the project folder under /wrkdir when the container is running.
```bash
singularity run --bind ${PWD}:/wrkdir /csc/mustjoki/singularity/ubuntu_seuratv5.sif
```

### Start the r-seurat-v5 virtual environment in Conda

Execute these commands to star the r-seurat-v5 virtual environment in Conda and to change directory to /wrkdir.
```bash
__conda_setup="$('/home/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate r-seurat-v5
cd /wrkdir/
```

You should now see the content of your project directory when you list the files.
```bash
ls -ltrha
```

### Run the R script inside the container. 

The example runs the scVI data integration method (https://satijalab.org/seurat/articles/seurat5_integration).

To run a different integration method, comment over the part that runs scVI and uncomment a different part that runs e.g. Harmony.

**Remember to also edit the variable names, i.e., "integrated.scvi", "umap.scvi", "scvi_clusters" in the Rscript.** 

For example, make these edits if you use Harmony instead of scVI:

- integrated.scvi -> harmony
- umap.scvi -> umap.harmony
- scvi_clusters -> harmony_clusters

As input, the script assumes an RData file (*raw_data_and_meta_data.RData*) with "raw_data" and "meta_data" objects.
Please see the beginning of this tutorial to learn how to extract them from a Seurat object.
The script also saves the results.

```bash
Rscript run_streamlined_integration_seurat_v5.R
```
While the script is running, you can leave at any time the screen session by executing "Ctrl+A+D" and come back by executing.

```bash
screen -r seurat_v5
```


### Exit the container, the computing node and the screen 

When you are done, exit the container.

```bash
exit
```

When you are done, exit the computing node.

```bash
exit
```

When you are done, terminate the screen.

```bash
exit
```

### Add the integration results to your Seurat object

Perform this step on your personal computer.
**If you have an older Seurat object (v4)**
```R
library(Seurat)

# This is the old object. Not the new one.
obj <- readRDS("my_seurat_v4_object.rds")

# This came as a result from the analysis.
cell_embeddings_scvi <- readRDS("seurat_v5_latents.rds")
#cell_embeddings_umap <- readRDS("seurat_v5_umap.rds") # If you want the umap embeddings as well

# Add the latents and UMAP as new slots to obj@reductions
obj[["integrated.scvi"]] <- CreateDimReducObject(embeddings = cell_embeddings_scvi) # If you want the umap embeddings as well
#obj[["umap.scvi"]] <- CreateDimReducObject(embeddings = cell_embeddings_umap)

# Perform clustering for the latent embeddings
# IMPORTANT: if you changed the number of dimensions in the integration step (R script), match these dimensions with them
# For example, if you used ndims=20 in scVI, set dims = 1:20 for FindNeighbors and RunUMAP.
obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "scvi_clusters")  
# Create a new UMAP embedding (you can use the one created by Seurat v5 too and just skip this, but make sure the dims was set correctly)
obj <- RunUMAP(obj,dims = 1:30,reduction.name = "umap.scvi",reduction = "integrated.scvi")


# Visualize the new UMAP embedding
DimPlot(obj,reduction="umap.scvi",group.by="scvi_clusters")
```

**Remember that if you use a different integration method instead of scVI, it's best to use different names for integrated.scvi and umap.scvi**.
**For example, if you use CCA, name them instead integrated.cca and umap.cca**.



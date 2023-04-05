# Velocyto

I have done cell velocity analyses with

1. [velocyto](http://velocyto.org/) (to get the spliced/unspliced reads)
2. [scVelo](https://github.com/theislab/scvelo) to get the downstream analyses

## Velocyto in Atlas:

This follows the basic [velocyto workflow](http://velocyto.org/)

### 1) Run samtools to get a sorted genotype. 

Notice that the INDIR is the CellRanger output, containing the outs/ folders etc.

``` 
#!/bin/bash
export PATH=/fs/vault/pipelines/t_receptor/bin/samtools-1.9/:$PATH
INDIR=/csc/mustjoki2/velocyto/data/nk_heme/Timepoints_plate1/outs/
samtools sort -t CB -O BAM -m 15000M -o $INDIR/cellsorted_possorted_genome_bam.bam $INDIR/possorted_genome_bam.bam

``` 

### 2) Run velocyto to get the unspliced / spliced reads.

```
#!/bin/bash
source /csc/mustjoki/anaconda3/bin/activate velocyto_envi
INDIR=/csc/mustjoki2/velocyto/data/nk_heme/Timepoints_plate1/
velocyto run10x $INDIR /csc/mustjoki2/velocyto/data/refdata-gex-GRCh38-2020-A/genes/genes.gtf
```

### 3) Download the .loom files and start the downstream analysis

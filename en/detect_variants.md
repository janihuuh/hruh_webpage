How to detect variants from scRNA-seq data with CellSNPLite

Jani Huuhtanen 3rd August 2023

This is based on CellSNPLite, which is faster than e.g., Vartrix. Try to familiarize yourself with the algorithms and the parameters.

1. Install envi:

```
#!/bin/bash
export PATH=$PATH:/csc/mustjoki/anaconda3/bin/

cd /csc/mustjoki/cellsnplite/
conda create --name your_name_here python=3.6

conda activate your_name_here
conda config --add channels bioconda
conda config --add channels conda-forge

conda install cellsnp-lite

conda deactivate
```

2. Use cellsnp-lite

```
#!/bin/bash
export PATH=$PATH:/csc/mustjoki/anaconda3/bin/
source /csc/mustjoki/anaconda3/bin/activate cellsnplite_test

## 1) cellsnp-lite
folder=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/ICAN_NK_AML/Batch1_120122/count_220217_A00464_0452_BHW5H3DRXY/NK_MOLM13/outs

BAM=$folder/possorted_genome_bam.bam
BARCODE=$folder/filtered_feature_bc_matrix/barcodes.tsv.gz
OUT_DIR=/csc/mustjoki/cellsnplite/results/nk_molm13_filt/
REGION_VCF=/csc/mustjoki/cellsnplite/apps/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz
mkdir $OUT_DIR
cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -R $REGION_VCF -p 5 --minMAF 0.1 --minCOUNT 20 --gzip
```

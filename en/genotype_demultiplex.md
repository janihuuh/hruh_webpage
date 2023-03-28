# Genotype based demultiplexing on scRNAseq data

Added 28th March 2023 by Jani Huuhtanen (jani.huuhtanen@helsinki.fi)

The workflow consits of three parts, 

1) Gather the SNPs from scRNAseq. This could be based on common SNPs or if you have some orthogonal information, then you can use them as well.
2) (Optional) Calculate the mitochondrial variants (based on mquad; other solutions exist)
3) Do "clustering" based on the SNPs with vireo

This is highly based on the [vireo-workflow](https://vireosnp.readthedocs.io/en/latest/) and the [MQuad-workflow](https://github.com/single-cell-genetics/MQuad) 


```
{
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

## 2) mquad for informative mito variants
INPUT_DIR=/csc/mustjoki/cellsnplite/results/nk_molm13_filt/
OUT_DIR=/csc/mustjoki/cellsnplite/results/nk_molm13_filt/mquad/
mkdir $OUT_DIR
mquad -c $INPUT_DIR -o $OUT_DIR -p 5

## 3) vireo
INPUT_DIR=/csc/mustjoki/cellsnplite/results/nk_molm13_filt/
OUT_DIR=/csc/mustjoki/cellsnplite/results/nk_molm13_filt/vireo/
n_donor=4
mkdir $OUT_DIR
vireo -c $INPUT_DIR -N $n_donor -o $OUT_DIR

}
```

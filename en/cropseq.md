Added April 17 2023 by Olli


# CROP-seq analysis

## 1.	Get sgRNA barcode UMI/read counts
⁃	FASTQ files of the targeted sgRNA amplification libraries are run through Cell Ranger count pipeline by the core
⁃	Extract UMI counts of guides associated with each cell using the get_barcodes.py script (https://github.com/shendurelab/single-cell-ko-screens, Hill et al., 2018).
⁃	Analysis scripts are in the following folder and result files can be stored here as well 

```bash
/csc/mustjoki2/scrna_seq/sc-ko-screen
```

-	Look at sudhl4 directory as an example
⁃	Make a new directory for your datasets to sc_ko_screen and in that directory do the following:
  ⁃	Make whitelist.txt which has a list of all used sgRNA sequences 
  ⁃	Make a new bash script for every sample as below (you can copy run_ctrl.sh and change the following):
      o	Change —input_bams to the path to your BAM file (use the directory with “_CR” in the name which indicates the targeted sgRNA site amplification)
      o	Change output filename (-o) to your preference e.g. “experimentname_cropseq_cr.txt”
  
  ```bash
newgrp sg_mustjoki;
export PYTHONPATH="";
export PYTHONPATH=/apps/python3.7.1/lib/python3.7/site-packages:$PYTHONPATH;
export PYTHONPATH=/apps/python3.7.1/Lib:$PYTHONPATH;
export PYTHONPATH=/fs/vault/pipelines/t_receptor/bin/sc-ko-screen-shendure/include/lib/python3.7/site-packages/:$PYTHONPATH;
/apps/python3.7.1/python /fs/vault/pipelines/t_receptor/bin/sc-ko-screen-shendure/get_barcodes.py --input_bams /fas/NGS/pipes/cellranger/fimm_sca_dufva/CROPseq_NK_SUDHL4_MM1S/Batch1-2_131120-271120/count_210301_A00464_0291_BHYG2NDSXY/SUDHL4_CROPseq/outs/possorted_genome_bam.bam -o sudhl4_cropseq.txt --whitelist whitelist.txt --search_seq GTGGAAAGGACGAAACACCG --all_reads --force_correction 2
```
      
⁃	If you want, you can also run the command on the transcriptome library (directory without “_CR”) as in the original CROP-seq paper but much fewer reads mapping to the sgRNAs will be found
⁃	Run the bash scripts in the queue


## 2.	Analysis in R

Create object

Move features, barcodes, matrix files from the “outs” directory and the barcode counting output (“experimentname_cropseq_cr.txt”) to your computer

Script example to generate Seurat object with sgRNA-cell pairing in /csc/mustjoki2/scrna_seq/sc-ko-screen/r_script_examples/run_makeObjectSUDHL4_singlet.R

⁃	Briefly, what the script does is the following:
⁃	Merge data from different samples (10x lanes) e.g. treated and untreated conditions (this may be useful to be able to compare treatments)
⁃	Basic scRNA-seq QC (check that thresholds for QC parameters look good for your data and adjust as needed)
⁃	Map the cells to sgRNAs and select cells with single sgRNA detected
⁃	To assign guides to cells, we have included in the analysis cells harboring sequences with > 10 UMI counts and accounting for > 50% of the UMI counts in the cell. Out of these, cells in which the second most frequent guide accounted for > 20% of the UMI counts were considered to express two guides and were removed from the analysis
⁃	Cluster cells and remove unwanted (e.g. NK cells/other cells not harboring any sgRNAs) or low quality clusters (this part needs you to select the wanted clusters by looking at differentially expressed genes/singleR annotations etc)
⁃	The resulting object should have only good quality cells expressing a single sgRNA


Run differential expression analyses

⁃	Script example to perform differential expression analysis using the Seurat object generated above: /csc/mustjoki2/scrna_seq/sc-ko-screen/r_script_examples/run_analyse_SUDHL4_singlet.R 
⁃	Briefly, what the script does is the following:
  ⁃	Compute DEG between each perturbation and control (you need to change condition names from e.g. “crop_seurat_nk_1_16” and “results_nk_1_16” to your conditions)
  ⁃	Plot volcano plots of DEG between each perturbation and control and dot plots of top DEG for each perturbation
  ⁃	Compute DEG and plot volcano plots between treatments using cells with control sgRNAs from each condition
  ⁃	Calculate enrichment/depletion of cells carrying different sgRNA between treatment conditions (using different sgRNAs as biological replicates)

Run mixscape analyses (Papalexi et al. Nature Genetics 2021)

⁃	Script example to perform differential expression analysis using the Seurat object /csc/mustjoki2/scrna_seq/sc-ko-screen/r_script_examples/run_mixscape_SUDHL4_singlet.R 
⁃	Script adapted from https://satijalab.org/seurat/articles/mixscape_vignette.html


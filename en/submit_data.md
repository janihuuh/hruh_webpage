# Submitting and publishing data

Added Jani Huuhtanen 5th April 2023

## TCR&beta; sequencing done at Adaptive 

This is quite easily done via them, and they are quite fast as well. You just need to contact techsupport@adaptivebiotech.com with the following information:

Spreadsheet format:

* name that the sample is in the Adaptive system immuneAccess
* any metadata for the sample

Text format:

* Title
* Abstract
* Primary investigator
* Authors
* Affiliations


## Where to submit single-cell RNA-sequencing data

The more the merrier, eh? Multiple places to upload the data. I think the University of Helsinki lawyers interpret the GDPR regulations as such that we cannot:

* Upload data to servers that are in the USA (ie., GEO, which the biggest database)

Thus, at least I have used:

* Raw data = [EGA](https://www.google.com/search?q=ega&oq=ega&aqs=chrome.0.69i59j35i39j0i67i650j0i512j46i175i199i512j0i512l2j69i60.335j0j7&sourceid=chrome&ie=UTF-8)
* Count matrices, TCRab-seq = [ArrayExpress](https://www.ebi.ac.uk/fg/annotare/)
* Seurat-objects, TCRab-seq, any other data = e.g., [Zenodo] 

For example in out [T-LGLL paper](https://www.nature.com/articles/s41467-022-29173-z#data-availability), we have used following text in data availability:

> The processed scRNA-sequencing data for both the T-LGLL and healthy samples generated in this study are available at ArrayExpression under accession code E-MTAB-11170. The raw scRNA-sequencing and bulk-RNA-sequencing are available in the European Genome-Phenome Archive under accession code EGAS00001005297. The TCRαβ-sequencing data, TCRβ-sequencing data, and Seurat-objects are available at Zenodo under: https://doi.org/10.5281/zenodo.4739231 [https://zenodo.org/record/4739231] with restricted access due to GDPR regulations and data can be accessed by placing a request via Zenodo. The publicly available scRNA+TCRαβ-sequencing and TCRβ-sequencing data used in this study are listed in Supplementary Data 1. Source data are provided with this manuscript. Source data are provided with this paper.

## Uploading data to EGA

Here you upload your .fastq/.bam files, i.e. the raw sequencing files. You need to upload the raw data as well if you are planning to upload the count matrices - then you need to link EGA and Annotare together.

Jason Theodoropoulos / Olli Dufva knows this

## Uploading data to Annotare

### 1. What to upload.

Here you upload the count matrices, i.e. the output from the CellRanger, i.e.:

* barcodes.tsv.gz
* features.tsv.gz
* matrix.mtx.gz

Preferably in this format, the .gz-compressed files are easier to transform.

When assigning data to ther matrix or barcode files should be assigned as the _processed_ _data_ _files_, and not e.g., raw data files (that are .fastq / .bam)

### 2. Follow the instructions

Essentially, you need to add here your meta data corresponding to samples and the methods of how you processed the samples. 

If you need help, consider checking what has already been included, e.g. in [here](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11170).

### 3. Link ArrayExpress to EGA

Send an email to Annotare saying that yoy have data now in Annotare and EGA and you need to link them. E.g., check here:

> to: annotare@ebi.ac.uk:
> 
> “Hi,

> We (Jason Theodoropoulos cc’d and I) are submitting data to ArrayExpress with annotare with my account (jani.huuhtanen@helsinki.fi) by the title “Single-cell profiling of the leukemic and non-leukemic immune cell compartments in CD8+ T-cell Large Granular Lymphocytic Leukemia”. All the processed data has been uploaded to Annotare but our raw data is in EGA with the study ID EGAS00001005297.
We would like to combine the EGA submission with the ArrayExpress submission in such way that the raw data would be in EGA with restricted access but processed data in ArrayExpress. The subjects Pt1 to Pt9 and HC1 to HC6 correspond to Patient1 to Patient9 and HC1 to HC6 in the EGA submission.

> Could you help us in this? Please let us know if you need any additional information.

> Thanks in advance,
> Jani Huuhtanen”

## Uploading data to Zenodo

Quite straightforward if you follow the instrictuions. Requires the least formatting. 

https://zenodo.org/


## Uploading code

Preferebly to [GitHub](https://github.com/), with accompanying DOI that you can cite. 

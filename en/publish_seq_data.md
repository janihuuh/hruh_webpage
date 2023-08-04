# Submitting and publishing scRNAseq data

Added Sofie Lundgren / August 3rd 2023

## What should be published and where? Who should have access to the data?

For every sample, there are raw data (.fastq-files) and processed data (*counts*, *matrix* and *barcodes.tsv* -files for gene expression and *clonotypes.csv*, *consensus\_annotations.csv* and *filtered\_consensus\_annotations.csv* -files for V(D)J expression). Many publishers request that the raw data will be published as well, but it is not ethical to allow unlimited download of the raw genomic sequencing data, as it is sensitive data. However, most researchers are using same pipelines (e.g. CellRanger) for raw data processing, and hence publishing the processed data (which is not sensitive) is the best way to support open science.

One good option is to publish the raw data in European Genome Archive (EGA), where the data will be stored on European servers but won’t be available for free download without a separate permission from a named Data access committee (our group). This way you can avoid immediate data privacy concerns. Simultaneously, the processed data can be published e.g. in ArrayExpress, where it can be freely downloaded by other reseachers.

## Tips for publishing your processed data in ArrayExpress

* Follow the instructions on ArrayExpress homepage (https://www.ebi.ac.uk/training/online/courses/array-express-quick-tour/submitting-data-to-arrayexpress/) 
* The submission is made through the Annotare -portal (https://www.ebi.ac.uk/fg/annotare/) 
* You can set the publishing date to be e.g. after 1 year, and after your article has been accepted for publication, the date can changed and data released.
* In Annotare, you’ll need to describe the sequencing protocol accurately (e.g. library layout, library selection, library strategy, cDNA, cell barcode, sample barcode read lengths etc). If you are not sure which are the right options for your dataset, please confirm the parameters from single-cell sequencing unit (Jenni Lahtela / Emma Saarinen / Anna Näätänen). Most information can be found in the iLab request (Comments) or data release email.
* If you are not uploading the raw data to ArrayExpress (which you should not do), you’ll get an error message while trying to submit your dataset, asking you to assign raw data files. At this point, you should send a request to ArrayExpress Helpdesk, asking them to override the validation error and complete the submission for you. In the message, you can tell that the raw data will be published in EGA but not in ArrayExpress due to data privacy concerns. Make sure that before asking the curator to complete the submission, you have filled all other necessary information correctly in the portal.

## How to upload your raw data to EGA
1. Send a submission request to EGA helpdesk (fill form on EGA webpage https://ega-archive.org/submission-form.php), to get an ega-box address for your data upload.

2. Encrypt your files using EGAcryptor
    * Use queues on FiMM server! (See the instructions from “How to Atlas”)
    
    ```java -jar /csc/mustjoki2/hemap\_immunology/ega\_submission\_hemap/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar -i /path/to/your/files/ -o /path/to/your/new/folder/encrypted/```

        * Input (-i) can be a parent directory containing all patient/sample folders, the encrypted files will have same directory structure
    * Note! There are many fastq files for each sample (e.g. 4 fastq files gene expression and as many for V(DJ) expression), don’t get confused.
    
3.	Create metadata files
    1. To add md5 checksums to metadata table and get all filenames and paths, first make a copy of the checksum files
        ```cp --parents -R ./\*/\*.md5 ./md5```
    2.  Move the “md5” folder to your computer and use attached [R script](scrnaseq_ega_submission_metadata_3_and_5_.R), to create metadata table from the filenames
        * !!! Please do not edit the shared r script but make a copy of it to your computer which you can edit as much as you want!
        * See How to atlas instructions fro moving files (scp -r command)
        * 5'VDJ data -> use [scrnaseq_ega_submission_metadata_5_and_VDJ.R](scrnaseq_ega_submission_metadata_5_and_VDJ.R)
        * Files that you need to run the script:
            * Sample template filled with your sample meta information (example: ega\_aa\_scrna\_samples.csv, template can also be found on ega webpage)
         * (Optional: manuscript\_sample\_IDs.csv)
3.	Upload all directories to EGA (this can take a couple of days if you have a lot of files)
    * Remember to remove identifiers from file names accordingly as original file names are visible for people downloading your data from EGA. 
        * e.g. ```find . -name ‘FHRB_1234*’ -exec rename -v FHRB_1234 Patient1 {} \``` changes all prefixes with the FHRB identifier to Patient1.
    * Use a screen (see the instructions from How to Atlas)
        * No need to use queues but via the screen the upload can continue even if your connection is lost
    * ```ncftpput -u ega-box-xxxx -p password -R ftp.ega.ebi.ac.uk ./ /path/to/your/new/folder/encrypted/``` 
        * If your upload is for some reason disconnected or interrupted, you can add the argument -z to the ncftpput command to continue uploading from where the last upload was interrupted.
    * Commands suggested on EGA website can be used if all files are in the same directory but don’t work recursively with subdirectories
4.	Upload metadata files in EGA portal to link files to samples and add metadata
    * Register samples: add the metadata file (e.g. ega\_aa\_submission\_metadata.csv)
    * Create 2 experiments, one with “library layout: PAIRED” (for read files) and other with “library layout: SINGLE” (for index files), see below
    * Link files and samples: add the files containing r and i files (e.g. ega\_submission\_5prime\_files\_md5\_i.csv and ega\_submission\_5prime\_files\_md5\_r.csv)
    * Note! Metadata for the index files and read files should be downloaded separately to different experiments
    * Metadata for index files (containing “\_i” in the filename): select “your\_SINGLE\_experiment” in Step 1 and “one Fastq file (Single)” -option in Step 2
    * 	Metadata for read files (containing “\_r” in the filename): “your\_PAIRED\_experiment” in Step 1 and “Two Fastq files (Paired)” -option in Step 2
5.	Create an “analysis file” for the samples, to avoid a validation error in the portal
    * E.g. list of samples or readme.txt file containing information of the dataset
    * These need to go through the same encryptor and upload -steps as the raw data files, and then you can link them manually to the samples in the portal
    * An example of naming a sample in the portal:
 
 ##### Samples -section

|        |           |  |
| ---------- |---------------------------------| ------------------------------|
| Alias      | AA3\_10\_2015\_BM\_S32\_L001\_R | Should be unique for each file|
| Title      | AA3\_10\_2015\_BM\_S32\_L001\_R      |  |
| Subject ID | AA3\_10\_2015\_BM      | Should be unique for each sample |

 ##### Link files and samples -section

|        |           |  |
| ---------- |---------------------------------| ------------------------------|
| Alias      | 0371f9d4-c63b-46… | Created by the submission system |
| Sample Alias      | AA3\_0\_2015\_BM\_S32\_L001\_R | Should be same as in the samples -section |
| Filename | encrypted/FHRB1641\_BM\_oct15\_TCR\_4/FHRB1641\_BM\_oct15\_TCR\_L001\_R\_fastq.gz | Should match the name of uploaded file |


An example of data access policy:

 |        |           | 
| ---------- |---------------------------------| 
| Select a DAC      | Hematology Research Unit Helsinki Data Access Committee | 
| Short description      | Hematology Research Unit Helsinki Data Access Policy | 
| Policy | The proposed research use must be consistent with the specific data use limitation for this study. Data requests will be evaluated based on the Helsinki University Hospital ethics committee permit concerning the samples. | 

# TCRseq basics

28th March 2023; added by Jani Huuhtanen (jani.huuhtanen@helsinki.fi)

## Intro

Most useful pages:

* [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/input.html)
* 

## Getting the data from Adaptive immuneAccess

If data is downloaded from [immuneAccess](https://clients.adaptivebiotech.com/login), you should Exprt with v2 

## 0) Change Adaptive to vdjtools format 

Init the variables used in the analysis; change to something that is suitable to your analysis
```
me=$(whoami)
files=/Users/$me/Dropbox/aplastic_anemia_tcr

vdj=$files/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjdb=$files/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=$files/data/selected_tcrb/databases/vdjdb_new
output_path=$files/results/

cd $files/data/unselected_tcrb/Bethesda/;
bethesda=$(ls -d "$PWD"/*);

```

Then convert with the vdjtools

```
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $bethesda data/unselected_tcrb/Bethesda/vdj/
```

## 1) Filter non-functional TCRs 

```
java -Xmx4G -jar $vdj FilterNonFunctional $bethesda data/unselected_tcrb/Bethesda/filtered/;
```

## 2) Downsample (if needed)

```
java -Xmx4G -jar $vdj DownSample --size 10000 $bethesda            data/unselected_tcrb/resampled/Bethesda/
```

## 3) Estimate diversities

```
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $melbourne      $files/results/diversity/raw/melbourne
```

## 4) Query against vdjdb

```
java -Xmx4G -jar $vdjdb --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $bethesda       $output_path/raw/bethesda/
```


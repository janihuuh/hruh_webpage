# TCRseq basics

28th March 2023; added by Jani Huuhtanen (jani.huuhtanen@helsinki.fi)

## Intro

Most useful pages:

* [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/input.html)
* [vdjdb](https://vdjdb.cdr3.net/)
* [screpertoire](https://github.com/ncborcherding/scRepertoire)
* [https://immunarch.com/](https://immunarch.com/)

Most of the TCR&beta;-seq lab in the data is coming from the Adaptive platform. If you don't have access to it, you should email tiina.hannunen@helsinki.fi to get the access to it. 

## Getting the data from Adaptive immuneAccess

If data is downloaded from [immuneAccess](https://clients.adaptivebiotech.com/login), you should Exprt with v2. 

You can of course browe around immuneAccess to see if there's something relevant to your topic and you can analyze that data.

## 0) Change Adaptive to vdjtools format 

Personally, I found the [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/input.html) a really useful programme with quick turnaround time and it covers all the basic commands one would need for even more advanced TCR-seq analyses. In order to work with vdjtools, you should change the data to the vdjtools format.

This code is intended to show how to init the variables used in the analysis used in bash; change it (e.g., folders to match your data, the path to match your vdjtools) to something that is suitable to your analysis.

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

Really important part, as we are usually studying only the functional TCR-repertoire. 

```
java -Xmx4G -jar $vdj FilterNonFunctional $bethesda data/unselected_tcrb/Bethesda/filtered/;
```

## 2) Downsample (if needed)

Most (if not all!) diveristy metrics (including the clonality!) are really sensitive to the sample size, i.e., the number of TCR reads per sample. There are two options you can do TCRb-seq with Adaptive, where the "Survey" is somehting that is used by the HRUH, resulting in 10 000 - 100 000 total reads, but the "Discovery" mode used by many other labs results in 1 000 000 - 100 000 000 reads per sample. That of course creates a huge bias.

One should perhaps have a look at the smallest amount of reads per sample in the project, and perhaps then assign a threshold of what to use. You could also e.g., omit all samples below a certain threshold, because it could imply that that sample isn't informative. Some values used by Jani include 10k, 20k 40k, and 100k.

```
java -Xmx4G -jar $vdj DownSample --size 10000 $bethesda data/unselected_tcrb/resampled/Bethesda/
```

## 3) Estimate diversities

Estimate all different diversities. Remember that clonality is not the same as Gini and it's not the same as Simpson - there metrices capture different types of diversity estimates and are different and should be analyzed as such!

```
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $bethesda $files/results/diversity/raw/melbourne
```

## 4) Query against vdjdb

This part is called the "hard matching" against the [VDJdb](https://github.com/antigenomics/vdjdb-db). You should familiarize yourself with the VDJdb, what does it include, and especially the confidence score and how they have calculated it:


score | description
------|----------------------
0     | Low confidence/no information - a critical aspect of sequencing/specificity validation is missing
1     | Moderate confidence - no verification / poor TCR sequence confidence
2     | High confidence - has some specificity verification, good TCR sequence confidence
3     | Very high confidence - has extensive verification or structural data


```
java -Xmx4G -jar $vdjdb --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $bethesda $output_path/raw/bethesda/
```


# GLIPH2

Added 29th March 2023 by Jani Huuhtanen (jani.huuhtanen@helsinki.fi)

The [GLIPH2 webpage](http://50.255.35.37:8080/tools) is really useful, and you can even do the analysis there. Although I would highly recommend downloading the command line interface tool. 

The most important thing in the CLI tool is to but the TCRa info as "", and not as NA or "NA" or anything else, that will end up in failure.

## Comments on GLIPH2:

* it is deterministic, i.e., the output is dependent on the TCRs you input it. If you subselect / add TCRs, don't expect the same output
* it will always cluster ~40% of the TCRs you inputted there
* most of the motifs can be non-specific to your question, thus it is advicable that you blacklist the motifs based on some orthogonal data, e.g., blacklist the motifs found on healthy TCRb-seq data from e.g., [this](https://www.nature.com/articles/ncomms15869)
* you cannot do antigen-specificity predictions with GLIPH2! If you input some TCRs with known specificity (e.g., TCRs targeting CMV), and you see that your TCRs cluster with those, you cannot call them CMV-specific! No one has done metrics on that, ie. what is the false-positive and false-negative rates for that type of analysis. For that type of analyses you should focus on supervised analyses e.g., [TCRGP](https://github.com/emmijokinen/TCRGP)
* However, GLIPH2 can be used to pre-select TCRs for e.g., TCRGP. Then you can "clear" your e.g., tetramer-specific signal with GLIPH2 to get rid of outlier-like TCRs that could be technical artefacts (e.g., autofluoresence)

## Ideas for GLIPH2:

* GLIPH2 is really useful, but you need to be creative with the data analysis. See e.g. 

* Largest GLIPH2 analysis to date [this](https://www.cell.com/immunity/pdfExtended/S1074-7613(21)00080-7)
* Antigen-drive (see [this](https://www.nature.com/articles/s41467-022-29173-z), especially Figure 3)

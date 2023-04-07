
# Running R scripts in job queue

Running R scripts in job queue

## 1.	Make an R script e.g. 

```bash
/projects/fimm_ngs_mustjoki/tcr/analysis.R
```

## 2.	Make a bash script e.g. “run_R.sh” containing the following:

```R
/apps/statistics2/R-3.5.1/bin/R --no-save < /projects/fimm_ngs_mustjoki/tcr/analysis.R >& messages.err
```

messages.err will contain the error, warning etc. messages from R.

## 3.	Submit the bash script to queue as above.

### To install R packages, specify a path in the project directory and save it into libpath variable:

```R
package_library <- "path/to/packages"
libpath <- .libPaths(package_library)
 
source("http://bioconductor.org/biocLite.R")
biocLite("scran", lib = package_library)
biocLite("scater", lib = package_library)
```

Or from CRAN: 
```R
install.packages(“igraph”, repos = "http://cran.us.r-project.org")
```

To use the installed packages, specify the path when loading library:

```R
library(DNAcopy, lib.loc =  package_library)
```

# Create conda environment

Modified by Jani Huuhtanen based on the info by Oscar Brück in 2020

## Conda environment - Why?

*TL; DR*: Python works best with virtual environments. Anaconda (or Conda) enables easier and smarter building of virtual environments than “virtualenv”.

*Long explanation*: “Anaconda” is a free and open-source distribution of the Python and R programming languages for scientific computing, that aims to simplify package management and deployment. Thus, Anaconda builds a virtual environment, where your R and Python packages work well together each time the user installs or updates packages. It is also easy to copy entire virtual environments to colleagues or to another computer or computational enviroment.

Another virtual environment manager, “virtualenv” (a.k.a. pip package manager), is popular among Python users but has important caveats. The big difference between Conda and the pip package manager is in how package dependencies are managed, which is a significant challenge for Python data science and the reason Conda exists. Pip installs all Python package dependencies required, whether or not those conflict with other packages you installed previously. So your working installation of, for example, Google Tensorflow, can suddenly stop working when you pip install a different package that needs a different version of the numpy library. More insidiously, everything might still appear to work but now you get different results from your data science, or you are unable to reproduce the same results elsewhere because you didn’t pip install in the same order.

Conda analyzes your current environment, everything you have installed, any version limitations you specify (e.g. you only want tensorflow >= 2.0) and figures out how to install compatible dependencies. Or it will tell you that what you want can’t be done. Pip, by contrast, will just install the thing you wanted and any dependencies, even if that breaks other things.

If you like to work with more graphical softwares, for example you prefer to work with RStudio rather than R, you can also download Anaconda Navigator in your computer. Anaconda Navigator is a desktop graphical user interface (GUI) included in Anaconda distribution that allows users to launch applications and manage conda packages, environments and channels without using command-line commands.

## For who is conda for?
* If you use R and have problems with meeting package requirements
* If you use R and need packages such as tensorflow or reticulate, which rely on Python
* If you use Python

## How?

First, install conda. 

```
conda env create --name mycondaenv -f=/csc/mustjoki/anaconda3/envs/conda_r3.6_python3.6.0_environment.yml
```

Activate conda environment with python 3.7

```
source /csc/mustjoki/anaconda3/bin/activate mycondaenv

```

After that install R and python packages that you need and are not included in the enviroment, for example

```
conda install -c r r-reticulate -y;
```

Install essential conda, R and python configurations

```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

Install R 3.5.1.

```
### CHANGE THIS IF YOU WANT ANOTHER R VERSION
conda install -c r r=3.5.1
conda install -c conda-forge r-essentials
```

After that install R and python packages you need

```
## Here is a list of R packages Oscar uses. Modify the list if you do not recognise or need a certain package.
## Most likely some packages become installed with another package (dplyr with tidyverse), but I didn't feel like optimizing the package list as you can run this in one loop
## The run will take most likely 15-60 min
conda install -c r r-reticulate -y;
conda install -c conda-forge keras -y;
conda install -c conda-forge r-magick -y;
conda install -c bioconda/label/in-conda-forge r-purrr -y;
conda install -c r r-devtools -y;
conda install -c r r-tidyverse -y;
conda install -c r r-dplyr -y;
conda install -c conda-forge xorg-libx11 -y;
conda install -c eumetsat fftw3 -y;
conda install -c anaconda cairo -y;
conda install -c conda-forge r-readxl -y;
conda install -c bioconda r-readr -y;
conda install -c r r-reshape2 -y;
conda install -c conda-forge r-uwot -y;
conda install -c conda-forge r-flexdashboard -y;
conda install -c r r-hmisc -y;
conda install -c conda-forge r-survminer -y;
conda install -c conda-forge r-openimager -y;
conda install -c r r-igraph -y;
conda install -c mdekstrand r-plotroc -y;
conda install -c r r-gplots -y;
conda install -c conda-forge r-circlize -y;
conda install -c r r-proc -y;
conda install -c grst phenograph -y;
conda install -c bioconda bioconductor-complexheatmap -y;
conda install -c conda-forge r-janitor -y;
conda install -c conda-forge r-tictoc -y;
conda install -c conda-forge r-corrplot -y;
conda install -c r r-glmnet -y;
conda install -c conda-forge r-keras -y;
conda install -c r rstudio -y;
conda install -c conda-forge libiconv -y; # was needed to open R
Then what?

```

## Activate the groups

Now you are ready to go!

If you need other packages, download them from https://anaconda.org/ by looking up a code snipet
For example, the code for R tidyverse is
```
conda install -c r r-tidyverse
```

Some packages do not get downloaded with conda. These will have to be downloaded the old way by opening R and running "install.packages" in R

```

# Open R
R
# Quit R (type this in R console)
quit("no")

```

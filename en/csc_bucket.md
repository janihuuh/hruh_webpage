Added April 17 2023 by Olli

# Downloading data from CSC bucket using command line

# 1. Download data using rclone to the FIMM cluster

Download data to our storage in FIMM cluster (/csc/mustjoki or /csc/mustjoki2) to be able to decrypt and extract compressed files. This is to check that the data can be decrypted and to be able to download e.g. scRNA-seq count matrices to your computer.

First create a directory for your data in /csc/mustjoki or /csc/mustjoki2 and run the following in that directory. Do this in a screen so that if your connection to the cluster breaks it doesn’t matter. Example for user dufvaoll with project project_2007414:

```bash
screen -S some_name_for_your_screen
source allas_conf -u dufvaoll -p project_2007414
rclone copy allas:nkheme_data_import_1 nkheme_data_import_1

# Exit screen (ctrl+A+D)
```

# 2. Create conda environment for using crytp4gh in the cluster

Buy Oscar Brück beer every time you use this part

```bash
/csc/mustjoki/anaconda3/bin/conda create --name mycondaenv python=3.7
source /csc/mustjoki/anaconda3/bin/activate mycondaenv
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
pip install crypt4gh
conda deactivate
```

# 3. Decrypt one file

Put your secret key file into the directory you created for the data in /csc/mustjoki (nkheme in the example). Go to the directory where your file is and run the following:

```bash
source /csc/mustjoki/anaconda3/bin/activate odcondaenv
crypt4gh decrypt --sk /csc/mustjoki2/nkheme/nkheme_crypt4gh.key < filename.tar.gz.c4gh > filename.tar.gz
```

# 4. Decrypt as batch if you have many samples

This is again good to do in a screen:

```bash
screen -r some_name_for_your_screen
source /csc/mustjoki/anaconda3/bin/activate mycondaenv
for f in ./*.c4gh; do crypt4gh decrypt --sk /csc/mustjoki2/nkheme/nkheme_crypt4gh.key < $f > "${f/%.tar.gz.c4gh/.tar.gz}"; done
```



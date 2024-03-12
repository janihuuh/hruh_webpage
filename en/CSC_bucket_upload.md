# Transfer of sensitive data to Allas

Written in January 2024 by Essi

## Create new screen 

```bash
screen -S <SOME_NAME_FOR_SCREEN>
```

## Install two Python packages (and conda environment if you don't have it yet)

**You can also skip this part and use someone else's Conda environment (if it includes the packages)**

I have the bash profile in homes/<USERNAME>/.bash_profile and I added PATH=$PATH:/csc/mustjoki/anaconda3/bin to it. Remember to source it (source /homes/<USERNAME>/.bash_profile)

```bash
/csc/mustjoki/anaconda3/bin/conda create --name <SOME_ENV_NAME> python=3.7
```

The default location for the new conda environment is (hidden) /homes/<USERNAME>/.conda. If you want to use a different location, do this: 

```bash
/csc/mustjoki/anaconda3/bin/conda create --prefix /csc/mustjoki/<YOUR_FOLDER>/<SOME_ENV_NAME> python=3.7
```

This takes some time (maybe 5-10 min) and it might look like nothing is happening

```bash
# Activate the environment
source /csc/mustjoki/anaconda3/bin/activate <SOME_ENV_NAME>


##Oscar suggested doing these but you might not need to (at least I got "Warning: 'bioconda' already in 'channels' list, moving to the top")
#conda config --add channels conda-forge
#conda config --add channels defaults
#conda config --add channels r
#conda config --add channels bioconda

# Install the needed Python packages
pip install crypt4gh
pip install s3cmd
```

## allas-cli-utils (to be able to use a-put, a-get etc.)

```bash
git clone https://github.com/CSCfi/allas-cli-utils
cd  allas-cli-utils/
export PATH=${PATH}:$(pwd)
```

Note: If you are in one of the fas2/NGS/... folders, you'll probably get "Permission denied". In that case, go to some of our own folders /csc/mustjoki/... until you reach the a-put step.
Alternatively, you can skip git clone and just export my path /csc/mustjoki/essi/allas-cli-utils

## Connect to allas

Connect to CSC Allas using swift based connection and choose our project

```bash
source allas_conf -u <CSC_USERNAME> -p project_2007704
```

Notice that the connection is open for 8 hours. I'm still not sure whether that can be easily extended so at least for now, it's best to transfer small enough subfolders (maybe max 400 GB). You probably don't want to store too huge objects in Allas anyway.

## Use a-put to encrypt and transfer the file to Allas

Francesca has `<OUR_KEY>`

Go to your folder and select a subfolder you want to transfer. It will be encrypted and compressed during the transfer. Think carefully what kinds of entities you want to transfer. Later you will retrieve one object (encrypted and compressed folder) at the time. The size should not be huge but on the other hand you don't want to have 1000 separate objects. I'd say objects of size 50-500 GB could be a good plan.

Important! Before you a-put a subfolder, think of an informative name for it. If your `<FOLDER_NAME>` is right now not very infomative, change it now: `mv old_folder_name new_folder_name`

You can use option --asis if you want to keep the folder structure such that it encrypts and compresses each file separately. I left it out because a-get is difficult for each file separately

```bash
a-put --tmpdir /homes/<USERNAME>/allas_tmp --encrypt c4gh --pk <OUR_KEY>.pub <FOLDER_NAME> -b Project_2007704_HRUH_pub/<NEW_ALLAS_FOLDER>
```

Substitute `<NEW_ALLAS_FOLDER>` with some informative project name. (Really new only if you are doing this for the first time for the project. In case you already have a project folder you want to add to, you can check its name with rclone lsd allas:Project_2007704_HRUH_pub/)

`<FOLDER_NAME>` is the informative subfolder name of the subfolder you are transferring. The end result will be `allas:Project_2007704_HRUH_pub/<NEW_ALLAS_FOLDER>/<FOLDER_NAME>.tar.c4gh`

## Check that the transfer happened

```bash
a-check --encrypt c4gh --pk <OUR_KEY>.pub <FOLDER_NAME> -b Project_2007704_HRUH_pub/<NEW_ALLAS_FOLDER>
```

List the content of your project folder: 

```bash
rclone ls allas:Project_2007704_HRUH_pub/<NEW_ALLAS_FOLDER>/
```

Does the size of the object look reasonable? It should be quite close to that of the original folder

## Get and decrypt the data

```bash
a-get --sk ../<OUR_KEY>.sec Project_2007704_HRUH_pub/<NEW_ALLAS_FOLDER>/<FOLDER_NAME>.tar.c4gh
```

Ask Francesca about the password

## When you are ready with the transfer and the checking

```bash
conda deactivate
```

Ctrl a d to exit the screen, then: 

```bash
screen -X -S <SOME_NAME_FOR_SCREEN> kill
```

**Importantly: Describe the data and the subfolder structure in the Google sheet. The link was shared through Slack (or you can ask Essi or anyone who has done this before)**

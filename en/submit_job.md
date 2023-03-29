# Submitting jobs to the job queue (adapted from Matti Kankainen)

Use the job queue for all longer running computations.

Let us imagine that I wanted to run command ‘ls –lha’ to see all permissions to current directory. Usually, this command could of course be executed on a command line but let us now imagine that this command would be a longer running computation, thus needed to execute through job queue.


#### 1. Create a shell script (.sh-ending file)  
```bash
/projects/fimm_ngs_mustjoki/tcr/run.sh 
```
You can go for example to directory /projects/fimm_ngs_mustjoki/tcr/ with command
cd /projects/fimm_ngs_mustjoki/tcr/. After that, you can make a shell script named run.sh with command vim run.sh. Now there should be a file run.sh in directory /projects/fimm_ngs_mustjoki/tcr/. In reality, this could be any directory, and shell script could have any name as long as it ends with .sh.


#### 2. Edit permissions to your script to be able to run it
At least I can skip this part -Timo
```bash
chmod +x run.sh 
```


#### 3. Write commands to shell script.
We want to execute command ```ls –lha``` so let us write it to our file.


#### 4. Submit job to queue
This can be done 3 ways


##### 4.1 Submit job with grun.py
```bash
grun.py -n xyz -q highmem.q -c "/projects/fimm_ngs_mustjoki/tcr/run.sh"
```

* -n xyz = xyz is a prefix for error and stdout files and the name for the run in the queue
* -c ”pathtoscript” = path to your script to be submitted to the queue
* -q highmem.q = queue where the job is submitted


##### 4.2 Submit job with qsub specifying commands in command line
```bash
qsub -N xyz -cwd -q all.q /projects/fimm_ngs_mustjoki/tcr/run.sh
```

* qsub = Command that submits your job to queue
* N xyz = xyz is a prefix for error and stdout files and the name for the run in the queue
* cwd = Current working directory, without this submitted command (ls -lha) will be executed in your home directory, and all reports are also created there. Now everything happens in directory you are currently in.
* q all.q = queue where the job is submitted. Different queues are listed later.
* /projects/fimm_ngs_mustjoki/tcr/run.sh = path to your script. If the script is located in the current directory, you can only type the file name. This should be the last argument.


##### 4.3 Submit job with qsub specifying commands in shell script
```bash
qsub run.sh
```
You can write all commands like -N xyz and -q all.q to your run.sh. This makes it more clear for at least me. Write commands to different lines, starting lines with #$. If you want to specify for instance queue in your file, you can write #$ -q all.q in the beginning of that file. Then you can submit your job by typing qsub run.sh to command line. Example run.sh file below.

```bash
#$ -N xyz
#$ -q all.q
#$ -cwd
#$ -e /projects/fimm_ngs_mustjoki/tcr/error_file.txt
#$ -o /projects/fimm_ngs_mustjoki/tcr/stdout_file.txt

ls -lha
```
-e and -o arguments are optional and specify STDOUT and error reports. More about them in the next section.


#### 5. Files created after submitting files


##### 5.1 Files created after grun.py
If you used grun.py, after the run has started, three files will appear to the directory you started the run in: 

* xyz.OU 	STDOUT of the program you are running, if not directed to other file by the program
* xyz.ER	error report of the program you are running, if not directed to other file by the program
* xyz.job		run start command


##### 5.2 Files created after qsub
If you submitted jobs using qsub, when the run has started, two files will appear to the directory you started the run in:

1. STDOUT file of the program you are running (everything command line would usually print). If you have not specified it, file name is something like xyz.o6213656, where numbers vary. Result of ls -lha should be found inside this file.

2. Error report of the program you are running. For example, if you have typed a no existing command, error text is found inside this file. If you have not specified error file name, it is something like xyz.e6213656, where numbers vary.

You can specify names of STDOUT and error files with -e and -o commands, as seen in the example file above.


#### 6. Monitoring your job


##### 6.1 Job status with
```bash
qstat
```

In qstat output:

	job-ID 	= ID of your job, can be used to terminate the job
	prior	= 0.5 for all
	name	= the job name you specified, xyz in example
	user	= you
state	= r (running), Eqw (error, job will not start), q (queued), qw (job in queue waiting to be scheduled)


##### 6.2 Job memory usage
To find the memory usage of your job run, use
```bash
qstat -j <job-id> | grep maxvmem
```


#### 7. Deleting your job

If you made a mistake and want to terminate the process

```bash
qdel JOBID
```
E.g. qdel 1205638


#### 8. Running multiple jobs

If you need to run several jobs, e.g. same script for 4 different datasets, create 4 scripts and start them separately. If your task takes long, it is useful to break it into multiple parts (if possible) to be run in parallel so that it will finish sooner.

```bash
grun.py -n xyz1 -q highmem.q -c "/projects/fimm_ngs_mustjoki/tcr/run1.sh"
grun.py -n xyz2 -q highmem.q -c "/projects/fimm_ngs_mustjoki/tcr/run2.sh"
grun.py -n xyz3 -q highmem.q -c "/projects/fimm_ngs_mustjoki/tcr/run3.sh"
grun.py -n xyz4 -q highmem.q -c "/projects/fimm_ngs_mustjoki/tcr/run4.sh"
```

#### 9. Different queues and memory usage
If you need more memory, use –q hugemem.q

Update 11/2019: New flags for grun.py

```bash

-M, --memory [value]     Reserved memory. This value is in GB, max of 50

                                      Defaults depend on queue, and is generally the old limits.

-C, --cores [value]         Replaces the --slots flag as slots are now unbound from CPU-cores.

                                      The amount of CPU-cores to reserve. Currently max of 4.

                                      Defaults depend on queue.

--dry-run                        Switch, creates the job-file but does not queue it to the cluster.

                                      For debugging or manual editing of the job-files before submitting.

-d, --debug                    Switch, prints debug information


For those of you who use qsub you can request hardware via,

-l h_rss=[memory]         In bytes, supports K M G as value multipliers (for example 4G)

-l res_cores=[value]       Number of cores to reserve.


The defaults per queue are as follows:

Queue:               Default cores/memory/max memory
all.q                    1Core/2GB/3GB

medmem.q         2Cores/5GB/10GB

highmem.q         3Cores/10GB/20GB

hugemem.q        4Cores/20GB/45GB

express.q            2Cores/2GB/2GB /Max 90minutes runtime.

interactive.q        2Cores/16GB/16GB

sisu.q                  1Core/7GB/NA

sisu_impute.q     1Core/1.7GB/NA

Queues not listed are for special use only and does not currently have any default limits.

The long term goal is to move away from the queues and instead only use resource requests for jobs.
This also ties into being able to run singularity containers in the cluster.
Pure docker containers due to its security limitations are not currently planned.
```
Some "hidden" nodes:

```
If you need even more memory (up to 256G), use –q seq_hugemem2.q
Working interactively in a node

If you want to work interactively without using the job queue, you can enter a node and work there without disturbing others in the main node.
ssh nameofnode
e.g. ssh compute-16

```

#### 10. More information and commands for qsub
More information and commands for qsub can be found in http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html.

## Running R scripts in job queue

Running R scripts in job queue

#### 1.	Make an R script e.g. 

```bash
/projects/fimm_ngs_mustjoki/tcr/analysis.R
```

#### 2.	Make a bash script e.g. “run_R.sh” containing the following:

```R
/apps/statistics2/R-3.5.1/bin/R --no-save < /projects/fimm_ngs_mustjoki/tcr/analysis.R >& messages.err
```

messages.err will contain the error, warning etc. messages from R.

#### 3.	Submit the bash script to queue as above.

##### To install R packages, specify a path in the project directory and save it into libpath variable:

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

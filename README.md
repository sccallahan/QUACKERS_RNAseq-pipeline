QUACKERS (an RNAseq pipeline)
================
Carson Callahan
Feb 2, 2021

Intro
-----

![Desperate times, desperate measures.](http://static.makeuseof.com/wp-content/uploads/2017/12/weird-programming-principles-670x335.jpg "Desperate times, desperate measures.")

The pipelines our lab typically uses for NGS have recently been broken by changes to our HPC. While running the tools in a standalone fashion isn't particularly complicated and ultimately results in the same output, allowing users to run the entire workflow using by supplying arguments to a single script is generally more convenient and allows less experienced users to run their workflows unassisted... plus it saves some time on the user end by not needing to babysit individual jobs.

With this in mind, I've created a small, simple pipline called QUACKERS - **QU**ickly **A**ssembled **C**ode that **K**nits **E**ssential **R**NAseq **S**cripts (though it's less "quickly assembled" as time goes on...). In short, users supply arguments to a parent script, which then passes these arguments to child scripts that submit jobs to the cluster. Currently, the workflow will align fastqs using STAR, generate a multiqc report using fastqc outputs (from bams) and STAR outputs, and generate a counts matrix using featureCounts and [some custom bits of code from Ming (Tommy) Tang.](https://divingintogeneticsandgenomics.rbind.io/post/merge-featurecount-table-from-rnaseq/)

Setting up the tools
--------------------

This pipeline relies on STAR for alignments and featureCounts for counting aligned reads. Both of these tools require external files, the paths of which must be provided (covered later in Running the Code). In particular, you need the .gtf files and index generated from STAR and a .gtf file for featureCounts (this should generally be the same as the one used to generate the STAR index).

Distribute workflow
-------------------

``` bash

# connect to Shark
ssh username@shark.mdanderson.edu # alternatively, `ssh shark` if on a hard connection

# make a folder to hold all the files - name it whatever you like
mkdir workdir
cd workdir

# clone the repo
source activate py351 # need this for `git clone` and `multiqc` later
git clone https://github.com/sccallahan/QUACKERS_RNAseq-pipeline
cd QUACKERS_RNAseq-pipeline
```

This should deposit all the scripts into your directory. Below is a list of scripts you should now have in the same (working) directory:

-   `beginQuacking.sh`
-   `quack_mergeFastq.sh`
-   `quack_STAR.sh`
-   `quack_fastqc.sh`
-   `quack_multiqc.sh`
-   `quack_mergeCounts.sh`

Fastq directory structure
-------------------------

The code assumes a directory structure as follows:

-   Fastq\_dir
    -   Sample1
        -   Sample1\_R1.fastq.gz
        -   Sample1\_R2.fastq.gz
    -   Sample2
        -   Sample2\_R1.fastq.gz
        -   Sample2\_R2.fastq.gz

and so on.

Running the Code
----------------

Simply run the parent script `beginQuacking.sh` In brief, this script is run with the following syntax:

`bash beginQuacking.sh [fastq directory] [index path] [genome gtf file] [sleep time]`

Here's a brief explanation of the 4 arguments:

-   fastq directory - a directory with the previously mentioned structure, containing only folders with fastqs
-   index path - the *folder* containing the STAR index files
-   genome gtf file - path to the gtf file(s) used to generate the STAR index, or the gtfs to be used on the fly during mapping
-   sleep time - amount of time to wait before checking for files during aligment, counting, etc. For example, 60 would tell the script to check for completion of alignment every 60 seconds and print a message to the terminal. In most cases, 60-180 seconds is probably best.

So, running the script for me would look something like this:

``` bash
bash beginQuacking.sh /rsrch2/headneck/sccallahan/fastqs/ /rsrch2/headneck/sccallahan/apps/annot/starbase/hg19/starindex /rsrch2/headneck/sccallahan/apps/annot/starbase/hg19/annot/file.gtf 60
```

This will submit all the jobs. The scripts will make folders in the current workding directory and place relevant outputs there. The created folders are:

-   logs - has .out and .err logs
    -   lsf\_logs - the files containing the job submission scripts
-   fastqc\_output - fastqc outputs
-   multiqc\_output - multiqc outputs
-   STAR\_output - STAR outputs
-   featureCounts\_output - raw counts files for each sample, as well as a merged table with the counts for all samples.

Job control
-----------

The main control of interest is killing the jobs. This can be done as follows:

``` bash
bkill ` bjobs -u sccallahan |grep RUN |cut -f1 -d" "` # replace my username with yours, can also change RUN to PEND or relevant status
```

Downstream processing
---------------------

For downstream work, you should first check the multiqc output to make sure your sequencing was of sufficient quality (uniquely mapped reads, duplication rates, etc.).

Next, read the raw counts table into R. You'll notice the column names are very messy, but this can be easy cleaned up by manually setting column names, or being clever with regex to do it for you. Now you should have a counts matrix that is ready for analysis with DESeq2 or other tools.

QUACKERS (an RNAseq pipeline)
================
Carson Callahan
July 2, 2019

Intro
-----

![Desperate times, desperate measures.](http://static.makeuseof.com/wp-content/uploads/2017/12/weird-programming-principles-670x335.jpg "Desperate times, desperate measures.")

The pipelines our lab typically uses for NGS have recently been broken by changes to our HPC. While running the tools in a standalone fashion isn't particularly complicated and ultimately results in the same output, allowing users to run the entire workflow using by supplying arguments to a single script is generally more convenient and allows less experienced users to run their workflows unassisted... plus it saves some time on the user end by not needing to babysit individual jobs.

With this in mind, I've created a small, very simple pipline called QUACKERS (QUickly Assembled Code that Knits Essential RNAseq Scripts). In short, users supply arguments to a parent script, which then passes these arguments to child scripts that submit jobs to Shark. Currently, the workflow will quality check fastqs, generate a multiqc report, align fastqs using STAR, then generate a counts matrix using featureCounts and [some custom bits of code from Ming (Tommy) Tang.](https://divingintogeneticsandgenomics.rbind.io/post/merge-featurecount-table-from-rnaseq/)

Currently, this pipeline only supports human (hg19). I hope to add support for mm9 and generation of bigwig files for vizualization at a later date.

Setting up the tools
--------------------

This pipeline relies on STAR for alignments and featureCounts for counting aligned reads. Both of these tools require external files, the paths of which must be provided (covered later in Running the Code). In particular, you need the .gtf files and index from STAR (these come pre-packged) and a genome annoation file for featureCounts (this is not provided, but can be attained from UCSC - I will also upload my hg19 copy to this repo). The scripts are currently set up to take the path to the STAR hg19 folder, then navigate to the .gtf files and starindex on its own (assuming you haven't modified the setup - you could change the `Quackling_STAR.sh` to whatever paths you need).

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
git clone https://gitlab.com/sccallahan/quackers---an-rnaseq-pipeline
cd quackers---an-rnaseq-pipeline
```

This should deposit all the scripts into your directory. Next, you need to create a sample sheet. This sheet MUST BE NAMED `sample_sheet_STAR.txt`, or the pipeline will fail. Alternatively, you could edit the `Quackling_STAR.sh` file to whatever you named your sample sheet. The sample sheet needs to be a tabbed texted file, with each paired end fastq belonging to the same row and separate columns. Here's an example:

``` bash
HN30_KD1_1.fq.gz    HN30_KD1_2.fq.gz
HN30_KD2_1.fq.gz    HN30_KD2_2.fq.gz
HN30_KD3_1.fq.gz    HN30_KD3_2.fq.gz
HN30_NT1_1.fq.gz    HN30_NT1_2.fq.gz
HN30_NT2_1.fq.gz    HN30_NT2_2.fq.gz
HN30_NT3_1.fq.gz    HN30_NT3_2.fq.gz
```

You should now have the following files in the same (working) directory:

-   `QUACKERS.sh`
-   `Quackling_fastqc.sh`
-   `Quackling_STAR.sh`
-   `sample_sheet_STAR.sh`

And you should know the paths to the hg19 folder for STAR and have a gene annoation file for featureCounts. You should also have all your fastq files in their own directory (named whatever you like).

**NB:** The pipeline cannot currently recursively search folders, i.e., you'll need to move each .fastq.gz into the same parent directory (some cores provide fastq.gz directly, some give you each folder containing each pair of fastq.gz as they come off the sequencing machine). I plan on adding this functionality when I get some time to work on it. In the mean time, I've uploaded a hack-y workaround called `extract_fastq_from_folder.sh`. Place this script in your fastq directory and run it - it will move into each folder and rsync the fastq.gz into the parent directory, provided your directory structure is:

-   Fastq\_dir
    -   Fastq\_sample\_1\_dir
        -   Fastq\_files
    -   Fastq\_sample\_2\_dir
        -   Fastq\_files

and so on. After the pipeline is done, you can freely delete the rsync'd fastqs to save disk space if you'd like.

Running the Code
----------------

Simply run the parent script `QUACKERS.sh` In brief, this script is run with the following syntax:

`bash QUACKERS.sh [path/to/fastqs] [fastq extension] [path/to/STAR/hg19/folder] [path/to/gene/annotation/file]`

Here's a brief explanation of the 4 arguments:

-   Path to fastqs - Point this to a folder containing only the fastqs of interest
-   fastq exntension - This takes everything after the "." in the file name. For example, my ".fq.gz" would be entered simply as "fq.gz"
-   path to STAR hg19 - FULL PATH path to the STAR hg19 folder. Example, mine is /rsrch2/headneck/sccallahan/apps/annot/starbase/hg19
-   path to annotation file - FULL PATH to the UCSC annoation file. Example, mine is /rsrch2/headneck/sccallahan/common\_files/genes.gtf

So, running the script for me would look something like this:

``` bash
bash QUACKERS.sh /rsrch2/headneck/sccallahan/testdir/testfqs fq.gz /rsrch2/headneck/sccallahan/apps/annot/starbase/hg19 /rsrch2/headneck/sccallahan/common_files/genes.gtf
```

This will submit all the jobs. The scripts will make folders in the current workding directory and place relevant outputs there. The created folders are:

-   logs - has .out and .err logs
-   fastqc\_output - outputs of JUST fastqc
-   multiqc\_output - multiqc output, the html file is of particular interest
-   STAR\_output - all the outputs from STAR, including BAMs
-   featureCounts\_output - raw counts files for each sample, as well as a merged table with the counts for all samples.

Job control
-----------

The main control of interest is killing the jobs. This can be done as follows:

``` bash
bkill ` bjobs -u sccallahan |grep RUN |cut -f1 -d" "` # replace my username with yours, can also change RUN to PEND or relevant status
```

Downstream processing
---------------------

For downstream work, you should first check the multiqc output to make sure your sequencing was of sufficient quality. Next, check the STAR and featureCounts outputs to see how many reads aligned to the genome and how many mapped to the transcriptome. Ideally, you'd have at least 18-20M PE reads mapping to the transcriptome.

Next, read the raw counts table into R. You'll notice the column names are very messy, but this can be easy cleaned up by manually setting column names, or being clever with regex to do it for you. I don't have too many samples, so it's easy to set manually.

``` r
library(tidyverse)

counts <- read.table("/path/to/raw_counts_table.txt", header = TRUE)
counts <- column_to_rownames(counts, var = "Geneid")
colnames(counts) <- c("KD1", "KD2", "KD3", "NT1", "NT2", "NT3") # you could also try some regex here, which is useful if you have tons of samples
```

Now you should have a counts matrix that is ready for analysis with DESeq2 or other tools.

To Do
-----

-   Add support for mm9
-   Include bigwig generation in workflow

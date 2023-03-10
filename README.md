
# Utility collection of code

This repository is a collection of utility scripts for performing various tasks of automation, stats collection, file merging and more. 

## Scripts are as follows:

1. DESeq2_automate.R<br>
2. edgeR.R<br>
3. discardoverview.py<br>
4. extract-resolved-raw-fastq.py<br>
5. convertbamstat2csv.py<br>
6. volcano_plots.R<br>
7. hpfinder.py <br>
8. tmux-session.sh <br>



### DESeq2 automation
#### Prerequisites
R 3.5 or newer<br>
DESeq2<br>
edgeR

This script generates tables and plots for DEGs. To run this script please provide a count and experiment file. Example datasets are also attached.Biomart file is also uploaded in this repo.

### edgeR 
To run paired sample analysis and find out DEGs.

### discardoverview
This scripts collects and combine the stats from bamstat.pl and enigmaff,enigmarr filtering. 

Its a globber script which scans for `*bwameth.bamtools_stats.txt` and `*bwameth.enigmaff.bamtools_stats.txt`
and writes a csv format output

### extract-resolved-raw-fastq
This script will read resolved fastq file from SSPW run directory and save the ids then it will look for original raw fastq files and just extract the sequences which were resolved. It will write a new pair of resolved fastq files.


### volcanoplots
To create volcano plots on DEG genes with FDR filtering and it will run on all the comparisons from edgeR.

### hairpin finder in fastq reads
To detect a hairpin seq or its fragments in fastq pair file

### to save and restore tmux sessions
To prevent tmux sessions disappear when we restart a machine
bash tmux-session save
bash tmux-session restore



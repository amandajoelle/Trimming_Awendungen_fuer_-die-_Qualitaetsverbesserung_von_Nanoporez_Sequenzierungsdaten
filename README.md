# nextflow-test

Pipeline to analyse and make a data quality control of an RNA-Genom using two different assembler frameworks

## First Step

### Install Conda
https://linuxhint.com/install-anaconda-ubuntu-22-04/ 
### Install Prinseq, Flye, Quast.
#### Configure conda
> conda config –add channels bioconda 
> conda config –add channels conda-forge 
> conda config –show channels  
#### Install Priseq
$ wget -N http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz

$ tar -zxvf prinseq-lite-0.20.4.tar.gz

$ cp -puv prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq-lite && chmod +x /usr/local/bin/prinseq-lite

$ cp -puv prinseq-lite-0.20.4/prinseq-graphs.pl /usr/local/bin/prinseq-graphs && chmod +x /usr/local/bin/prinseq-graphs

#### Install Flye and Quast

> conda create –n genome_assembly flye  quast  
#### Activate the environment

> conda env list 

>conda activate genome_assembly 

### Use
    - Prinseq 
> prinseq-lite -fastq <filename>.<fastq | fastq.gz> -out_format 1 

EXAMPLE : prinseq-lite.pl -fastq 1_control_psbA3_2019_minq7.fastq -out_format 1


    - Flye 
> flye --nano-raw <filename>.<fasta> -o assembly –g 5.6m -t 10 –i 2 

EXAMPLE :

 flye --pacbio-raw 1_control_psbA3_2019_minq7_prinseq_good_OqHZ.fasta --out-dir out_p --threads 4


    - Quast 

> quast <filename>.<fasta> <filename>.<fasta> -o <directory-name>

EXAMPLE : quast -o quast_report E.coli_PacBio_40x.fasta prinseq_good_OqHZ.fasta

    
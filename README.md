# Canu+HISEA

This is modified `(Canu+HISEA)` assembly pipeline. [HISEA](https://github.com/lucian-ilie/HISEA) is an efficient all-vs-all long read aligner for SMRT sequencing data. Its algorithm is designed to produce highest alignment sensitivity among others. The [HISEA](https://github.com/lucian-ilie/HISEA) program has been integrated in Canu pipeline in order to produce better assemblies. For details of HISEA program, please see [HISEA](https://github.com/lucian-ilie/HISEA)

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://www.nanoporetech.com/products-services/minion-mki)).

This is modified Canu assembly pipeline `(Canu+HISEA)` which runs in four steps, similar to originial Canu pipeline:

* Detect overlaps in high-noise sequences using [HISEA](https://github.com/lucian-ilie/HISEA) (This is new in our pipeline)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

The evaluation of HISEA genome assembly is performed using 30X and 50X sub-sampled data extracted from original datasets downloaded from [Pacific Biosciences DevNet Datasets](https://github.com/PacificBiosciences/DevNet/wiki/Datasets). The 30X and 50X coverage datasets were sampled using the utility fastqSample available from the Canu pipeline.

The HISEA configuration files used for 30X Canu assembly pipeline can be downloaded from the table below. For 50X, same configuration file was used with modification of __corHiseaSensitivity__ parameter set to __normal__.

Genome | Configuration File Links
:--- | :--- 
E.coli | [E.coli configuration](http://www.csd.uwo.ca/faculty/ilie/HISEA/conf_files/ecoli_conf.txt) 
S.cerevisiae | [S.cerevisiae configuration](http://www.csd.uwo.ca/faculty/ilie/HISEA/conf_files/scerevisiae_conf.txt)
C.elegans | [C.elegans configuration](http://www.csd.uwo.ca/faculty/ilie/HISEA/conf_files/celegans_conf.txt)
A.thaliana | [A.thaliana configuration](http://www.csd.uwo.ca/faculty/ilie/HISEA/conf_files/arabidopsis_conf.txt)
D.melanogaster | [D.melanogaster configuration](http://www.csd.uwo.ca/faculty/ilie/HISEA/conf_files/drosophila_conf.txt)

A detailed comparison of HISEA with other leading programs can be found in HISEA paper (Nilesh Khiste and Lucian Ilie). Below are some plots showing NGA50 results for Canu+HISEA and Canu+MHAP.

## NGA50 Comparisons 

<img src="http://www.csd.uwo.ca/faculty/ilie/HISEA/images/nga50_ecoli.jpg" width="480" height="288" alt="E.coli"> | <img src="http://www.csd.uwo.ca/faculty/ilie/HISEA/images/nga50_scerevisiae.jpg" width="480" height="288" alt="S.cerevisiae"> 
--- | --- 
<img src="http://www.csd.uwo.ca/faculty/ilie/HISEA/images/nga50_celegans.jpg" width="480" height="288" alt="C.elegans"> | <img src="http://www.csd.uwo.ca/faculty/ilie/HISEA/images/nga50_Arabidopsis.jpg" width="480" height="288" alt="A.thaliana">
<img src="http://www.csd.uwo.ca/faculty/ilie/HISEA/images/nga50_droso.jpg" width="480" height="288" alt="D.melanogaster"> |


## Build:

    git clone https://github.com/lucian-ilie/Canu_HISEA.git
    cd Canu_HISEA/src
    make -j <number of threads>

## Run:

Brief command line help:

    ../<achitechture>/bin/canu
    

Full list of parameters:

    ../<architecture>/bin/canu -options

```javascript 
HISEA specific parameters for configuration file:

<tag>HiseaBlockSize
Chunk of reads that can fit into 1GB of memory. Combined with memory to compute the size
of chunk the reads are split into.

<tag>HiseaMerSize
Use k-mers of this size for detecting overlaps.

<tag>HiseaMemory
Memory size per block.

<tag>HiseaSensitivity
Either normal, high, or low

Here is an example of a dummy configuration file, config.txt:

corOverlapper=hisea
corHiseaMerSize=16
corHiseaSensitivity=high
corHiseaMemory=200
corHiseaBlockSize=20000
corOvlRefBlockSize=20000
useGrid=0
```        
    
## Learn:

For usage specific to HISEA configuration, please look at our [webpage](http://www.csd.uwo.ca/faculty/ilie/HISEA/). The [quick start](http://canu.readthedocs.io/en/stable/quick-start.html) will get you assembling quickly, while the [tutorial](http://canu.readthedocs.io/en/stable/tutorial.html) explains things in more detail. 

## Citation:

Coming soon - to be updated
 

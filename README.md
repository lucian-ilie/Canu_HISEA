# Canu

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://www.nanoporetech.com/products-services/minion-mki)).

This is modified Canu assembly pipeline which runs in four steps, similar to originial Canu pipeline:

* Detect overlaps in high-noise sequences using [HISEA](https://github.com/lucian-ilie/HISEA) (This is new in our pipeline)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

## Build:

    git clone https://github.com/lucian-ilie/Canu_HISEA.git
    cd Canu_HISEA/src
    make -j <number of threads>

## Run:

Brief command line help:

    ../<achitechture>/bin/canu
    

Full list of parameters:

    ../<architecture>/bin/canu -options
    
## Learn:

For usage specific to HISEA configuration, please look at our [webpage](http://www.csd.uwo.ca/faculty/ilie/HISEA/). The [quick start](http://canu.readthedocs.io/en/stable/quick-start.html) will get you assembling quickly, while the [tutorial](http://canu.readthedocs.io/en/stable/tutorial.html) explains things in more detail. 

## Citation:

 

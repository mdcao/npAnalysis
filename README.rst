=======================================================
Real-time analysis of Oxford Nanopore MinION sequencing
=======================================================


This directory contains information for setting up real-time analysis
of Oxford Nanopore sequencing data, as described in the paper:

Real-time strain typing and analysis of antibiotic resistance potential using Nanopore MinION sequencing

Minh Duc Cao, Devika Ganesamoorthy, Alysha Elliott, Huihui Zhang, Matthew Cooper, Lachlan Coin

bioRxiv doi: http://dx.doi.org/10.1101/019356

---------------------
Software installation
---------------------

The streamline agorithms and auxiliary programs for setting up the pipeline are provided 
in japsa package (https://github.com/mdcao/japsa). In addition, the following 
free-software dependencies are required:

1. Java >=1.7
2. bwa >=7.10 (0.7.10-r858-dirty recommended)
3. R with rJava and MultinomialCI installed
4. kalign2 (http://msa.sbc.su.se/cgi-bin/msa.cgi)
5. HDF5 (https://hdfgroup.org/HDF5/release/obtain5.html) -- optional, needed for translating native fast5 files to fastq (npReader (https://github.com/mdcao/npReader).
 
++++++++++++++++++++++++ 
Quick installation guide
++++++++++++++++++++++++ 
::
   git clone https://github.com/mdcao/japsa
   cd japsa
   make install \
     [INSTALL_DIR=~/.usr/local \]
     [MXMEM=7000m \]
     [SERVER=true \]
     [JLP=/usr/lib/jni:/usr/lib/R/site-library/rJava/jri]



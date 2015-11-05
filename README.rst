#######################################################
Real-time analysis of Oxford Nanopore MinION sequencing
#######################################################

.. sectnum::

This directory contains information for setting up real-time analysis
of Oxford Nanopore sequencing data, as described in the paper:

Real-time strain typing and analysis of antibiotic resistance potential using Nanopore MinION sequencing

Minh Duc Cao, Devika Ganesamoorthy, Alysha Elliott, Huihui Zhang, Matthew Cooper, Lachlan Coin

bioRxiv doi: http://dx.doi.org/10.1101/019356

=====================
Software installation
=====================

The streamline agorithms and auxiliary programs for setting up the pipeline are provided 
in japsa package (https://github.com/mdcao/japsa). In addition, the following 
free-software dependencies are required:

1. Java >=1.7
2. bwa >=7.10 (0.7.10-r858-dirty recommended)
3. R with rJava and MultinomialCI installed
4. kalign2 (http://msa.sbc.su.se/cgi-bin/msa.cgi)
5. HDF5 (https://hdfgroup.org/HDF5/release/obtain5.html) -- optional, needed for translating native fast5 files to fastq (npReader (https://github.com/mdcao/npReader).
 
------------------------
Quick installation guide
------------------------
::

   git clone https://github.com/mdcao/japsa
   cd japsa
   make install \
     [INSTALL_DIR=~/.usr/local \]
     [MXMEM=7000m \]
     [SERVER=true \]
     [JLP=/usr/lib/jni:/usr/lib/R/site-library/rJava/jri]

(Note: the directives in squared brackets are optional. If you use those, remove the brackets and set the values accordingly to suit your computer settings).

This will install japsa according the directives:

* *INSTALL_DIR*: specifies the directory to install japsa, make sure you have write permission to this directory
* *MXMEM*: specifies the default memory allocated to the java virtual machine
* *SERVER*: specifies whether to launch the java virtual machine in server mode
* *JLP*: points to where HDF libraries and JRI are installed (e.g, /usr/local/lib:/usr/lib/R/site-library/rJava/jri). The path to HDF is only needed for creating the pipelines to analyse directly from fast5 files or simultaneously  from the MinION sequencing. For your convenience, add INSTALL_DIR/bin to your PATH envirenment e.g., ::

   export PATH=~/.usr/local/bin:$PATH

For more detailed information for installing japsa, please refer to  Japsa installation guide on
http://japsa.readthedocs.org/en/latest/install.html

==================
Databases and data
==================

The analyses described in the paper require accessing to some pre-processed databases. We make avaibale these 
databases on http://genomicsresearch.org/public/researcher/npAnalysis/ (and a backup storage 
https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/).
Setting up these databases for use as follows. You can choose to download the databases relevant to your desired analyese.

--------------------------------
Bacterial species identification
--------------------------------

We pre-compile the database of all bacterial genomes obtained from NCBI genbank, with the 
addition of two K. quasipneumoniae strains (to be updated in the manuscript). Download the 
database (~2.8GB), and make an bwa index of the database as follows.::

   wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/SpeciesTyping.tar.gz
   tar zxvf SpeciesTyping.tar.gz
   cd SpeciesTyping/Bacteria/
   gunzip genomeDB.fasta.gz
   bwa index genomeDB.fasta

Note that it might take a while to build the bwa index for this 9G-base database.

-----------------------
Strain typing with MLST
-----------------------

The database for MLST typing for three species,  K. pneumoniae, E. coli and
S. aureus are make avaibale. Download (208KB) and unzip them::

   wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/MLST.tar.gz
   tar zxvf MLST.tar.gz


--------------------------------------------
Strain typing with gene presence and absence
--------------------------------------------

The database for gene presence and absence strain typing for K. pneumoniae, E. coli and
S. aureus can be obtained as follows::

  wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/StrainTyping.tar.gz
  tar zxvf StrainTyping.tar.gz

------------------------------  
Resistance gene identification
------------------------------

A database of antibiotic resistance gene obtained from resFinder (https://cge.cbs.dtu.dk/services/ResFinder/) and pre-processed and provided from::

  wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/ResGene.tar.gz
  tar zxvf ResGene.tar.gz
  
  
======================================
Setting up real-time analysis pipeline
======================================

The framework makes use of the interprocess communication mechanism pipe (https://en.wikipedia.org/wiki/Pipeline_(Unix) ) as well as network channels to set up the real-time pipeline. You can prepare one or more analyses to run in real-time.












#######################################################
Real-time analysis of Oxford Nanopore MinION sequencing
#######################################################

.. sectnum::

This directory contains information for setting up real-time analysis
of Oxford Nanopore sequencing data, as described in the paper:

Streaming algorithms to identify pathogens and antibiotic
resistance potential from real-time MinION sequencing

Cao, M. D., Ganesamoorthy, D., Elliott, A. G., Zhang, H., Cooper, M. A., & Coin, L. J. M. (2016). Streaming algorithms for identification of pathogens and antibiotic resistance potential from real-time MinIONTM sequencing. GigaScience, 5(1), 32. http://doi.org/10.1186/s13742-016-0137-2


=====================
Software installation
=====================

The streamline agorithms and auxiliary programs for setting up the pipeline are provided
in japsa package (https://github.com/mdcao/japsa). In addition, the following
free-software dependencies are required:

1. Java >=1.8
2. bwa >=7.10 (0.7.10-r858-dirty recommended)
3. kalign2 (http://msa.sbc.su.se/cgi-bin/msa.cgi)
4. HDF5 (https://hdfgroup.org/HDF5/release/obtain5.html) -- optional, needed for translating native fast5 files to fastq (npReader (https://github.com/mdcao/npReader).

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
* *JLP*: points to where HDF libraries are installed (e.g, /usr/local/lib). The path to HDF is only needed for creating the pipelines to analyse directly from fast5 files or simultaneously  from the MinION sequencing. For your convenience, add INSTALL_DIR/bin to your PATH envirenment e.g.,::

   export PATH=~/.usr/local/bin:$PATH

For more detailed information for installing japsa, please refer to  Japsa installation guide on
http://japsa.readthedocs.org/en/latest/install.html

==================
Databases and data
==================

The analyses described in the paper require accessing to some pre-processed databases. We make avaibale these
databases on  http://data.genomicsresearch.org/Projects/npAnalysis/ (and a back up storage at https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/)
Setting up these databases for use as follows. You can choose to download the databases relevant to your desired analyese.

--------------------------------
Bacterial species identification
--------------------------------

The species identification pipeline requires a database of genomes of
interest which is simply the concatenation of all genomes in fasta format.
Prepare an index file which specifies the species of each sequence in the
database. For example:

Content of genomeDB.fasta::

  >NC_0000011 Chromosome of species Genus1 species1
  ACGTACGTACGT
  >NC_00000012 Plasmid 1  of species Genus1 species1
  ACGTACGTACGT
  >NC_00000013 Plasmid 2  of species Genus1 species1
  ACGTACGTACGT
  >NC_00000021 Chromsome of species Genus1 species2
  ACGTACGTACGT
  >NC_00000031 Chromsome of species Genus2 species3
  ACGTACGTACGT


Content of speciesIndex::

  Genus1_species1 >NC_0000011 Chromosome of species Genus1 species1
  Genus1_species1 >NC_00000012 Plasmid 1  of species Genus1 species1
  Genus1_species1 >NC_00000013 Plasmid 2  of species Genus1 species1
  Genus1_species2 >NC_00000021 Chromsome of of species Genus1 species2
  Genus2_species3 >NC_00000031 Chromsome of of species Genus2 species3


Finally, build a bwa index of the database::

  bwa index genomeDB.fasta


We pre-compile the database of all bacterial genomes obtained from NCBI genbank, with the
addition of two K. quasipneumoniae strains (to be updated in the manuscript). Download the
database (~2.8GB), and make an bwa index of the database as follows.::

   wget http://data.genomicsresearch.org/Projects/npAnalysis/SpeciesTyping.tar.gz
   tar zxvf SpeciesTyping.tar.gz
   gunzip SpeciesTyping/Bacteria/genomeDB.fasta.gz
   bwa index SpeciesTyping/Bacteria/genomeDB.fasta

Note that it might take a while to build the bwa index for this 9Gb database.

-----------------------
Strain typing with MLST
-----------------------


The database for MLST typing for three species,  K. pneumoniae, E. coli and
S. aureus are make avaibale. Download (208KB) and unzip them::

   wget http://data.genomicsresearch.org/Projects/npAnalysis/MLST.tar.gz
   tar zxvf MLST.tar.gz


--------------------------------------------
Strain typing with gene presence and absence
--------------------------------------------

The database for gene presence and absence strain typing for K. pneumoniae, E. coli and
S. aureus can be obtained as follows::

  wget http://data.genomicsresearch.org/Projects/npAnalysis/StrainTyping.tar.gz
  tar zxvf StrainTyping.tar.gz

------------------------------
Resistance gene identification
------------------------------

A database of antibiotic resistance gene obtained from resFinder (https://cge.cbs.dtu.dk/services/ResFinder/) and pre-processed and provided from::

  wget http://data.genomicsresearch.org/Projects/npAnalysis/ResGene.tar.gz
  tar zxvf ResGene.tar.gz


======================================
Setting up real-time analysis pipeline
======================================

The framework makes use of the `interprocess communication mechanism pipe <https://en.wikipedia.org/wiki/Pipeline_(Unix)>`_ as well as network channels to set up the real-time pipeline. The japsa package provides `jsa.util.streamServer <http://japsa.readthedocs.org/en/latest/tools/jsa.util.streamServer.html>`_ and `jsa.util.streamClient <http://japsa.readthedocs.org/en/latest/tools/jsa.util.streamClient.html>`_ to facilitate setting a pipeline distributed on a computer cluser. You can prepare one or more analyses to run in real-time.

For bacterial species typing::

   jsa.util.streamServer -port 3456 \
     | bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 SpeciesTyping/Bacteria/genomeDB.fasta - 2> /dev/null \
     | jsa.np.rtSpeciesTyping -bam - -index SpeciesTyping/Bacteria/speciesIndex --read 50 -time 60 -out speciesTypingResults.out 2>  speciesTypingResults.log &

This will create a pipeline to identify species which reports every 60 seconds, with at least 50 more reads from the last report. The pipeline waits for input on port 3456 for incoming data.


For strain typing gene presence/absense for K. pneumoniae::

   jsa.util.streamServer -port 3457 \
     | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a StrainTyping/Klebsiella_pneumoniae/geneFam.fasta - 2> /dev/null \
     | jsa.np.rtStrainTyping -bam -  -geneDB StrainTyping/Klebsiella_pneumoniae/ -read 0 -time 20 --out kPStrainTyping.dat 2>  kPStrainTyping.log &

You can run strain typing pipelines for other species (e.g., E. coli and S. aureus)
if you have reasons to believe the sample may contain these species. If these pipeline
run on the same computer, make sure they listen to different ports.

For strain typing with MLST::

   jsa.util.streamServer -port 3458 \ 
     | bwa mem -t 8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y MLST/Klebsiella_pneumoniae/bwaIndex/genes.fasta - \
     | jsa.np.rtMLST -bam - -mlst MLST/Klebsiella_pneumoniae/ -read 1000 -time 600  --out KpMLST.dat &

Again, you can set up MLST for E. coli and/or S. aureus as well. However, due to high error rate of the current 
Oxford Nanopore sequencing, this analysis may require a large amount of data. The presence/absence analysis above is recommended.


For resistance gene identification::

   jsa.util.streamServer -port 3459 \ 
     | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a ResGene/resFinder/DB.fasta - 2> /dev/null \
     | jsa.np.rtResistGenes -bam - -score=0.0001 -time 120 -read 50 --resDB  ResGene/resFinder/  -tmp _tmp_ -o resGene.dat -thread 4  2> resGene.log &


You can run these sub-pipeline on one computer (they have to listen on different port) or over a number of computer. You can even split a sub-pipeline to run over two computers. For example, you can run the gene resistance analysis on one computer::

   jsa.util.streamServer -port 3460 \ 
    | jsa.np.rtResistGenes -bam - -score=0.0001 -time 120 -read 50 --resDB  ResGene/resFinder/ -tmp _tmp_ -o resGene.dat -thread 4  2> resGene.log &

and run bwa on another::

   jsa.util.streamServer -port 3461 \ 
    | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a ResGene/resFinder/DB.fasta - 2> /dev/null \
    | jsa.util.streamClient -input - -server computer1:3460

which listens for streaming data in fastq format from port 3461, aligns to the resistance gene database, and forwards the alignments in sam format the resistance gene analysis via the network.

In these sub-pipelines, you may want to modify the parameter -port for  jsa.util.streamServer and -t for bwa to suit your computer systems.

Once these `daemons <https://en.wikipedia.org/wiki/Daemon_(computing)>`_ are ready for their analyses, you can start npReader to streamline data into the integrated pipeline::

   jsa.np.npreader -GUI -realtime -folder <DownloadFolder> -fail -output data.fastq -stream server1:port1,server2:port2,server3:port3
 
in which the -folder parameter specifies the downloads folder from the Metrichor base-calling, and the -stream parameter lists the computer addresses and port numbers that the analyses are listening on. At this point, you can start the MinION and Metrichor to start the real-time analyse.

=======================
Retro-realtime analysis
=======================

Note npreader no longer support extracting time information. Please use an earlier version for this.

If your data have been sequenced, and depending on what processing steps have been done.

* If your data have not been base-cased, you can start the pipeline as above, and run Metrichor for base-calling your   data.

* If your data have been base-called, and are still in fast5 format, you can run npReader as above to stream data to    the pipeline.

* If your data have been converted to fastq format, you can run jsa.util.streamClient to stream to the pipeline::

    jsa.util.streamClient -input reads.fastq -server server1:port1,server2:port2,server3:port3
  
* If you want to emulate the timing of your sequenced data, first convert the data to fastq format and extract the timing information (make sure parameter -time is turned on)::

   jsa.np.npreader -folder <downloads> -fail -number -stat -time -out dataT.fastq
  
Next sort the reads in the order they were generated::
  
   jsa.seq.sort -i dataT.fastq -o dataS.fastq --sortKey=timestamp
  
Finally, stream the data using jsa.np.timeEmulate::
  
   jsa.np.timeEmulate -input dataS.fastq -scale 1 -output - |jsa.util.streamClient -input - -server  server1:port1,server2:port2,server3:port3

You can crease the value in -scale to test higher throughput.

We provides the data from our four MinION runs in fastq format, sorted in the order
of sequencing (key=cTime). To re-run our analyses, set up the analysis pipeline as above,
and then stream our data through the pipeline, eg.,::

   wget http://data.genomicsresearch.org/Projects/npAnalysis/data.tar.gz
   tar zxvf data.tar.gz
   jsa.np.timeEmulate -input data/nGN_045_R7_X4S.fastq -scale 120 -output - |jsa.util.streamClient -input - -server  server1:port1,server2:port2,server3:port3

===================
Data from the study
===================

The sequencing data for the experiments in the paper have been deposited
to ENA under Accession `PRJEB14532 <http://www.ebi.ac.uk/ena/data/view/PRJEB14532>`_.


======================
Further documentations
======================

More details of usage of the discussed programs are provided in `ReadTheDocs for Japsa <http://japsa.readthedocs.org/en/latest/>`_. More specificially:

* `npReader <http://japsa.readthedocs.org/en/latest/tools/jsa.np.npreader.html>`_
* `jsa.util.streamServer <http://japsa.readthedocs.org/en/latest/tools/jsa.util.streamServer.html>`_
* `jsa.util.streamClient <http://japsa.readthedocs.org/en/latest/tools/jsa.util.streamClient.html>`_
* `jsa.np.filter <http://japsa.readthedocs.org/en/latest/tools/jsa.np.filter.html>`_
* `jsa.np.rtSpeciesTyping <http://japsa.readthedocs.org/en/latest/tools/jsa.np.rtSpeciesTyping.html>`_
* `jsa.np.rtStrainTyping <http://japsa.readthedocs.org/en/latest/tools/jsa.np.rtStrainTyping.html>`_
* `jsa.np.rtMLST <http://japsa.readthedocs.org/en/latest/tools/jsa.np.rtMLST.html>`_
* `jsa.np.rtResistGenes <http://japsa.readthedocs.org/en/latest/tools/jsa.np.rtResistGenes.html>`_

=======
Contact
=======
Minh Duc Cao -- m.cao1@uq.edu.au





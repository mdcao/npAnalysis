
#Real-time analysis of Oxford Nanopore MinION sequencing


This directory contains information for setting up real-time analysis 
of Oxford Nanopore sequencing data, as described in the under-review 
manuscript:

Real-time strain typing and analysis of antibiotic resistance potential using 
Nanopore MinION sequencing

Minh Duc Cao, Devika Ganesamoorthy, Alysha G. Elliott, Huihui Zhang, Matthew 
A. Cooper and Lachlan Coin


##Software installation

The streamline agorithms and auxiliary programs for setting up the pipeline 
are provided in japsa package (https://github.com/mdcao/japsa). In addition, the 
following free-software dependencies are required:

 1. Java >=1.7
 2. bwa >=7.10 (0.7.10-r858-dirty recommended)
 3. R with rJava and MultinomialCI installed
 4. kalign2 (http://msa.sbc.su.se/cgi-bin/msa.cgi)
 5. HDF5 (https://hdfgroup.org/HDF5/release/obtain5.html)
 -- optional, needed for translating native fast5 files to fastq

###Quick installation guide

```
git clone https://github.com/mdcao/japsa
cd japsa
make install \
   [INSTALL_DIR=~/.usr/local \]
   [MXMEM=7000m \]
   [SERVER=true \]
   [JLP=/usr/lib/jni:/usr/lib/R/site-library/rJava/jri]
```
 
This will install japsa according the directives:

* *INSTALL_DIR*: specifies the directory to install japsa, e.g, ~/.usr/local. Make sure
 you have write privilege to this folder. Also add INSTALL_DIR/bin (e.g., )
 to your PATH envirenment e.g., 
```
export PATH=~/.usr/local/bin:$PATH
```
* *MXMEM*: specifies the default memory allocated to the java virtual machine
* *SERVER*: specifies whether to launch the java virtual machine in server mode
* *JLP*: points to where HDF libraries and JRI are installed e.g,
 /usr/local/lib:/usr/lib/R/site-library/rJava/jri. The path to HDF is only needed for
 creating the pipelines to analyse directly from fast5 files or simultaneously
 from the MinION sequencing

For more detailed information for installing japsa, please refer to  Japsa installation guide on
http://japsa.readthedocs.org/en/latest/install.html

##Databases and data

Each of the analyses requires a pre-processed database. We make avaibale these 
databases on http://genomicsresearch.org/public/researcher/npAnalysis/
(and a backup storage 
https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/).
Setting up these databases for use as follows:

###Bacterial species database

We make available the database of all bacterial genomes obtained from NCBI genbank, 
with the addition of two K. quasipneumoniae strains (see the manuscript). Download
the database (~2.8GB), and make an bwa index of the database as follows.
```
wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/SpeciesTyping.tar.gz
tar zxvf SpeciesTyping.tar.gz
cd SpeciesTyping/Bacteria/
gunzip genomeDB.fasta.gz
bwa index genomeDB.fasta.gz
```
Note that it might take a while to build the bwa index for this 9G-base database.

###Strain typing with MLST
The database for MLST typing for three species,  K. pneumoniae, E. coli and
S. aureus are make avaibale. Download (208KB) and unzip them:
```
wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/MLST.tar.gz
tar zxvf MLST.tar.gz
```

###Strain typing with gene presence and absence
The database for gene presence and absence strain typing for K. pneumoniae, E. coli and
S. aureus can be obtained as follows:
```
wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/StrainTyping.tar.gz
tar zxvf ResGene.tar.gz
```

###Resistance gene identification
```
wget https://swift.rc.nectar.org.au:8888/v1/AUTH_15574c7fb24c44b3b34069185efba190/npAnalysis/ResGene.tar.gz
tar zxvf ResGene.tar.gz
```

##Setting up pipeline
The framework makes use of the interprocess communication mechanism pipe as well as network channels to set up the real-time pipeline. You can prepare one or more analyses to run in real-time.

For bacterial species typing,

```
jsa.util.streamServer -port 3456 \
  | bwa mem -t 10 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 <path>/SpeciesTyping/Bacteria/genomeDB.fasta - 2> /dev/null \
  | jsa.np.rtSpeciesTyping -bam - -index <path>/SpeciesTyping/Bacteria/speciesIndex --read 50 -time 60 -out speciesTypingResults.out 2>  speciesTypingResults.log 
```
This will create a pipeline to identify 
species which reports every 60 seconds, with at least 50 more reads 
from the last report. The pipeline waits for input on port 3456 for incoming data.


For strain typing gene presence/absense for K. pneumoniae  
```
jsa.util.streamServer -port 3457 \
  | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a ${base}/geneFam.fasta - 2> /dev/null \
  | jsa.np.rtStrainTyping -bam -  -geneDB <path>/StrainTyping/Klebsiella_pneumoniae/ -read 0 -time 20 --out kPStrainTyping.dat 2>  kPStrainTyping.log 
```
You can run strain typing pipelines for other species (e.g., E. coli and S. aureus) 
if you have reason to believe the sample may contain these species. If these pipeline
run on the same computer, make sure they listen to different ports.


```
jsa.util.streamServer -port 3458 \ 
  | bwa mem -t 8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y <path>/MLST/Klebsiella_pneumoniae/bwaIndex/genes.fasta - \
  |  jsa.np.rtMLST -bam - -mlst <path>/MLST/Klebsiella_pneumoniae/ -read 1000 -time 600  --out KpMLST.out
```

Again, you set up MLST for E. coli and/or S. aureus as well.

For resistance gene identification:
```
jsa.util.streamServer -port 3459 \ 
  | bwa mem -t 2 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -Y -K 10000 -a <path>/ResGene/resFinder/DB.fasta - 2> /dev/null \
  | jsa.np.rtResistGenes -bam - -score=0.0001 -time 120 -read 50 --resDB  <path>/ResGene/resFinder/  -tmp _tmp_ -o resGene.dat -thread 4  2> resGene.log 
```

You may want to modify the parameter -port for  jsa.util.streamServer 
and -t for bwa to suit your computer system. 

Once these programs are ready for their analysis, you can start npReader to 
streamline data into 

```
jsa.np.f5reader -GUI -realtime -folder <DownloadFolder> -fail -output data.fastq -stream server1:port1,server2:port2,server3:port3
```
in which the -folder parameter specifies the downloads folder from the Metrichor
base-calling, and the -stream parameter
lists the computer address and port number that the analyses are listening.

Once this is done, you can start the MinION and Metrichor to start the analyses.


In case the MinION run and base-calling have been done, you can also run
npReader to convert the data in fast5 format to fastq and streamline in
to the pipelines. If the data have also been convereted to fastq format,
you can run
```
jsa.util.streamClient -input <data.fastq> -server server1:port1,server2:port2,server3:port3
```

###Retro-realtime analysis

First, convert the data to fastq format and extract the timing information (make sure parameter -time is turned on):
```
jsa.np.f5reader -folder <downloads> -fail -number -stat -time -out dataT.fastq
```
Next sort the reads in the order they were generated:
```
jsa.seq.sort -i dataT.fastq -o dataS.fastq --sortKey=timestamp
```
Finally, stream the data using jsa.np.timeEmulate:
```
jsa.np.timeEmulate -input dataS.fastq  -scale 1  -output - |jsa.util.streamClient -input - -server  server1:port1,server2:port2,server3:port3
```
You can crease the value in -scale to test higher throughput.

We provides the data from our four MinION runs in fastq format, sorted in the order
of sequencing (key=cTime). To re-run our analyses, set up the analysis pipeline as above,
and then stream our data through the pipeline, eg.,
```
jsa.np.timeEmulate -input nGN_045_R7_X4S.fastq  -scale 100 -key=cTime -output - |jsa.util.streamClient -input - -server server1:port1,server2:port2,server3:port3
```





##Contact

Please contact Minh Duc Cao (minhduc.cao@gmail.com) for further information.


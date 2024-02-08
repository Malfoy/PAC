# PAC
PAC stand for Partitionned Aggregated bloom Comb trees.

PAC is a scalable structure able to list datasets that may contain a query sequence among thousands of indexed files.

PAC is an AMQ (Approximate membership Query) structure based on bloom filters like Cobs, or Sequence Bloom tree tools (HowDeSBT,AllsomeSBT,SSBT,...)

A PAC index can be built from thousands datasets quickly and with low ressource memory and disk footprint and can deliver very high query trhoughput.

# Table of contents
1. [Command line options](#cmo)
2. [Installation](#install)
3. [Index overview](#ov)
4. [Credits](#cr)


## Command line options <a name="cmo"></a>


### 1. Indexing  :

<ins>Default construction:

`./PAC -f test/fof.txt -d index_test`

PAC will construct an index in the folder index_test using all the datasets containend in the file of file fof.txt.
Each line will be considered as a separated file to index.
PAC can handle FASTA, MultiLine FASTA  and FASTQ files gziped or not.
You can easyly create a file of files using commands like 'readlink -f my_folder'


<ins>Construction with custom parameters:

`./PAC -f test/fof.txt -d index_test -k 21`

PAC will index 21mers

`./PAC -f test/fof.txt -d index_test -b 1000000000` 

Indexed Bloom filters will be sized to the first power of two above 1000000000 (1,073,741,824).
Larger  Bloom filters decrease false positive rate but increase the index size.


`./PAC -f test/fof.txt -d index_test -u` 

PAC will only index kmers seen at least two times in a given dataset (kmers seen only once are filtered)


`./PAC -f test/fof.txt -d index_test -p 5`

PAC will split the index into 4^p files inside the index_test folder.
Please check that you can open enough file with the ulimit command ('ulimit -n 5000' allow 5000 opened files)
More partitions usually speed up and reduce the memory usage during construction and query operations but create more files.

`./PAC -f test/fof.txt -d index_test -c 16`

PAC will use 16 threads to accelarate the index construction.
PAC is able to efficiently use a high amount of threads to divide the construction time but the memory usage can increase linearly with the amount of threads used.

`./PAC -f test/fof.txt -d index_test -e 8`

PAC will use 8 bits integers to encode its aggragated bloom filters, this way PAC can index up to 256 files.
Bydefault PAC use 16 bits integers and can index up to 65,636 files, if you need to index more files use -e 32



### 2. Query  :
<ins>Default query:

`./PAC -l index_test -q test/GCF_000864885.1_ViralProj15500_genomic.fna.gz`

PAC will load the index inside the index_test folder and query each line of the GCF_000864885.1_ViralProj15500_genomic.fna.gz fasta file.
The output will look like this
```
    >NC_004718.3 SARS coronavirus Tor2, complete genome
    29719 50 701
    
```
Meaning that the query share 29719 kmers with the first entry, 50 with the second and 701 with the third.

<ins>Multithread query:
  
`./PAC -l index_test -q test/GCF_000864885.1_ViralProj15500_genomic.fna.gz -c 16`
  
PAC will use 16 threads to perform the queries.


### 3. Update index :
Add a new file of file to a previously built index.

`./PAC -f test/fof_add.txt -l index_test -d index_updated  `


## Installation <a name="install"></a>

### Requirements


* A modern, C++17 ready compiler such as `g++` version 8 or higher.
* `zlib` to be already installed on your system (on Debian/Ubuntu it corresponds to the packages `zlib1g-dev`).

### Single user installation

To download and compile `PAC`  use the
following commands:

```sh
git clone --depth 1 --recursive https://github.com/Malfoy/PAC.git ;

cd PAC ;

make ;


```
## Index overview <a name="ov"></a>
TBD

## Credits <a name="cr"></a>

For further information our preprint "Scalable sequence database search using Partitioned Aggregated Bloom Comb-Trees"
 can be consulted in [Biorxiv](https://www.biorxiv.org/content/10.1101/2022.02.11.480089v2)



Authors:
----------------


* Camille MARCHET   <camille.marchet@univ-lille.fr>
* Antoine LIMASSET  <antoine.limasset@univ-lille.fr>


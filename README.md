# metaprob2


## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [General usage](#general)
  - [Getting help](#help)
  - [Citing metaprob2](#cite)
- [Developers' Guide](#dguide)
- [Limitations](#limit)

## <a name="uguide"></a>Users' Guide
Sequencing technologies allow the sequencing of microbial communities directly from the environment without prior culturing. One of the major problems when analyzing a microbial sample is to taxonom- ically annotate its reads to identify the species it contains. Taxonomic analysis of microbial communities requires reads clustering, a process referred to as binning. The major problems of metagenomics reads bin- ning are the lack of taxonomically related genomes in existing reference databases, the uneven abundance ratio of species, and sequencing errors. 
In this paper we present MetaProb 2 an unsupervised binning method based on reads assembly and probabilistic k-mers statistics. The novelties of MetaProb 2 are the use of minimizers to efficiently assemble reads into unitigs and a community detection algorithm based on graph modularity to cluster unitigs and to detect representative unitigs. The effectiveness of MetaProb 2 is demonstrated in both simulated and synthetic datasets in comparison with state-of-art binning tools such as MetaProb, Abun- danceBin, Bimeta and MetaCluster.
---

## <a name="install"></a>Installation

In order to work, MetaProb2 needs 3 pieces of sotfware:

1. Minimap2, available  at https://github.com/lh3/minimap2; 
You can also use: git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)
2. Miniasm, available at https://github.com/lh3/miniasm;
You can also use: git clone https://github.com/lh3/miniasm  && (cd miniasm  && make)
3. MetaProb, available at https://bitbucket.org/samu661/metaprob/src/master/;
Once downloaded, go to  MetaProb/Release/ and then use: make all.

You need gcc and zlib to install Minimap2 and Miniasm; You also need Boost (https://www.boost.org/users/download/) and Eingen (http://eigen.tuxfamily.org/index.php?title=Main_Page) libraries in the local directory to use MetaProb. 
You need to install scikit-network for python3 (https://scikit-network.readthedocs.io/en/latest/).
Please follow the guides provided at the links above to correctly install the tools.
---

##  <a name="general"></a>General Usage

MetaProb2 is a metagemomic binning tool that uses mapping and assembly software together with a novel metagenomic community detection script to improve the results of MetaProb https://academic.oup.com/bioinformatics/article/32/17/i567/2450796 , https://bitbucket.org/samu661/metaprob/src/master/. 
We hereby provide the python3 and shell scripts to make the pipeline work without much effort.
The usage of the 3 different softwares is well explained at the liks provided above. For non-custom usage, you can simply downlaod the softwares, the libraries, the 2 python3 scripts provided in this repository and the shell script.
Set the parameters in the script: directories, files name and number of expected cluster. You can choose to provide no expected number of clusters: if so, please remove "-numSp $num" from MetaProb input.


---

## <a name="help"></a>Getting help


---

## <a name="cite"></a>Citing metaprob2

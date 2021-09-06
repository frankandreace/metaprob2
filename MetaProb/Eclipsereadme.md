#MetaProb

##Abstract
###Motivation
Sequencing technologies allow the sequencing of microbial communities
directly from the environment without prior culturing. Taxonomic analysis of microbial communities,
a process referred to as binning, is one of the most challenging task when analyzing metagenomic
reads data. The major problems are the lack of taxonomically related genomes in existing reference
databases, the uneven abundance ratio of species, and the limitations due to short read lengths and
sequencing errors.

###Results
MetaProb is a novel assembly-assisted tool for un-supervised metagenomic binning. The
novelty of MetaProb derives from solving a few important problems: how to divide reads into groups
of independent reads, so that l-mer frequencies are not overestimated; how to convert l-mer counts
into probabilistic sequence signatures, that will correct for variable distribution of l-mers, and for
unbalanced groups of reads, in order to produce better estimates of the underlying genome statistic.
We show that MetaProb is effective for both simulated and real datasets. It can accurately (with
F-measures of 87 for short reads and 97 long reads) and efficiently bin short and long reads with
varying abundance ratios.

#Installation

##Install dependencies
Install in local directory the following library:  
  
[Boost library](http://www.boost.org/users/download/)  
[Eingen library](http://eigen.tuxfamily.org/index.php?title=Main_Page)  

##Option Project Eclipse
Link library:  
+ GCC C++ Linker -> Libraries -> -l pthread boost_system  
+ GCC C++ Linker -> miscellaneus -> -fopenmp  
  
CPP compiler options: -std=c++0x -O0 -floop-parallelize-all -g3 -Wall -c -fmessage-length=0 -fopenmp  
Dialect: -std=c++0x  
Optimization: -O3 -floop-parallelize-all  

##Download Metaprob
Download MetaProb at:  
[MetaProb v2](https://bitbucket.org/samu661/metaprob/downloads/MetaProb_v2.tar.gz)  

##Compilation:
Open terminal and go to MetaProb/Release/ and then:  
make all  

###Option compilation
Y = Number of core  
-jY

#Algorithm Option
In the following paragraph is described the input file's structure and the parameters available in the MetaProb algorithm.  

##File accepted and structure
File accepted have the following structure:  
Structure file .fna example:  
> \>IDENTIFICATION  
> ATAATTGGCAAGTGTTTTAGTCTTAGAGAGATTCTCTAAGTCTAACTTGAACTCAATTTGGAATCATTTCCCAATTTTTA

Structure .fastq exemple:  
> @IDENTIFICATION #0/1  
> CCCATGCCTTTAGCCAAATTCACGGTTTGATCACCCCTAAAACCAGCCAATATACCGAAGTGGAAGCCAGCATAAATGGCCTCAATATTACCGAAATGGAT  
> +  
> HBIIIIIIIHHDIHIGIIGGIHIIGIDIIIIBIHI@IIH@HIIHIIF5IIHEII>BDAHIBIEDBEIDG@HAEH*I@AEI=#CE?G17EEDHDEB@@?#8B  
  
In this NEW VERSION, paired-end reads are passed to the algorithm in separeted file in which the reads are paired 
in the same order in which are writen. So we raccomend to control the paired-end read if they are paired in the
correct manner.

##Parameter
**-si** File path single-end reads  
**-pi** File paths paired-end reads   
**-dirOutput** Path directory output files. Default: output/  
**-graphType** 0 = pair read in Paired-End read with same ID, 1 = All read Single-End, 2 = All read Single-End and then union groups with Paired-End info. Default: 0  
**-numSp** Number expected specie in file. Default: 2  
**-q** Size of q-mer used to create graph adiacences: Default: 30  
**-m** Threshold of shared q-mer to create graph adiacences. Default: 5  
**-ssize** Max Seed size in each group. Default: 9000  
**-lmerFreq** Size of L-mer used to compute feature vector. Default: 4  
**-feature** Feature used to compute. Default: 1  
**-mg** Only group mode activated (output only the groups created). Default: Not Active  
**-eK** Estimate K for K-means algorithm. Default: Active unless you specify -numSp  

##Advanced Option
**-iterMaxKmeans** Max iteration of Kmeans algorithm. Default: 100  
**-timeMaxKmeans** Max time in seconds of Kmeans algorithm. Default: 3600  

#Feature available
There are several features available to describe the information contained in the groups, but two of this are the best:    
+ *For Short Read (bp < 300)*	-> **-feature 1** = NORM_D2star_All_Read_Prob_Lmer_Euclidian  
+ *For Long Read (bp > 300)*	-> **-feature 2** = NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian  
  
**-feature 1** = NORM_D2star_All_Read_Prob_Lmer_Euclidian, //NORM_D2star_All_Read_Prob_Kmer + Euclidian norm  
**-feature 2** = NORM_D2star_Group_Prob_Bernulli_Lmer_Euclidian, //NORM_D2star_Group_Prob_Bernulli_Kmer + Euclidian norm  
**-feature 3** = NORM_BiMeta, //Bimeta's paper norm  
**-feature 4** = NORM_D2star_Group_Prob_Lmer, //probability of L-mer in a group, given seed read L-mer's count vector  
**-feature 5** = NORM_D2star_All_Read_Prob_Lmer, //probability of L-mer in collection, given seed read L-mer's count vectors  
**-feature 6** = NORM_D2star_Group_Prob_Bernulli_Lmer, //probability of L-mer with Bernulli model in a group  
**-feature 7** = NORM_D2star_All_Read_Prob_Bernoulli, //probability of L-mer with Bernulli model in collection  
**-feature 8** = NORM_D2star_All_Seed_Read_Prob_Bernoulli, //probability of L-mer with Bernulli model in collection seed read  
**-feature 9** = NORM_BiMeta_Euclidian, //NORM_Size_Seed + Euclidian norm  
**-feature 10** = NORM_D2star_Group_Prob_Lmer_Euclidian, //NORM_D2star_Group_Prob_Kmer + Euclidian norm  
**-feature 11** = NORM_D2star_All_Read_Prob_Bernoulli_Euclidian, //NORM_D2star_Prob_Bernoulli_All_Read + Euclidian norm  
**-feature 12** = NORM_D2star_All_Seed_Read_Prob_Bernoulli_Euclidian, //NORM_D2star_Prob_Bernoulli_All_Seed_Read + Euclidian norm  

#Run
Calls algorithm where is compiled:  
./MetaProb -si ../TestInputFile/long_example_1.fna -numSp 2 -feature 10 -m 45  
./MetaProb -pi ../TestInputFile/short_example_1.fna.1 ../TestInputFile/short_example_1.fna.2 -numSp 2  
./MetaProb -pi ../TestInputFile/short_example_2.fna.1 ../TestInputFile/short_example_2.fna.2 -numSp 2 -feature 9  

##Best Parameter
+ *For Short Read (bp < 300)*	-> **-feature 1 -m 5**  
+ *For Long Read (bp > 300)*	-> **-feature 2 -m 45**  
  
#Test File
##Single-End Reads:
[Single-End_dataset_part1](https://bitbucket.org/samu661/metaprob/downloads/single_end_dataset.part1.rar)  
[Single-End_dataset_part2](https://bitbucket.org/samu661/metaprob/downloads/single_end_dataset.part2.rar)  

##Paired-End Reads split file
[Paired-End_split_dataset_part01](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part01.rar)  
[Paired-End_split_dataset_part02](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part02.rar)  
[Paired-End_split_dataset_part03](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part03.rar)  
[Paired-End_split_dataset_part04](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part04.rar)  
[Paired-End_split_dataset_part05](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part05.rar)  
[Paired-End_split_dataset_part06](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part06.rar)  
[Paired-End_split_dataset_part07](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part07.rar)  
[Paired-End_split_dataset_part08](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part08.rar)  
[Paired-End_split_dataset_part09](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part09.rar)  
[Paired-End_split_dataset_part10](https://bitbucket.org/samu661/metaprob/downloads/paired_end_dataset_splitted.part10.rar)  
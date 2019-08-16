# Bit-mapping
This is an experimental tool for mapping RNA-seq reads to reference transcriptome implemented by C++ in Linux.

## Introduction

Bit-mapping accelerates mapping speed by reducing the dimension of the reads. It uses learning to hash algorithm, spherical hashing, to compute the binary hash codes of reads, and map reads to reference transcriptome by their hash codes rather than the original data.

## Installation

Just clone or download this respository, and compile it by

> bash compile.sh

## Dependencies

- [Cereal](https://github.com/USCiLab/cereal) 

- [libdivsufsort](https://github.com/y-256/libdivsufsort)

- openmp is supported by compiler

Please ensure that they are installed in your machine before installing bit-mapping.

## Usage

There are several important parameters in bit-mapping:

|Parameter|Function|Default value|
|---|---|---|
|-k|set the length of k-mer|3|
|-r|set the file name of the reference transcriptome | |
|-s|set the length of the skip part|6|
|-m|set the length of the region-searching part|10|
|-g|set the length of the ignored part|0|
|-p|set thread|4|
|-I|set the first mates of the paired-end reads| |
|-i|set the second mates of the paired-end reads| |
|-o|set the file name of the output sam file| |
|-c|set the length of the hash code|37|
|-t|set the tolerance of similarity|2|


To build the index of the reference transcriptome, please run the following command:
> bmap index -r ref.fa -c max_code_len

The max_code_len is the maximum length of hash code length of reference transcriptome, which is specified by users. Its value is 128 bits in our example index data. When building the index of the reference, users have to specify the maximum length of the hash code. Once the index is produced, it can be used to other datasets with the same settings. Note that the given length of the hash codes cannot exceed that of the index, or you should re-produce the index for longer length. Please make sure the index folder is under your current working directory.

To run the mapping phase, run the command as follows:

>bmap mapping -c code_len -t tolerance -p thread_num -I read_1.fastq -i read_2.fastq -o res.sam

This is the command for paired-end reads. read_1.fastq is the first mates of the paired-end reads, read_2.fastq is the second mates paired-end reads, and res.sam is the output sam file. code_len is the specified length of hash code of reads. For n bp reads, the maximum value of code_len is:

> n - length of region-searching part - length of skipping part - length of ignored part - length of k-mer + 1

## Simulated data

The simulated data with 250 bp reads and the index of the reference transcriptome are contained in the fold named sample_data. The maximum length of hash code is 128 in this index.
The reference transcriptome used in this example is [homo_sapiens](ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz) in Ensembl.


## Note

We only focus on optimising the mapping speed, since the speed of data reading and writing is effected by the manner of implementation. For example:

> hash code length:107  
  read dim:125  
  total reads:2000000  
  Parse Fastq File Finished(6.31085 seconds)  
  Loading Spherical Hashing Finished (0.000178 seconds)  
  Total Hashing Time: (13.544875 seconds)  
 [Building BooPHF]  100  %   elapsed:   0 min 1  sec   remaining:   0 min 0  sec  
  Computing Hashcode Of Ref Kmer Using BBhash Finished (1.118350 seconds)  
  Analysis of read region ...  
  Analysis Finished(1.25309 seconds)  
  Load Suffix Array Finished (5.880918 seconds)  
  Load Suffix Array Region Finished (8.534228 seconds)  
  Loading Reference Code Time: (22.036242 seconds)  
  Resizing Ref Hashcode Finished: (176.343352 seconds)  
  Total Mapping Time: (22.877301 seconds)  
  Load Ref Info Finished (9.036697 seconds)  
  Load Code Bucket Info Finished (21.639003 seconds)  
  Analyse Results Finished (3.275507 seconds)  
  Merge Result Finished (2.416870 seconds)  
  both  hamming distances are larger than 10:4313  
  half mapping:70315  
  total reads:2000000  
  Parse Fastq File Finished(10.7545 seconds)  
  Generate SAM File Finished (42.518203 seconds)  
  Total Running Time (355.328327 seconds)  

This is an example log when running bit-mapping. The crucial time of mapping process includes the total hashing time (line 6) and the total mapping time (line 15). You can get the time of mapping process by summing the hashing time and the mapping time.

## Authors

The bit-mapping algorithm is developed by Xueting Yu and Xuejun Liu. The software is mainly implemented by Xueting Yu.

## Contact information

For any query, please contact Xueting yu via xt_yu@nuaa.edu.cn or Xuejun Liu via xuejun.liu@nuaa.edu.cn.
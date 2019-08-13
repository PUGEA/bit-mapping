#pragma once

#ifndef _FASTX_PARSER_
#define _FASTX_PARSER_

#include <fstream>
#include "Points.h"

class FastxParser{

private:
	std::ifstream seqfile__;
	std::ifstream seqfile2__;
	/*std::ofstream out_ref_file__;
	std::ofstream out_read_file_1__;
	std::ofstream out_read_file_2__;

	int klen__;
	int kmer_count__;*/
	int dim;

	enum filetype__{DONE,FASTA, FASTQ,ERROR};
	enum readtype__{SINGLE,PAIR};

public:
	FastxParser();
	FastxParser(int dim);
	void initialize(int dim);
	// Adapted from
	// https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
	// Don't modify the qual
	void reverse_complete(std::string& seq,std::string& rev_Read );

	void reverse_complete_ictoic(std::string& seq,std::string& rev_Read);

	void reverse_complete_ictos(REAL_TYPE *seq,std::string& rev_Read);

	void reverse_complete(REAL_TYPE* seq,REAL_TYPE* rev_Read);

	filetype__ get_file_type(std::ifstream &file);

	void open_seqfile(std::string filename);

	void open_seqfile(std::string filename1,std::string filename2);

	void read_fq_oneseq(std::ifstream &seqfile,std::string &name, std::string &seq, std::string &qual);

	void  read_fa_oneseq(std::ifstream &seqfile,std::string &name, std::string &seq);


	std::string merge_ref_seq(std::string filename,int seg_len);

	//directly store the reads to the files
	size_t store_reads(Points &read_1_buff,Points &read_2_buff,std::string read_file_1,std::string read_file_2);
};

#endif
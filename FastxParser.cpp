#include <iostream>
#include <fstream>
#include <sstream>//stringstream
#include <vector>
#include <limits>
#include <cstdint>
#include <chrono>
#include <cereal/archives/binary.hpp>
#include "Stopwatch.hpp"

#include "Utils.hpp"
#include "Points.h"
#include "FastxParser.hpp"

FastxParser::FastxParser(){}
FastxParser::FastxParser(int dim)
{
	this->dim = dim;
}
void FastxParser::initialize(int dim){
	this->dim = dim;
}
// Adapted from
// https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
// Don't modify the qual
void FastxParser::reverse_complete(std::string& seq,std::string& rev_Read ) {

    rev_Read.resize(seq.length(), 'A');
    int32_t end = seq.length()-1, start = 0;
    //rev_Read[end] = '\0';
    //qualWork[end] = '\0';
    while (start < end) {
        rev_Read[start] = (char)rc_table[(int8_t)seq[end]];
        rev_Read[end] = (char)rc_table[(int8_t)seq[start]];
        ++ start;
        -- end;
    }
    // If odd # of bases, we still have to complement the middle
    if (start == end) {
        rev_Read[start] = (char)rc_table[(int8_t)seq[start]];
        // but don't need to mess with quality
        // qualWork[start] = qual[start];
    }
    //std::swap(seq, rev_Read);
    //std::swap(qual, qualWork);
}

void FastxParser::reverse_complete_ictoic(std::string& seq,std::string& rev_Read) {
    rev_Read.resize(seq.length(), 'A');
    int32_t end = seq.length()-1, start = 0;
    while (start < end) {
        rev_Read[start] = rc_ictoic_table[(int8_t)seq[end]];
        rev_Read[end] = rc_ictoic_table[(int8_t)seq[start]];
        ++ start;
        -- end;
    }
    // If odd # of bases, we still have to complement the middle
    if (start == end) {
        rev_Read[start] = rc_itoi_table[(int8_t)seq[start]];
        // but don't need to mess with quality
    }
}

void FastxParser::reverse_complete_ictos(REAL_TYPE *seq,std::string& rev_Read) {
    rev_Read.resize(dim, 'A');
    int32_t end = dim - 1, start = 0;
    while (start < end) {
        rev_Read[start] = (char)rc_ictos_table[(int8_t)seq[end]];
        rev_Read[end] = (char)rc_ictos_table[(int8_t)seq[start]];
        ++ start;
        -- end;
    }
    // If odd # of bases, we still have to complement the middle
    if (start == end) {
        rev_Read[start] = rc_ictos_table[(int8_t)seq[start]];
    }
}

void FastxParser::reverse_complete(REAL_TYPE* seq,REAL_TYPE* rev_Read) {
    int32_t end = dim - 1, start = 0;
    while (start < end) {
        rev_Read[start] = rc_itoi_table[(int8_t)seq[end]];
        rev_Read[end] = rc_itoi_table[(int8_t)seq[start]];
        ++ start;
        -- end;
    }
    // If odd # of bases, we still have to complement the middle
    if (start == end) {
        rev_Read[start] = rc_itoi_table[(int8_t)seq[start]];
    }
}

FastxParser::filetype__ FastxParser::get_file_type(std::ifstream &file){
	
	switch(file.peek()) {
	    case EOF: return DONE;
	    case '>': return FASTA;
	    case '@': return FASTQ;
	    default: return ERROR;
	}
	
}

void FastxParser::open_seqfile(std::string filename){
	seqfile__.open(filename,std::ifstream::binary);//|std::fstream::in);
	if(!seqfile__.is_open()){
		std::cerr<<"fail to open fasta/fastq files:" << filename << std::endl;
		exit(0);
	}
}

void FastxParser::open_seqfile(std::string filename1,std::string filename2){
	seqfile__.open(filename1,std::ifstream::binary);//,std::fstream::in);
	seqfile2__.open(filename2,std::ifstream::binary);//,std::fstream::in);
	if(!seqfile__.is_open()){
		std::cerr<<"fail to open fasta/fastq files:" << filename1 << std::endl;
		exit(0);
	}
	if(!seqfile2__.is_open()){
		std::cerr<<"fail to open fasta/fastq files:" << filename2 << std::endl;
		exit(0);
	}
}

void FastxParser::read_fq_oneseq(std::ifstream &seqfile,std::string &name, std::string &seq, std::string &qual){
	std::string line;
//	std::string seq;
	//ignore @
	char tmp;

	seqfile.get(tmp);
	std::getline(seqfile,line);
	//use sstream to split the name and description
	std::stringstream ss(line);
	name.clear();
	ss >> name;

	while(seqfile.peek() != '+' && seqfile.peek() != EOF){	
		std::getline(seqfile,line);
		seq.append(line);
	}

	//ignore the description
	seqfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	qual.clear();

	while(qual.size() < seq.length() && seqfile.good()) {
	      std::getline(seqfile, line);
	      qual.append(line);
	}
}

void FastxParser::read_fa_oneseq(std::ifstream &seqfile,std::string &name, std::string &seq){
	std::string line;
//	std::string seq;
	int len;
	//ignore >
	seqfile.get();
	std::getline(seqfile,line);
	std::stringstream ss(line);
	name.clear();
	ss >> name;

//	seqint.clear();
	while(seqfile.peek() != '>' && seqfile.peek() != EOF){
		std::getline(seqfile,line);
		seq.append(line);
	} 
//	trans_seq2array(seq,seqint);
}


std::string FastxParser::merge_ref_seq(std::string filename,int seg_len){
	open_seqfile(filename);
	filetype__ filetype = get_file_type(seqfile__);

	size_t size;
	std::string name;
	std::string qual;
	std::string seq,all_seq;
	int len = 0;

	std::ofstream train_data(INPUT_REF_FILE_NAME);
	int train_num = 0;
	int start = Parameter::region_searching;
	int cond = dim - Parameter::ignore - seg_len + 1; //BCODE_LEN * seg_len + Parameter::region_searching;

	std::vector<int> loc_to_ref;
	std::vector<std::string> ref_name;
	std::vector<size_t> ref_start;

	srand (time(NULL));
	auto t1 = std::chrono::high_resolution_clock::now();
	for (size  = 0;seqfile__.peek() != EOF;size++){
		seq.clear();
		read_fa_oneseq(seqfile__,name,seq);
	//	std::cout << seq << std::endl;
		// generate training data
		if(seq.length() > dim)
		{
			if (train_num < NUM_TRAIN_SAMPLES)
			{
				train_num++;
				int rand_loc = rand() % (seq.length() - dim + 1);
				std::string train_read = seq.substr(rand_loc,dim);
				
				while(true)
				{
					if (start < cond)
					{
						for (int j = 0; j < seg_len; ++j)
						{
							train_data << stoic_table[(int8_t)train_read[start + j]];
						//	std::cout << d[k][j];
						}
						train_data << std::endl;
						//std::cout << std::endl;
						//start += seg_len;
						start += 1;
					}else
					{
						start = Parameter::region_searching;
						break;
					}		
				}
			}

		}
		ref_name.push_back(name);
		ref_start.push_back(len);
		len += seq.size() + 1;
		loc_to_ref.resize(len,size);
		for (int i = 0; i < seq.size(); ++i)
		{
			all_seq.push_back((char)stoic_table[(int8_t)seq[i]]);
		}
		all_seq.append("9");
		
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	double fastxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	std::cout << "- Total Transcriptomes:" << size << std::endl;
	std::cout << "- The Whole Transcriptomes' Length:" << len << std::endl;
	std::cout << "- Time Of Merging Ref Fasta (" << fastxparser_time << " seconds)" << std::endl;

	// save loc_to_ref info to disk
	{
		std::ofstream loc_to_ref_file(MERGE_REF_POS_FILE);
		cereal::BinaryOutputArchive ar(loc_to_ref_file);
		ar(CEREAL_NVP(loc_to_ref));
	}

	// save merged ref to disk
	{
		std::ofstream merged_ref_file(MERGE_REF_SEQ_FILE);
		cereal::BinaryOutputArchive ar(merged_ref_file);
		ar(CEREAL_NVP(all_seq));
	}

	{
		std::ofstream rname_file(REF_NAME_FILE);
		cereal::BinaryOutputArchive ar(rname_file);
		ar(CEREAL_NVP(ref_name));
	}
	{
		std::ofstream rstart_file(MERGE_REF_START_FILE);
		cereal::BinaryOutputArchive ar(rstart_file);
		ar(CEREAL_NVP(ref_start));
	}

	t2 = std::chrono::high_resolution_clock::now();
	fastxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	std::cout << "- Time Of Saving Reference Position In Merged Ref:" << fastxparser_time << std::endl;
	seqfile__.close();
	train_data.close();
	return all_seq;
}

//directly store the reads to the files
size_t FastxParser::store_reads(Points &read_1_buff,Points &read_2_buff,std::string read_file_1,std::string read_file_2){
	system("mkdir tmp");
	open_seqfile(read_file_1,read_file_2);

	size_t size  = 0;
	std::vector<std::string> name_vec;
	std::string name;
	std::string fqual;
	std::string squal;
	std::string first;
	std::string second;

	std::pair<std::string,std::string> seqpair;
	std::pair<std::string,std::string> qualpair;

	filetype__ filetype = get_file_type(seqfile__);
	filetype__ filetype2 = get_file_type(seqfile2__);
	std::string rc_read;
	
	Stopwatch T0("");
    T0.Reset();     T0.Start();
	if (filetype == FASTQ && filetype2 == FASTQ){
		// count size of reads
		for (;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;size++){
			std::getline(seqfile__,first);
		}
		size /= 4;
		std::cout << "- total reads:" << size << std::endl;
		seqfile__.seekg(0,std::ios::beg);
		read_1_buff.Initialize(size,dim);
		read_2_buff.Initialize(size,dim);
		
		for (int i = 0;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;i++){
			first.clear();
			second.clear();
			name.clear();
			read_fq_oneseq(seqfile__,name,first,fqual);
			read_fq_oneseq(seqfile2__,name,second,squal);
			name_vec.push_back(name);
			for (int j = 0; j < dim; ++j){
				read_1_buff.d[i][j] = stoi_table[(int8_t)first[j]];
				read_2_buff.d[i][j] = stoi_table[(int8_t)second[j]];
			}
		}
	}
	T0.Stop();
	std::cout << "- Parse Fastq File Finished(" << T0.GetTime() << " seconds)" << std::endl;
	
	seqfile__.close();
	seqfile2__.close();

	//store reads name
	{
		std::ofstream reads_name(PAIR_1_NAME_FILE);
		cereal::BinaryOutputArchive ar(reads_name);
		ar(name_vec);
	}
	return size;
}




/*

void get_one_kmer(std::string seq,int loc, bool if_out){
	std::string kmer;
	int ik = 0;
	kmer.resize(klen__);
	for (int i = loc; i < loc + klen__; ++i){
		kmer[ik] = seq[i];
		ik++;
		if (if_out){
			out_ref_file__ << (char)stoic_table[(int8_t)seq[i]];
		}
		
	}
	if (if_out){
		out_ref_file__ << std::endl;	
	}
	
	kmer_count__++;																																																																																																																																																																																																																																																										
}


// directly store the ref kmers to the file
void store_ref_kmers(std::string filename, int klen){
	open_seqfile(filename);
	filetype__ filetype = get_file_type(seqfile__);

	size_t size;
	std::string name;
	std::string qual;
	std::string seq;
	int len;

	klen__ = klen;
	kmer_count__ = 0;

	out_ref_file__.open("data.txt",std::ofstream::binary|std::ofstream::out);
	if (!out_ref_file__.is_open())
	{
		std::cerr << "can't open out_ref_file: data.txt" << std::endl;
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	for (size  = 0;seqfile__.peek() != EOF;size++){
//	for (size = 0;size < TEST_REF;size++){
		read_fa_oneseq(seqfile__,name,seq);
		std::cout << "parse transcriptome:" << name << std::endl;
		len = seq.size();
		for (int loc = 0;loc < (len - klen + 1);loc++){	
			get_one_kmer(seq,loc,true);
		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	double fasxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	std::cout << "total transcriptomes:" << size << std::endl;
	std::cout << "total kmers:" << kmer_count__ << std::endl;
	std::cout << "parser time of fasta:" << fasxparser_time << std::endl;

	seqfile__.close();
	out_ref_file__.close();
}



//directly store the reads to the files
void store_reads(char* filename,int klen){
	open_seqfile(filename);

	klen__ = klen;

	out_read_file_1__.open("dataset/ref_sim.txt",std::ofstream::out);
	if (!out_read_file_1__.is_open()){
		std::cerr << "can't open out_read_file: read_1.txt" << std::endl;
	}

	size_t size  = 0;
	std::vector<std::string> name_vec_1;
	std::string name_1;
	std::string fqual;
	std::string squal;
	std::string first;
	std::pair<std::string,std::string> seqpair;
	std::pair<std::string,std::string> qualpair;

	filetype__ filetype = get_file_type(seqfile__);
	if (filetype == FASTQ){
		for (;seqfile__.peek() != EOF;size++){
	//	for (;size < TEST_READ;size++){
			first.clear();
			name_1.clear();
			read_fq_oneseq(seqfile__,name_1,first,fqual);

			name_vec_1.push_back(name_1);

			for (int i = 0; i < klen__; ++i){
				out_read_file_1__ << (char)stoic_table[(int8_t)first[i]];
			}
			out_read_file_1__ << std::endl;
		}
	}
	std::cout << "- total reads:" << size << std::endl;
	seqfile__.close();
	out_read_file_1__.close();

	//store reads name
	{
		std::ofstream reads_name_1(PAIR_1_NAME_FILE);
		cereal::BinaryOutputArchive ar(reads_name_1);
		ar(name_vec_1);
	}
}

//store data to RAM
void read_fastx(char* filename,singleSeqList *seqlist){	
	open_seqfile(filename);
	filetype__ filetype = get_file_type(seqfile__);

	size_t size;
	std::string name;
	std::string qual;
	std::string seq;
	if (filetype == FASTQ){

		for (;seqfile__.peek() != EOF;size++){
	//	for(size = 0;size < TEST_READ;size++){
			read_fq_oneseq(seqfile__,name,seq,qual);
			seqlist->seq.push_back(std::move(seq));
			seqlist->name.push_back(std::move(name));
			seqlist->qual.push_back(std::move(qual));
		}
		seqfile__.close();
	}else if(filetype == FASTA){
		auto t1 = std::chrono::high_resolution_clock::now();
		for (size  = 0;seqfile__.peek() != EOF;size++){
	//	for (size = 0;size < TEST_REF;size++){
			read_fa_oneseq(seqfile__,name,seq);
			seqlist->seq.push_back(std::move(seq));
			seqlist->name.push_back(std::move(name));
		}
		auto t2 = std::chrono::high_resolution_clock::now();
		double fasxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
		std::cout << "Fastx Parser Time:" << fasxparser_time << std::endl;
	}
	seqfile__.close();
}

//store data to RAM
void read_fastx(char* filename1, char* filename2,pairSeqList *seqlist){	
	open_seqfile(filename1,filename2);

	size_t size  = 0;
	std::string name;
	std::string fqual;
	std::string squal;
	std::string first;
	std::string second;
	std::pair<std::string,std::string> seqpair;
	std::pair<std::string,std::string> qualpair;

	filetype__ filetype = get_file_type(seqfile__);
	filetype__ filetype2 = get_file_type(seqfile2__);
	if (filetype == FASTQ && filetype2 == FASTQ){
		for (;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;size++){
	//	for (;size < TEST_READ;size++){
			read_fq_oneseq(seqfile__,name,first,fqual);
			read_fq_oneseq(seqfile2__,name,second,squal);
			seqpair.first = first;
			seqpair.second = second;
			qualpair.first = fqual;
			qualpair.second = squal;
			seqlist->seq.push_back(std::move(seqpair));
			seqlist->name.push_back(std::move(name));
			seqlist->qual.push_back(std::move(qualpair));

		}
	}else if(filetype == FASTA && filetype2 == FASTA){
		for (size  = 0;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;size++){
	//	for (size = 0;size < TEST_REF;size++){
			read_fa_oneseq(seqfile__,name,seqpair.first);
			read_fa_oneseq(seqfile2__,name,seqpair.second);
			seqlist->name.push_back(std::move(name));
			seqlist->seq.push_back(std::move(seqpair));
		}
	}else{
		std::cerr << "the types of two files are not the same. " << std::endl;
		seqfile__.close();
		seqfile2__.close();
		return ;
	}
	seqfile__.close();
	seqfile2__.close();
}*/
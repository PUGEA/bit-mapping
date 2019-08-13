#ifndef _KMER_UTILS_
#define _KMER_UTILS_

#include<stdio.h>
#include<iostream>
#include<vector>
#include<cstdint>
#include <falconn/lsh_nn_table.h>
#include <chrono>
#include"fastxParser.hpp"
#include"matio.h" 


//#define TRANSCRIPT_LIST std::vector<seqInfo>
//#define KMER_LIST std::vector<refKmerInfo>
#define GET_ARRAY_LEN(array,len) {len = (sizeof(array) / sizeof(array[0]));}

struct refKmerList{	
	std::vector<point> seq;
	std::vector<std::string > tname;
	std::vector<int> location;
};


class kmerUtils{
private:
	//----------
	std::ofstream outfile;
	//----------

	refKmerList klist;
	
	size_t klen;

	int whether_kmer_exist(point new_kmer){
		int index = 0;
		int i = 1;
		for (auto &kmer : klist.seq){
			for (i = 0; i < klen; ++i){
				if (kmer[i] != new_kmer[i]){
					break;
				}
			}
			if (i == klen - 1){
				return index;
			}
			index++;
		}
		return -1;
	}

	point get_one_kmer(point seq,int loc){
		point kmer;
		int ik = 0;
		kmer.resize(klen);
		for (int i = loc; i < loc + klen; ++i){
			kmer[ik] = seq[i];
			ik++;
		}
		return kmer;																																																																																																																																																																																																																																																												
	}

	void genkmer_from_one_transcript(point seq,std::string name){	
		point kmer;
		int len = seq.size();
		std::string tname;
		int location;

	//	int count_same_kmer[len - klen + 1] = {0};

		for (int loc = 0;loc < (len - klen + 1);loc++){	
			kmer = get_one_kmer(seq,loc);
		//	int index = whether_kmer_exist(kmer);
			
			klist.seq.push_back(std::move(kmer));
			//tname.push_back(std::move(name));
			klist.tname.push_back(std::move(name));
			//location.push_back(std::move(loc));
			klist.location.push_back(std::move(loc));
		/*	if (index == -1){
				
			}else{
			//	klist.tname[index].push_back(std::move(name));
			//	klist.location[index].push_back(std::move(loc));
			//	count_same_kmer[index]++;
			}*/

	/*		for (int i = 0; i < kmer.size(); ++i)
			{
				std::cout << kmer[i];
			}
			std::cout << std::endl;*/

			
		}
	/*	std::cout << "same kmer count:" << std::endl;
		for (int i = 0; i < len-klen+1; ++i)
		{
			std::cout << count_same_kmer[i] << " ";
		}
		std::cout << std::endl;*/
	}

public:

	kmerUtils(){}

	~kmerUtils(){}

	void set_kmer_length(int len){
		klen = len;
	}

	refKmerList get_kmer_list(){																									
		return klist;
	}

	void genkmers_from_all_transcripts(singleSeqList *tlist){
		//--------------
		outfile.open("data.txt",std::ofstream::binary|std::ofstream::out);
		if (!outfile.is_open())
		{
			std::cerr << "can't open outfile data.txt" << std::endl;
		}
		//--------------

		std::vector<point> transcripts = tlist->seq;
		std::vector<std::string> name = tlist->name;
		int tind;
		auto t1 = std::chrono::high_resolution_clock::now();
		for (tind = 0; tind < tlist->seq.size(); ++tind){
			if (transcripts[tind].size() >= klen){
				genkmer_from_one_transcript(transcripts[tind],name[tind]);
			}	
		}
		auto t2 = std::chrono::high_resolution_clock::now();
		double split_kmer_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
		std::cout << "split_kmer_time:" << split_kmer_time << std::endl;

		std::cout << "ref kmer count:" << klist.seq.size() << std::endl;

		//-----------
		outfile.close();
		//-----------
	}
/*
	// TODO SINGLE_SEQ_LIST preadlist
	void gen_matfile(PAIR_SEQ_LIST preadlist){
		mat_t *matfp;
		matvar_t *matvar;
		size_t refnum = klist.size();
		size_t readnum = preadlist.size();
		size_t refdims[2] = {refnum,klen};
		size_t readdims[2] = {readnum,klen};
		int8_t xtrain[refnum][klen];
		int8_t xtest1[readnum][klen];
		int8_t xtest2[readnum][klen];

		for (int i = 0; i < refnum; ++i){
			for (int j = 0; j < klen; ++j){
				xtrain[i][j] = klist[i].seq[j];
			}
		}

		for(int i = 0;i < readnum;i++){
			for (int j = 0; j < klen; ++j){
				xtest1[i][j] = preadlist[i].first.seq[j];
				xtest2[i][j] = preadlist[i].second.seq[j];
			}
		}
		//create ref matfile
		matfp = Mat_CreateVer("ref_Kmer_Mat.mat",NULL,MAT_FT_DEFAULT);
		if (matfp == NULL){
			std::cerr << "Error creating MAT file: ref_kmers_mat.mat." << std::endl;
		}
		matvar = Mat_VarCreate("xTrain",MAT_C_INT8,MAT_T_INT8,2,refdims,xtrain,0);
		if (matvar == NULL){
			std::cerr << "Error creating variable for xTrain." << std::endl;
		}else{
			Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
			Mat_VarFree(matvar);
		}
		Mat_Close(matfp);
		std::cout << "generate RefMatfile successfull." <<  klist.size() << std::endl;

		//create read matfile
		matfp = Mat_CreateVer("read_Mat.mat",NULL,MAT_FT_DEFAULT);
		if (matfp == NULL){
			std::cerr << "Error creating MAT file: read_Mat.mat." << std::endl;
		}
		matvar = Mat_VarCreate("xTest1",MAT_C_INT8,MAT_T_INT8,2,readdims,xtrain,0);
		if (matvar == NULL){
			std::cerr << "Error creating variable for xTest1." << std::endl;
		}else{
			Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
			Mat_VarFree(matvar);
		}
		matvar = Mat_VarCreate("xTest2",MAT_C_INT8,MAT_T_INT8,2,readdims,xtrain,0);
		if (matvar == NULL){
			std::cerr << "Error creating variable for xTest2." << std::endl;
		}else{
			Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
			Mat_VarFree(matvar);
		}
		std::cout << "generate ReadMatfile successfull." <<  preadlist.size() << std::endl;
	}*/
};

#endif
#pragma once
#include "algorithm"
#include "Mapping.hpp"

/*#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"*/

struct res_analysis
{
    int ref_of_read;
    int pos_of_ref_of_read_1;
};

struct res_analysis_intersection:res_analysis
{
    int pos_of_ref_of_read_2;

    void equal_to_pair_1(res_analysis & rai)
    {
    	ref_of_read = std::move(rai.ref_of_read);
    	pos_of_ref_of_read_1 = std::move(rai.pos_of_ref_of_read_1);
    	pos_of_ref_of_read_2 = 0;
    }
    void equal_to_pair_2(res_analysis & rai)
    {
    	ref_of_read = std::move(rai.ref_of_read);
    	pos_of_ref_of_read_2 = std::move(rai.pos_of_ref_of_read_1);
    	pos_of_ref_of_read_1 = 0;
    }
};

struct SAM_format{
	int oidx;

	std::string qname;
	std::vector<std::string> rname;
	std::string cigar;
	std::string rnext;
	std::string seq;
	std::string qual;
	std::vector<int> flag;
	std::vector<int> pos;
	std::vector<int> mapq;
	std::vector<int> pnext;
	std::vector<int> tlen;
	template<typename OStream>
    friend OStream &operator<<(OStream &os, const SAM_format &s){	
    	os << s.qname << '\t' // QNAME
			<< s.flag[s.oidx] << '\t' // FLAGS
			<< s.rname[s.oidx] << '\t' // RNAME
			<< s.pos[s.oidx] << '\t' // pos (1-based)
			<< s.mapq[s.oidx]  << '\t' // MAPQ
			<< s.cigar << '\t' // CIGAR
			<< s.rnext  << '\t' // MATE NAME
			<< s.pnext[s.oidx] << '\t' // MATE pos
			<< s.tlen[s.oidx] << '\t' // TLEN
			<< s.seq << '\t' // SEQ
			<< "*\t" << '\n';// QSTR
    }

};

class SAMwriter
{
private:
//	std::shared_ptr<spdlog::logger> samlog;

	std::vector<std::vector<int>> ref_of_merged_res;
	std::vector<std::vector<int>> pos_1_of_ref_of_merged_res;
	std::vector<std::vector<int>> pos_2_of_ref_of_merged_res;
	std::vector<int> res_from_which_pair;

    std::vector<bool> is_read_1_rev;
    std::vector<bool> is_read_2_rev;
	int dim;
	int read_size;

	int flag_1[4] = {4,8,8,2};
	int flag_2[2] = {32,16};

	FastxParser fparser;

public:
	SAMwriter(int dim,int read_size);

	void Transfer_Info_From_Mapping(std::vector<bool> &is_read_1_rev,std::vector<bool> &is_read_2_rev);

	static bool greater_comp(res_analysis a,res_analysis b){
		return a.ref_of_read < b.ref_of_read;
	}

	void Analyse_Result(int pair,std::vector<size_t> &code_bucket,std::vector<int> &loc_to_ref, std::vector<size_t> &ref_start,
						std::vector<std::vector<res_analysis>> &ra_of_read,std::vector<int> &mres_min_dis,region_profile *rpro);

    void Analyse_Result_Pair(region_profile &rp);

    void Res_intersection(std::vector<res_analysis> ra_of_read_1,std::vector<res_analysis> ra_of_read_2,std::vector<res_analysis_intersection> &intersection);
   
	void Merge_Result(std::vector<std::vector<res_analysis>> &ra_of_read_1,std::vector<std::vector<res_analysis>> &ra_of_read_2,
					std::vector<int> &mres_1_min_dis,std::vector<int> &mres_2_min_dis);


	void Set_Seq_Of_SAM(bool is_read_rev,SAM_format &read,REAL_TYPE *seq);
	
	void Set_Flag(SAM_format &read,int idx,int is_rc,int rdx);

	void SetUmapped(SAM_format &read, SAM_format &next_read, int size, int i);

	void Generate_SAM(Points &read_1_buff,Points &read_2_buff,std::string sam_file);	
};
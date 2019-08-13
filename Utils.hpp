#ifndef _UTILS_
#define _UTILS_
#include "Parameters.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/access.hpp> 

static constexpr int8_t rc_table[128] = {
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // A C G  to  T G C
        78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // T U to A
        78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // a c g to T G C
        78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78  // t u to A, 78 is N
};
/*
static constexpr int rc_table[128] = {
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 52, 78, 56, 78,  78,  78, 54, 78, 78, 78, 78,  78, 78, 48, 78, // A C G  to  T G C
        78, 78,  78, 78,  50, 50, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // T U to A
        78, 52, 78, 56, 78,  78,  78, 54, 78, 78, 78, 78,  78, 78, 48, 78, // a c g to T G C
        78, 78,  78, 78,  50, 50, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78  // t u to A, 78 is N
};*/

static constexpr REAL_TYPE rc_itoi_table[9] = {
        0, 0, 4, 0, 2, 0, 8, 0, 6
};

static constexpr char rc_itoic_table[9] = {
        0, 0, 52, 0, 50, 0, 56, 0, 54
};

static constexpr char itos_table[9] = {
        78, 78,  65, 78,  84,  78,  67, 78, 71

};

static constexpr char rc_ictos_table[9] = {
        78, 78,  84, 78,  65,  78,  71, 78,  67

};

static constexpr char rc_ictoic_table[80] = {
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  52, 78,  50,  78,  56, 78,  54, 78, 78, 78,  78, 78, 78, 78 

};

//change base to int  ATCG to 1234 character
// A=2=50 T=4=52 C=6=54 G=8=56 
static constexpr char stoic_table[128] = {
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 50, 0, 54, 0,  0,  0, 56, 0, 0, 0, 0,  0, 0, 48,0, 
        0, 0,  0, 0,  52, 52, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 50, 0, 54, 0,  0,  0, 56, 0, 0, 0, 0,  0, 0, 48, 0, 
        0, 0,  0, 0,  52, 52, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0  
};

//change base to REALTYPE.  ATCG to 1234 number
static constexpr REAL_TYPE stoi_table[128] = {
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 2, 0, 6, 0,  0,  0, 8, 0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  4, 4, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 2, 0, 6, 0,  0,  0, 8, 0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  4, 4, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0  
};

//change int character to base  character
static constexpr int8_t ictos_table[80] = {
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, 
        78, 78,  65, 78,  84,  78,  67, 78,  71, 78, 78, 78,  78, 78, 78, 78 

};


//change int character to int
static constexpr int8_t ictoi_table[80] = {
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  0, 0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
        0, 0,  2, 0,  4,  0,  6, 0,  8, 0, 0, 0,  0, 0, 0, 0 
};

static constexpr char itoic_table[80] = {
        48, 0,  50, 0,  52,  0,  54, 0,  56, 57
};

class SArray
{
public:
    int *con_ptr;
    int size;
    std::vector<int> con;
    SArray(){};
    ~SArray(){};

    void tran2vec()
    {
        for (int i = 0; i < size; ++i)
        {
            con.push_back(std::move(con_ptr[i]));
        }
    }

    template<class Archive>
    void serialize(Archive &ar)
    {
        ar(con);
    }
};

struct region_profile
{
    std::vector<uint64_t> rkmer_idx;
    std::vector<size_t> region_start_idx;
    std::vector<size_t> region_end_idx;

    //record the start idx of a bucket in file according to an unique code  after inducing the region size
    std::vector<std::vector<size_t>> code_bucket_idx;

    template<class Archive>
    void serialize(Archive &ar)
    {
        ar(rkmer_idx,region_start_idx,region_end_idx,code_bucket_idx);
    }
};

struct mapped_res
{
    std::vector<std::vector<int>> min_code_idx;
    std::vector<int> min_dis;

    template<class Archive>
    void serialize(Archive &ar)
    {
        ar(min_code_idx,min_dis);
    }
};

struct singleSeqList{	
	std::vector<std::string> seq;
	std::vector<std::string> name;
	std::vector<std::string> qual;
};

struct pairSeqList{	
	std::vector<std::pair<std::string,std::string> > seq;
	std::vector<std::string> name;
	std::vector<std::pair<std::string,std::string> > qual;
};

#endif 
#pragma once

#include <divsufsort.h> 

#include "BooPHF.h"
#include "Utils.hpp"
#include "FastxParser.hpp"
#include "BinaryHash.h"


// custom hash for BBhash
class Custom_string_Hasher
{
public:
    // the class should have operator () with this signature :
    uint64_t operator ()   (std::string key, uint64_t seed=0) const
    {
        

        uint64_t hash  =  hash_fn(key);
        hash ^= seed;
        
        return hash;
    }
    
     std::hash<std::string> hash_fn;
};

class Mapping
{
private:
    // then tell BBhash to use this custom hash : (also appears below, line 104)
    typedef boomphf::mphf< std::string, Custom_string_Hasher > sa_boophf_t;
    // BBhash perfect hash function
    sa_boophf_t * sa_bphf = NULL;

    const int buffer_size = 500000;

    std::vector<std::string> region_kmer;

    std::vector<std::vector<int>> ref_of_read_1;
    std::vector<std::vector<int>> ref_of_read_2;

    std::vector<std::string> rcode_of_read_1;
    std::vector<std::string> rcode_of_read_2;

    std::vector<uint64_t> forward_read_1_region;  
    std::vector<uint64_t> rc_read_1_region;    
    std::vector<uint64_t> forward_read_2_region;  
    std::vector<uint64_t> rc_read_2_region;  

    std::vector<bitset<BCODE_LEN>> read_code_1;
    std::vector<bitset<BCODE_LEN>> read_code_2;
    std::vector<bitset<BCODE_LEN>> rev_read_code_1;
    std::vector<bitset<BCODE_LEN>> rev_read_code_2;

    std::vector<std::string> ref_name;
    std::vector<int> loc_to_ref;

    std::string ref_string;
    SArray sarry;
    SphericalHashing src_sh;
    FastxParser fparser;

public:
    //-------- preserve after mapping------

    std::vector<bool> is_read_1_rev;
    std::vector<bool> is_read_2_rev;
    region_profile *rpro;

    int read_size;
    int dim;
    Mapping(int dim);

    Mapping(int dim,region_profile &rpro);
    ~Mapping();

    void Learn_Spherical_Hashing(Points &buff,int code_len,int seg_len);

    void Load_Spherical_Hashing(int rsize,int code_len,int seg_len);

    void Compute_All_Ref_Code();

    void Merge_Duplicated_In_Region();

    void Load_SA(int seg_len);

    // change the library of constructing sa_ptr, then set idx_name to "tmp/all_ref_seq.bin"
    void Suffix_Array(int seg_len);

    void REAL_TYPE_to_String(std::string &seq,REAL_TYPE *d,bool is_rc);
    
    void Get_Read_Region(Points &read_1_buff,Points &read_2_buff);

    void Resize_Hashcode(std::vector<bitset<BCODE_LEN>> &ref_code_vec);

    std::pair<int,std::vector<int>> Mapping_Process(size_t read_region, std::vector<bitset<BCODE_LEN>> &reduced_region_code,bitset<BCODE_LEN> bCodeRead,int pair);

    void Hashing_Reads(Points &read_1_buff,Points &read_2_buff);

    int Hash_Mapping_with_SA();

    void Load_Ref_Info();

    void Output_Result(int pair,Points &read_ptr);
};
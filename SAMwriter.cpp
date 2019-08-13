#include <iostream>
#include <fstream>

#include "algorithm"
#include "SAMwriter.hpp"
SAMwriter::SAMwriter(int dim,int read_size)
{
	this->dim = dim;
	this->read_size = read_size;
	fparser.initialize(dim);
};

void SAMwriter::Transfer_Info_From_Mapping(std::vector<bool> &is_read_1_rev,std::vector<bool> &is_read_2_rev)
{
	this->is_read_1_rev = std::move(is_read_1_rev);
	this->is_read_2_rev = std::move(is_read_2_rev);
}

void SAMwriter::Analyse_Result(int pair,std::vector<size_t> &code_bucket,std::vector<int> &loc_to_ref, std::vector<size_t> &ref_start,
					std::vector<std::vector<res_analysis>> &ra_of_read,std::vector<int> &mres_min_dis,region_profile *rpro)
{
    std::ifstream tmp_loc;
    std::ifstream tmp_dis;

    std::ifstream read_region_file;
    if (pair == PAIR_1)
    {
        tmp_loc.open(PAIR_1_LOC_FILE);
        tmp_dis.open(PAIR_1_DIS_FILE);
        
        read_region_file.open(PAIR_1_RES_REGION_FILE);
    }else if (pair == PAIR_2)
    {
        tmp_loc.open(PAIR_2_LOC_FILE);
        tmp_dis.open(PAIR_2_DIS_FILE);
        read_region_file.open(PAIR_2_RES_REGION_FILE);
    }

    mapped_res mres;
    // load mapping result from disk using cereal
    {
        cereal::BinaryInputArchive ar_loc(tmp_loc);
        cereal::BinaryInputArchive ar_dis(tmp_dis);
        ar_loc(CEREAL_NVP(mres.min_code_idx));
        ar_dis(CEREAL_NVP(mres.min_dis));
    };

	std::vector<uint64_t> read_region;
    {
        cereal::BinaryInputArchive ar_read(read_region_file);
        ar_read(read_region);
    }

    Stopwatch T0("");
    T0.Reset();     T0.Start();

    size_t ref_loc = 0;
    int loc_size = 0;
    int read_size = mres.min_code_idx.size();
    //std::cout << read_size << std::endl;
    ra_of_read.resize(read_size);
    std::vector<res_analysis> ra_vec;
    res_analysis ra;

#ifdef USE_PARALLELIZATION
    #pragma omp parallel for private(ra_vec,loc_size,ref_loc,ra) num_threads(Parameter::thread)
#endif
    for (unsigned int qIndex = 0;qIndex < read_size;++qIndex)
    {      
        if (mres.min_code_idx[qIndex][0] >= 0)
        {
        	for (int j = 0; j < mres.min_code_idx[qIndex].size(); ++j)
        	{
                if (mres.min_code_idx[qIndex][j] == rpro->code_bucket_idx[read_region[qIndex]].size() - 1)
                {
                	if (read_region[qIndex] + 1 < rpro->code_bucket_idx.size() - 1)
                    {
                    	loc_size = rpro->code_bucket_idx[read_region[qIndex] + 1][0] 
                        	- rpro->code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j]];
                    }else
                    {
                    	loc_size  = 1;
                    }
                }else
                {
                    loc_size = rpro->code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j] + 1] 
                        - rpro->code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j]];
                }
                for (unsigned int i = 0; i < loc_size; ++i)
                {
                    ref_loc = code_bucket[(rpro->code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j]] + i)];
                    if (ref_loc > 0)
                    {
                    	ra.ref_of_read = loc_to_ref[ref_loc];
                    	ra.pos_of_ref_of_read_1 = ref_loc - ref_start[loc_to_ref[ref_loc]] + 1;
                        ra_vec.push_back(ra);
                    }
                }
        	}
            ra_of_read[qIndex] = ra_vec;
            ra_vec.clear();
        }
    }

    T0.Stop();
    mres_min_dis = std::move(mres.min_dis);
}

void SAMwriter::Analyse_Result_Pair(region_profile &rp)
{
	region_profile *rpro = &rp;

	std::vector<size_t> code_bucket;
	std::vector<int> loc_to_ref;
	std::vector<size_t> ref_start;

	std::vector<int> mres_1_min_dis,mres_2_min_dis;
	std::vector<std::vector<res_analysis>> ra_of_read_1,ra_of_read_2;

	Stopwatch T0("");
    T0.Reset();     T0.Start();

    {
        std::ifstream loc_to_ref_file(MERGE_REF_POS_FILE);
        cereal::BinaryInputArchive ar_ref_pos(loc_to_ref_file);
        ar_ref_pos(loc_to_ref);
    }

    {
    	std::ifstream ref_start_file(MERGE_REF_START_FILE);
    	cereal::BinaryInputArchive ar(ref_start_file);
    	ar(ref_start);
    }
    T0.Stop();
    printf("- Load Ref Info Finished (%f seconds)\n",T0.GetTime() );

	T0.Reset();     T0.Start();
    std::ifstream bucket_file("bin/code_bucket.bin");
    bucket_file.seekg(0,bucket_file.end);
    int bucket_size = bucket_file.tellg() / sizeof(size_t);
    code_bucket.resize(bucket_size);
    bucket_file.seekg(0,bucket_file.beg);
    bucket_file.read(reinterpret_cast<char*>(&code_bucket[0]),bucket_size * sizeof(size_t));
//    MemoryMapped code_bucket("bin/code_bucket.bin");

    T0.Stop();
	printf("- Load Code Bucket Info Finished (%f seconds)\n",T0.GetTime() );

	T0.Reset();     T0.Start();
	Analyse_Result(PAIR_1,code_bucket,loc_to_ref,ref_start,ra_of_read_1,mres_1_min_dis,rpro);
	Analyse_Result(PAIR_2,code_bucket,loc_to_ref,ref_start,ra_of_read_2,mres_2_min_dis,rpro);

	bucket_file.close();
//	code_bucket.close();
	T0.Stop();
	printf("- Analyse Results Finished (%f seconds)\n",T0.GetTime() );

	Merge_Result(ra_of_read_1,ra_of_read_2,mres_1_min_dis,mres_2_min_dis);
}

void SAMwriter::Res_intersection(std::vector<res_analysis> ra_of_read_1,std::vector<res_analysis> ra_of_read_2,std::vector<res_analysis_intersection> &intersection)
{
	res_analysis_intersection rai;
	int idx = 0;
	for (int i = 0; i < ra_of_read_1.size(); ++i)
	{
		for (; idx < ra_of_read_2.size(); ++idx)
		{
			if (ra_of_read_1[i].ref_of_read == ra_of_read_2[idx].ref_of_read)
			{
				rai.ref_of_read = ra_of_read_1[i].ref_of_read;
				rai.pos_of_ref_of_read_1 = ra_of_read_1[i].pos_of_ref_of_read_1;
				rai.pos_of_ref_of_read_2 = ra_of_read_2[idx].pos_of_ref_of_read_1;
				intersection.push_back(rai);
			}else if (ra_of_read_1[i].ref_of_read < ra_of_read_2[idx].ref_of_read)
			{
				break;
			}
		}
	}
}
void SAMwriter::Merge_Result(std::vector<std::vector<res_analysis>> &ra_of_read_1,std::vector<std::vector<res_analysis>> &ra_of_read_2,
				std::vector<int> &mres_1_min_dis,std::vector<int> &mres_2_min_dis)
{	
	Stopwatch T0("");
    T0.Reset();     T0.Start();

    int read_size = ra_of_read_1.size();
    ref_of_merged_res.resize(read_size);
    pos_1_of_ref_of_merged_res.resize(read_size);
    pos_2_of_ref_of_merged_res.resize(read_size);
    res_from_which_pair.resize(read_size);

    std::vector<res_analysis_intersection> intersection;
    std::vector<res_analysis_intersection>::iterator it;
    std::vector<int> pos_1;
    std::vector<int> pos_2;
    std::vector<int> ref;

    int half_right = 0;
    int error = 0;

#ifdef USE_PARALLELIZATION
    #pragma omp parallel for private(pos_1,pos_2,intersection,it,ref) reduction(+:half_right,error) num_threads(Parameter::thread)
#endif
    for (int i = 0; i < read_size; ++i)
    {
    	ref.clear();
        pos_1.clear();
        pos_2.clear();
        intersection.clear();

        sort(ra_of_read_1[i].begin(),ra_of_read_1[i].end(),greater_comp);
        sort(ra_of_read_2[i].begin(),ra_of_read_2[i].end(),greater_comp);
         Res_intersection(ra_of_read_1[i],ra_of_read_2[i],intersection);

        if (intersection.size() == 0)
        {
        	intersection.clear();
        	if (mres_1_min_dis[i] <= Parameter::tolerance && mres_2_min_dis[i] >= mres_1_min_dis[i])
        	{
        		intersection.resize(ra_of_read_1[i].size());
        		for (int j = 0; j < ra_of_read_1[i].size(); ++j)
        		{
        			intersection[j].equal_to_pair_1(ra_of_read_1[i][j]);

        		}
        		half_right++;
        		res_from_which_pair[i] = 1;
        	}else if (mres_2_min_dis[i] <= Parameter::tolerance && mres_1_min_dis[i] >= mres_2_min_dis[i])
        	{
        		intersection.resize(ra_of_read_2[i].size());
        		for (int j = 0; j < ra_of_read_2[i].size(); ++j)
        		{
        			intersection[j].equal_to_pair_2(ra_of_read_2[i][j]);
        		}
        		half_right++;
        		res_from_which_pair[i] = 2;
        	}else if (mres_1_min_dis[i] > Parameter::tolerance && mres_2_min_dis[i] > Parameter::tolerance)
        	{
        		error++;
        		res_from_which_pair[i] = 0;
        		pos_1.push_back(0);
        		pos_2.push_back(0);
        		// TODO
        		ref.push_back(-1);
        	}
        }else
        {
        	res_from_which_pair[i] = 3;
        }
        for (it = intersection.begin(); it != intersection.end(); ++it)
        {
            ref.push_back(std::move((*it).ref_of_read));
        	pos_1.push_back(std::move((*it).pos_of_ref_of_read_1));
        	pos_2.push_back(std::move((*it).pos_of_ref_of_read_2));
        }

    	ref_of_merged_res[i] = std::move(ref);
    	pos_1_of_ref_of_merged_res[i] = std::move(pos_1);
    	pos_2_of_ref_of_merged_res[i] = std::move(pos_2);
    }
    T0.Stop();

    printf("- Merge Result Finished (%f seconds)\n",T0.GetTime() );

    std::cout << "- both  hamming distances are larger than "<< Parameter::tolerance << ":" << error << std::endl;
    std::cout << "- half mapping:" << half_right << std::endl;
}


void SAMwriter::Set_Seq_Of_SAM(bool is_read_rev,SAM_format &read,REAL_TYPE *seq)
{
	std::string rev_read;
	std::string seq_str;
	seq_str.resize(dim);
	for (int i = 0; i < dim; ++i)
	{
		seq_str[i] = itos_table[(int8_t)seq[i]];
	}
	if (is_read_rev)
    {
        fparser.reverse_complete(seq_str,rev_read);
    //    std::cout << seq_str << std::endl << rev_read << std::endl;
        seq_str = rev_read;
    }	
    read.seq = seq_str;
}

void SAMwriter::Set_Flag(SAM_format &read,int idx,int is_rc,int rdx)
{	
	int flag = 0;
	flag += flag_1[res_from_which_pair[idx]];

	if (read.tlen[rdx] < 0)
	{
		flag += 64;
	}else if (read.tlen[rdx] > 0)
	{
		flag += 128;
	}

	flag += flag_2[is_rc];
	if (is_rc)
	{
		flag += 16;
	}else
	{
		flag += 32;
	}
	read.flag.push_back(flag);
}

void SAMwriter::SetUmapped(SAM_format &read, SAM_format &next_read, int size, int i)
{
	next_read.qname = read.qname;
	read.rname.push_back("*");
	next_read.rname.push_back("*");
	read.pos.push_back(0);
	next_read.pos.push_back(0);
	read.pnext.push_back(0);
	next_read.pnext.push_back(0);
	read.tlen.push_back(0);
	next_read.tlen.push_back(0);

	Set_Flag(read,i,is_read_1_rev[size + i],0);
	Set_Flag(next_read,i,is_read_2_rev[size + i],0);

	read.rnext = "*";
	next_read.rnext = "*";
	read.mapq.push_back(255);
	next_read.mapq.push_back(255);


	read.cigar = "*";
	next_read.cigar = "*";
	read.tlen.push_back(0);
	next_read.tlen.push_back(0);
}

void SAMwriter::Generate_SAM(Points &read_1_buff,Points &read_2_buff,std::string sam_file)
{
	std::ofstream samfile(sam_file);

	std::vector<std::string> read_name;
    {
        std::ifstream read_name_file(PAIR_1_NAME_FILE);
        cereal::BinaryInputArchive ar(read_name_file);
        ar(read_name);
    }

    std::vector<std::string> ref_name;
	{
        std::ifstream ref_name_file(REF_NAME_FILE);
        cereal::BinaryInputArchive ar_ref_name(ref_name_file);
        ar_ref_name(ref_name);
    }

//    read_1_buff.Initialize_MemoryMapped(USED_READ_FILE_NAME_1,read_size,dim);
//	read_2_buff.Initialize_MemoryMapped(USED_READ_FILE_NAME_2,read_size,dim);

    Stopwatch T2("");
    T2.Reset();     T2.Start();

    std::vector<SAM_format> read,next_read;
    std::vector<int> ref;
    std::string dim_str = std::to_string(dim);
    int ref_size = 0;
    int buffer_size = READ_BUFFER_SIZE;
    REAL_TYPE *read_1,*read_2;
    int pos_1 = 0,pos_2 = 0;
    std::ofstream out("read_test.txt");
    bool test = false;

    for (int size = 0; size < read_size; size += buffer_size)
    {
    	if (read_size - size > READ_BUFFER_SIZE)
    	{
    		buffer_size = READ_BUFFER_SIZE;
    	}else
    	{
    		buffer_size = read_size - size;
    	}
    	read.clear();
	    next_read.clear();
    	read.resize(buffer_size);
    	next_read.resize(buffer_size);
	#ifdef USE_PARALLELIZATION
        #pragma omp parallel for private(ref,ref_size,read_1,read_2,pos_2,pos_1) num_threads(Parameter::thread)
    #endif
	    for (int i = 0; i < buffer_size; ++i)
	    {
	    	read_1 = new REAL_TYPE [dim];
	    	read_2 = new REAL_TYPE [dim];
	    	read_1_buff.A_Read(size + i,read_1);
	    	read_2_buff.A_Read(size + i,read_2);

	    	read[i].qname = std::move(read_name[size + i]);

	    	if (res_from_which_pair[size + i] == 1)
	    	{
	    		read[i].cigar = dim_str + "M";
		    	Set_Seq_Of_SAM(is_read_1_rev[size + i],read[i],read_1);
		    	ref = std::move(ref_of_merged_res[size + i]);
		    	if (ref.size() > 100)
		    	{
		    		ref_size = 100;
		    	}else
		    	{
		    		ref_size = ref.size();
		    	}
	    		read[i].rnext = "=";
		    	next_read[i].rnext = "=";
		    	for (int j = 0; j < ref_size; ++j)
		    	{
		    		read[i].tlen.push_back(0);
		    		read[i].pnext.push_back(0);
		    		pos_1 = pos_1_of_ref_of_merged_res[size + i][j] - Parameter::skip;
		    		if (pos_1 <= 0)
		    		{
		    			if (ref_size > 1)
		    			{
		    				continue;
		    			}else
		    			{
		    				res_from_which_pair[size + i] = 0;
		    				SetUmapped(read[i], next_read[i], size, i);
		    				break;
		    			}
		    		}
	    			read[i].pos.push_back(pos_1);

		    		read[i].mapq.push_back(1);
		    		next_read[i].mapq.push_back(1);

		    		Set_Flag(read[i],i,is_read_1_rev[size + i],j);
		    		if (ref[j] >= 0)
		    		{
		    			read[i].rname.push_back(ref_name[ref[j]]);
		    		}else
		    		{
		    			read[i].rname.push_back("*");
		    		}
		    	//	samfile << read;
		    	}
	    	}else if (res_from_which_pair[size + i] == 2)
	    	{
	    		read[i].cigar = dim_str + "M";
	    		Set_Seq_Of_SAM(is_read_2_rev[size + i],read[i],read_2);
	    		ref = std::move(ref_of_merged_res[size + i]);
	    		if (ref.size() > 100)
		    	{
		    		ref_size = 100;
		    	}else
		    	{
		    		ref_size = ref.size();
		    	}
		    	read[i].rnext = "=";
		    	next_read[i].rnext = "=";
		    	for (int j = 0; j < ref_size; ++j)
		    	{
		    		read[i].tlen.push_back(0);
		    		read[i].pnext.push_back(0);
		    		pos_1 = pos_2_of_ref_of_merged_res[size + i][j] - Parameter::skip;
		    		if (pos_1 <= 0)
		    		{
		    			if (ref_size > 1)
		    			{
		    				continue;
		    			}else
		    			{
		    				res_from_which_pair[size + i] = 0;
		    				SetUmapped(read[i], next_read[i], size, i);
		    				break;
		    			}
		    		}
	    			read[i].pos.push_back(pos_1);

		    		read[i].mapq.push_back(1);
		    		next_read[i].mapq.push_back(1);

	    			Set_Flag(read[i],i,is_read_2_rev[size + i],j);
		    		if (ref[j] >= 0)
		    		{
		    			read[i].rname.push_back(ref_name[ref[j]]);
		    		}else
		    		{
		    			read[i].rname.push_back("*");
		    		}
		    	//	samfile << read;
		    	}
	    	}else if (res_from_which_pair[size + i] == 3)
	    	{
	    		read[i].cigar = dim_str + "M";
	    		next_read[i].qname = read[i].qname;
		    	next_read[i].cigar = dim_str + "M";

	    		Set_Seq_Of_SAM(is_read_1_rev[size + i],read[i],read_1);
	    		Set_Seq_Of_SAM(is_read_2_rev[size + i],next_read[i],read_2);
	    		ref = std::move(ref_of_merged_res[size + i]);
	    		read[i].rnext = "=";
		    	next_read[i].rnext = "=";
		    	for (int j = 0; j < ref.size(); ++j)
		    	{
		    		pos_1 = pos_1_of_ref_of_merged_res[size + i][j] - Parameter::skip;
		    		pos_2 = pos_2_of_ref_of_merged_res[size + i][j] - Parameter::skip;
		    		if (pos_1 <= 0 || pos_2 <= 0)
		    		{
		    			if (ref_size > 1)
		    			{
		    				continue;
		    			}else
		    			{
		    				res_from_which_pair[size + i] = 0;
		    				SetUmapped(read[i], next_read[i], size, i);
		    				break;
		    			}
		    		}
		    		/*if (pos_1 <= 0 )
		    		{
		    			pos_1 = pos_1_of_ref_of_merged_res[size + i][j];
		    		}
		    		if (pos_2 <= 0 )
		    		{
		    			pos_2 = pos_2_of_ref_of_merged_res[size + i][j];
		    		}*/
		    		read[i].pos.push_back(pos_1);
		    		next_read[i].pos.push_back(pos_2);

		    		read[i].pnext.push_back(next_read[i].pos[j]);
		    		next_read[i].pnext.push_back(read[i].pos[j]);

		    		read[i].mapq.push_back(1);
		    		next_read[i].mapq.push_back(1);

		    		if (read[i].pos[j] > next_read[i].pos[j])
		    		{
		    			read[i].tlen.push_back(-(read[i].pos[j] + dim - next_read[i].pos[j]));
		    			next_read[i].tlen.push_back(-read[i].tlen[j]);
		    		}else if (read[i].pos[j] < next_read[i].pos[j])
		    		{
		    			read[i].tlen.push_back(next_read[i].pos[j] + dim - read[i].pos[j]);
		    			next_read[i].tlen.push_back(-read[i].tlen[j]);
		    		}else if (read[i].pos[j] == next_read[i].pos[j])
		    		{
		    			read[i].tlen.push_back(dim);
		    			next_read[i].tlen.push_back(-dim);
		    		}
	
		    		Set_Flag(read[i],i,is_read_1_rev[size + i],j);
		    		Set_Flag(next_read[i],i,is_read_2_rev[size + i],j);

	    			read[i].rname.push_back(ref_name[ref[j]]);
	    			next_read[i].rname.push_back(ref_name[ref[j]]);
		    		
		    	//	samfile << read;
		    	//	samfile << next_read;
		    	}
	    	}else if (res_from_which_pair[size + i] == 0)
	    	{
	    		SetUmapped(read[i], next_read[i], size, i);
	    	}
	    	delete read_1;
	    	delete read_2;
	    }

	    for (int i = 0; i < buffer_size; ++i)
	    {
	    //	std::cout <<size+i << std::endl;
	    	for (int j = 0; j < read[i].rname.size(); ++j)
	    	{
	    		if (!test)
	    		{
	    			out << read[i].qname << std::endl;
	    			test = true;
	    		}
	    		read[i].oidx = j;
	    		samfile << read[i];
		    	if (res_from_which_pair[size + i] == 3 || res_from_which_pair[size + i] == 0)
		    	{
		    		next_read[i].oidx = j;
		    		//std::cout << next_read[i].tlen[j] << std::endl;
		    		samfile << next_read[i];
		    	}
	    	}
	    	test = false;
	    }
	}
    samfile.close();
 
	read_2_buff.ReleaseMem();
	read_1_buff.ReleaseMem();
	system("rm tmp/*");
	T2.Stop();
	printf("- Generate SAM File Finished (%f seconds)\n",T2.GetTime());
}
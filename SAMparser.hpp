#include <iostream>
#include <fstream>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include "fastxParser.hpp"
#include "Utils.hpp"

/*#define REAL_TYPE int
#define REF_LOCATION std::std::vector<int>

*/

class SAMparser
{
public:
	singleSeqList *reflist ;
	int read_index;
	int dim;
	int nP;
	std::string seq_1;
	std::string seq_2;
	std::vector<int> loc_1;
	std::vector<int> loc_2;
	std::vector<int> ref_1;
	std::vector<int> ref_2;
	bool again;


	SAMparser(){
		reflist = new singleSeqList();
		read_index = 0;
		dim = 0;
		nP = 0;
		seq_1 = "";
		seq_2 = "";
		again = false;
	}

	void initialize(int _nP){
		nP = _nP;
		dim = DIM;
	/*	rinfo.nP = nP;
		rinfo.dim_ = dim;
		rinfo.read_seq = new REAL_TYPE * [ nP ];
		for(int i=0;i<nP;i++)
		{
			rinfo.read_seq[i] = new REAL_TYPE [ dim ];
		}
		rinfo.read_name = new std::string [nP];
		rmref = new read_mapped_ref[nP];*/
	}


//-------------parse sam files-----------------------------
	int get_refseq(std::string ref_name){
		for (int i = 0; i < reflist->name.size(); ++i){
			if (reflist->name[i] == ref_name){
				return i;
			}
		}
		return -1;
	}

	std::string get_One_Of_Sam(std::ofstream &out_res,std::ifstream &sam,std::string last_read_name){
		std::string line;
		char tmp;
		if (!sam.is_open()){
			std::cerr << "can't open sam file." <<std::endl;
		}
		//skip head
		sam.get(tmp);
		while(tmp == '@'){
			std::getline(sam,line);
			sam.get(tmp);
		}
		sam.seekg(-1,std::ios_base::cur);

		std::getline(sam,line);
		std::stringstream ss(line);
		std::string read_name;
		std::string ref_name;
		int ref_loc;
		int sam_index = 0;
		std::string read_seq;

		ss >> read_name;
		sam_index++;
	//	std::cout << read_name << std::endl;
		for (int i = 0; i < 2; ++i){
			//ref name that one read_idx mapped to
			ss >> ref_name;
			sam_index++;
		}
		//the loc of ref that one read_idx mapped to
		ss >> ref_loc;
		//the index in sam start at 1
		ref_loc--;
		sam_index++;
		while(ss.peek() != EOF){
			ss >> read_seq;
			sam_index++;
			if (sam_index == 10){
				break;
			}
		}
		if (again && read_seq != "*")
		{
			seq_2 = read_seq;
			again= false;
		}
		if (last_read_name == "")
		{
			last_read_name = read_name;
			if (read_seq != "*")
			{
				seq_1 = read_seq;
				again = true;
			}
		}
		if (read_name != last_read_name){
			out_res << ">" << last_read_name << std::endl;
			out_res << "=" << seq_1 << std::endl;
			for (int i = 0; i < ref_1.size(); ++i)
			{
				out_res << "+";
				for (int j = 0; j < dim; ++j){
					out_res << reflist->seq[ref_1[i]][loc_1[i] + j];
				}
				out_res << std::endl;
			}
			ref_1.clear();
			loc_1.clear();

			out_res << ">" << last_read_name << std::endl;
			out_res << "=" << seq_2 << std::endl;
			for (int i = 0; i < ref_2.size(); ++i)
			{
				out_res << "+";
				for (int j = 0; j < dim; ++j){
					out_res << reflist->seq[ref_2[i]][loc_2[i] + j];
				}
				out_res << std::endl;
			}
			ref_2.clear();
			loc_2.clear();

			if (read_seq != "*")
			{
				seq_1 = read_seq;
				again = true;
			}

		/*	if (sam_index < 12){
				if (sam_index == 10 && read_seq != "*"){
					out_res << "=" << read_seq << std::endl;		
				}else{
					std::getline(sam,line);
					ss << line;
					while(ss.peek() != EOF){
						ss >> read_seq;
						sam_index++;
						if (sam_index == 10 && read_seq != "*"){
							out_res << "=" << read_seq << std::endl;
							break;
						}
					}
				}
			}*/
		}	

		if (ref_name != "*")
		{
			int rind = get_refseq(ref_name);
			if (read_seq == seq_1)
			{
				ref_1.push_back(rind);
				loc_1.push_back(ref_loc - SKIP);
			}else
			{
				ref_2.push_back(rind);
				loc_2.push_back(ref_loc - SKIP);
			}

			//print the ref seq that this read_idx mapped to
			//out_res << "+";
			
		/*	for (int i = 0; i < dim; ++i){
				out_res << reflist->seq[rind][ref_loc + i];
			}
			out_res << std::endl;*/
		}
		return read_name;

	/*	//add info to read_mapped_ref
		rmref[read_index].read_name = read_name;
		rmref[read_index].ref_name.push_back(ref_name);
		rmref[read_index].locs_of_ref.push_back(ref_loc);


		//add new read_idx to structure
		rinfo.read_name[read_index] = read_name;
		if (sam_index < 12){
			if (sam_index == 10 && read_seq != '*'){
				//TODO
				for (int i = 0; i < rinfo.dim; ++i){
					rinfo.read_seq[read_index][i] = (REAL_TYPE)( read_seq[k] - 48)
				}
				read_index++;

			}else{
				std::getline(sam,line);
				ss << line;
				while(ss.peek() != EOF){
					ss >> read_seq;
					sam_index++;
					if (sam_index == 10) && read_seq != '*'{
						for (int i = 0; i < rinfo.dim; ++i){
							rinfo.read_seq[read_index][i] = (REAL_TYPE)( read_seq[i] - 48)
						}
						read_index++;
						break;
					}
				}
			}
			
		}*/
		
	}

	void get_Sam(char *ref_file,char *sam_file,char *sam_res_file,int _nP,int _dim){
		//extract ref seq to RAM
		read_fastx(ref_file,reflist);

		std::ifstream sam;
		std::ofstream out_res;
		sam.open(sam_file);
		out_res.open(sam_res_file);
		if (!sam.is_open())
		{
			std::cerr << "can't open sam file" << std::endl;
		}

		std::string last_read_name = "";
		read_index = 0;
		nP = _nP;
		dim = _dim;
		int count = 1;
		while(count <= 100){
			last_read_name = get_One_Of_Sam(out_res,sam,last_read_name);
			count++;
		}
		sam.close();
		out_res.close();
		std::cout << "Parsing Sam Finished ." << std::endl;
	}

};
#include "BinaryHash.h"
#include "Mapping.hpp"
#include "FastxParser.hpp"
#include "SAMwriter.hpp"

#ifdef USE_PARALLELIZATION
#include <omp.h>
#include <iostream>
#include <sstream>
#include<unistd.h>
#endif


int main(int argc, char const *argv[])
//int main()
{	
	std::string par;
	/*Parameter::dim = 48;
	Parameter::bcode_len = 30;
	Parameter::tolerance = 4;
	Parameter::type = "mapping";
	Parameter::read_file_1 = "../SpheriHash/dataset/srrdata/SRR5337025_paired_1.fastq";
	Parameter::read_file_2 = "../SpheriHash/dataset/srrdata/SRR5337025_paired_2.fastq";
	Parameter::sam = "res/sim_100_res.sam";
	Parameter::transcripts_file_name = "/home/yxt/Documents/work/RNA-seq/reference/transcripts/Homo_sapiens.GRCh38.cdna.all.fa";*/
	for (int i = 0; i < argc; ++i)
	{
		par = argv[i];
		//std::cout << par << std::endl;
		if (par.compare("index") == 0 || par.compare("mapping") == 0)
		{
			Parameter::type = argv[i];
		}else if (par.compare("-k") == 0)
		{
			Parameter::kmer = std::atoi(argv[i + 1]);
		}else if (par.compare("-r") == 0)
		{
			Parameter::transcripts_file_name = argv[i + 1];
		}else if (par.compare("-s") == 0)
		{
			Parameter::skip = std::atoi(argv[i + 1]);
		}else if (par.compare("-m") == 0)
		{
			Parameter::region_searching = std::atoi(argv[i + 1]);
		}else if (par.compare("-g") == 0)
		{
			Parameter::ignore = std::atoi(argv[i + 1]);
		}else if (par.compare("-p") == 0)
		{
			Parameter::thread = std::atoi(argv[i + 1]);
		}else if (par.compare("-I") == 0)
		{
			Parameter::read_file_1 = argv[i + 1];
		}else if (par.compare("-i") == 0)
		{
			Parameter::read_file_2 = argv[i + 1];
		}else if (par.compare("-o") == 0)
		{
			Parameter::sam = argv[i + 1];
		}else if (par.compare("-d") == 0)
		{
			Parameter::dim = std::atoi(argv[i + 1]);
		}else if (par.compare("-c") == 0)
		{
			Parameter::bcode_len = std::atoi(argv[i + 1]);
			Parameter::dim = (Parameter::bcode_len + Parameter::region_searching + Parameter::skip + Parameter::ignore) + Parameter::kmer - 1;
		}else if (par.compare("-t") == 0)
		{
			Parameter::tolerance = std::atoi(argv[i + 1]);
		}
	}

	//Parameter::bcode_len = (Parameter::dim - Parameter::region_searching - Parameter::skip - Parameter::ignore) - Parameter::kmer + 1; /// 
	std::cout << "- hash code length:"<< Parameter::bcode_len << std::endl;
	std::cout << "- read dim:"<< Parameter::dim << std::endl;

	Stopwatch T0("");
	T0.Reset();     T0.Start();

	region_profile rpro;

	if (Parameter::type.compare("index") == 0)
	{
		Mapping *map = new Mapping(Parameter::dim,rpro);
		std::cout << "- dim:" << Parameter::dim << std::endl;
		std::cout << "- segment length:" << Parameter::kmer << std::endl;

		map->Suffix_Array(Parameter::kmer);

		T0.Stop();
		printf("- Index Finished (%f seconds)\n\n\n\n",T0.GetTime() );
	}else if (Parameter::type.compare("mapping") == 0)
	{		
    	Points read_1_buff;
    	Points read_2_buff;
    	FastxParser fparser(Parameter::dim);

    	int read_size = fparser.store_reads(read_1_buff, read_2_buff, Parameter::read_file_1, Parameter::read_file_2);

    	Mapping *map = new Mapping(Parameter::dim, rpro);

		map->Load_Spherical_Hashing(read_size,Parameter::bcode_len,Parameter::kmer);


		map->Hashing_Reads(read_1_buff,read_2_buff);
		map->Get_Read_Region(read_1_buff,read_2_buff);

		read_2_buff.ReleaseMem();
		read_1_buff.ReleaseMem();
		
		map->Load_SA(Parameter::kmer);

		map->Hash_Mapping_with_SA();

		/*
		map->Load_Ref_Info();
		map->Output_Result(PAIR_1,read_1_buff);
		map->Output_Result(PAIR_2,read_2_buff);*/

		
	 	SAMwriter *sp = new SAMwriter(Parameter::dim, read_size); 
	 	sp->Transfer_Info_From_Mapping(map->is_read_1_rev,map->is_read_2_rev);//,map->read_1_buff,map->read_2_buff);
	 	delete map;

		sp->Analyse_Result_Pair(rpro);

		read_size = fparser.store_reads(read_1_buff, read_2_buff, Parameter::read_file_1, Parameter::read_file_2);
		
		sp->Generate_SAM(read_1_buff,read_2_buff,Parameter::sam);
		T0.Stop();
		printf("- Total Running Time (%f seconds)\n",T0.GetTime() );
	}
		
	return 0;
}
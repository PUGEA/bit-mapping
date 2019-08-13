#pragma once
#include "string.h"
#include <iostream>

namespace Parameter
{
	extern std::string type;
	extern int bcode_len;
	extern std::string transcripts_file_name;
	extern int tolerance;
	extern int ignore;
	extern int thread;

	extern int skip;
	extern int region_searching;
	extern int kmer;
	extern std::string read_file_1;
	extern std::string read_file_2;
	extern std::string sam;
	extern int dim;
	extern int filter;
};
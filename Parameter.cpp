#include "Parameter.hpp"

namespace Parameter
{
	std::string type = "index";
	int bcode_len = 37;
	std::string transcripts_file_name = "../reference/transcripts/Homo_sapiens.GRCh37.cdna.all.fa";
	int tolerance = 2;
	int ignore = 0;
	int thread = 4;

	int skip = 6;
	int region_searching = 10;
	int kmer = 3;
	std::string read_file_1;
	std::string read_file_2;
	std::string sam;
	int dim = 50;
	int filter = 2;
};
#include "bisearch.h"

#include <boost/program_options.hpp>
#include <fstream>

#include "nucleotide.h"

namespace po = boost::program_options;
using namespace std;

enum{
	SUCCESS = 0,
	ERROR = 1,
};

int main(int argc, char* argv[]){

	string inputfile;
	int product_len_min=100;
	int product_len_max=600;
	int primer_len_min=20;
	int primer_len_max=35;
	PCRCondition cond=PCRCondition();
	float max_tm_diff=8.0f;
	float max_met_tm_diff = 2.5f;
	int score_threshold=100;
	int max_results=200;

	try{

		po::options_description desc("Options");
		desc.add_options()
		("input", po::value<string>(&inputfile)->required(),"input file name")
		("help,h", "help message")
		("product_len_min", po::value<int>(&product_len_min), "minimum pcr product size")
		("product_len_max", po::value<int>(&product_len_max), "maximum pcr product size")
		("primer_len_min", po::value<int>(&primer_len_min), "minimum primer length")
		("primer_len_max", po::value<int>(&primer_len_max), "maximum primer length")
		("max_tm_diff", po::value<float>(&max_tm_diff), "maximum Tm differences between fw and rv primer")
		("max_met_tm_diff", po::value<float>(&max_met_tm_diff), "maximum Tm differences between methyl and unmethyl version.")
		("score_threshold", po::value<int>(&score_threshold), "threshold value for primer pair score.")
		("max_results", po::value<int>(&max_results), "maximum number of results");

		po::positional_options_description positional;
		positional.add("input", 1);


		po::variables_map vm;
		po::store(po::command_line_parser(argc,argv).options(desc).positional(positional).run(), vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			return SUCCESS;
		}

		try{
			po::notify(vm);
		}
		catch(po::required_option& e){
			cout << e.what();
			cout << desc << "\n";
			return ERROR;
		}
		catch(po::error& e){
			cout << e.what();
			cout << desc << "\n";
			return ERROR;
		}
	}
	catch(std::exception& e) 
	{ 
    	std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << std::endl; 
    	return ERROR; 
	}

	std::ifstream f(inputfile);
	if(!f.is_open()){
		std::cerr << "No such file: " << inputfile << std::endl;
	}
	string seq;
	string line;
	while(f.good()){
		getline(f, line);
		if(line[0]=='>')
			continue;
		seq += line;
	}

	try{
		bisearch(seq.c_str(), cout, 
			product_len_min, product_len_max, primer_len_min, primer_len_max,
			cond, max_tm_diff, max_met_tm_diff, score_threshold, max_results);
	}
	catch(std::exception& e) 
	{ 
    	std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << std::endl; 
    	return ERROR; 
	}

	return SUCCESS;
};


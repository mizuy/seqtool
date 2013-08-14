#include "melt_temp.h"
#include "nucleotide.h"

#include <boost/program_options.hpp>
#include <fstream>

const PCRCondition defaultPCRCondition(0.0, 0.05f, 0.0015, 0.000001);
const PCRCondition myPCRCondition(0.033, 0.0, 0.0065, 0.0000005);

namespace po = boost::program_options;
using namespace std;

enum{
	SUCCESS = 0,
	ERROR = 1,
};

int main(int argc, char* argv[]){
    string sequence;

	try{
		po::options_description desc("Options");
		desc.add_options()
          ("sequence", po::value<string>(&sequence)->required(),"primer sequence");

		po::positional_options_description positional;
		positional.add("sequence", 1);

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
    	std::cerr << "During commandline options analysis, unhandled Exception reached the top of main: " 
              << e.what() << std::endl; 
    	return ERROR; 
	}


	try{
      melt_temp::MeltTemp mt(myPCRCondition);
      std::vector<char> seq;
      for(char i : sequence){
        seq.push_back(nucleotide :: base_int(i));
      }
      
      float tm = mt.seq_tm(seq);
      std::cout << "Sequence: " << sequence <<  std::endl;
      std::cout << "Tm = " <<  tm << " degree" <<  std::endl;
	}
	catch(std::exception& e) 
	{ 
    	std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << std::endl; 
    	return ERROR; 
	}

	return SUCCESS;
};


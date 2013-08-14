
#include "bisearch.h"

#include <boost/program_options.hpp>
#include <fstream>

#include "nucleotide.h"

// Na 33mM, K 0mM, Mg 6.5mM, Primer 0.5uM
const PCRCondition defaultPCRCondition(0.0, 0.05f, 0.0015, 0.000001);
const PCRCondition myPCRCondition(0.033, 0.0, 0.0065, 0.0000005);


namespace po = boost::program_options;
using namespace std;

enum{
	SUCCESS = 0,
	ERROR = 1,
};

//http://stackoverflow.com/questions/6240950/platform-independent-dev-null-in-c
template<typename Ch, typename Traits = std::char_traits<Ch> >
struct basic_nullbuf : std::basic_streambuf<Ch, Traits> {
     typedef std::basic_streambuf<Ch, Traits> base_type;
     typedef typename base_type::int_type int_type;
     typedef typename base_type::traits_type traits_type;

     virtual int_type overflow(int_type c) {
         return traits_type::not_eof(c);
     }
};

// convenient typedefs
typedef basic_nullbuf<char> nullbuf;
typedef basic_nullbuf<wchar_t> wnullbuf;

// buffers and streams
// in some .h
//extern std::ostream cnull;
//extern std::wostream wcnull;

// in a concrete .cpp
nullbuf null_obj;
wnullbuf wnull_obj;
std::ostream cnull(&null_obj);
std::wostream wcnull(&wnull_obj);

int main(int argc, char* argv[]){
    string inputfile;
    string outputfile;
    string logfile;
	int product_len_min=100;
	int product_len_max=250;
	int primer_len_min=20;
	int primer_len_max=35;
    float max_tm_diff=4.0f;
    //float max_tm_diff=8.0f;
	float max_met_tm_diff = 2.5f;
    int max_cpg_in_primer=1;
	int score_threshold=100;
	int max_results=200;

	try{
		po::options_description desc("Options");
		desc.add_options()
        ("input", po::value<string>(&inputfile)->required(),"input file name (fasta)")
        ("output", po::value<string>(&outputfile),"output file name (.seqv), default is standard output.")
        ("logfile", po::value<string>(&logfile),"logging file name.")
		("help,h", "show help message")
		("product_len_min", po::value<int>(&product_len_min), "minimum pcr product size")
		("product_len_max", po::value<int>(&product_len_max), "maximum pcr product size")
		("primer_len_min", po::value<int>(&primer_len_min), "minimum primer length")
		("primer_len_max", po::value<int>(&primer_len_max), "maximum primer length")
		("max_tm_diff", po::value<float>(&max_tm_diff), "maximum Tm differences between fw and rv primer")
        ("max_met_tm_diff", po::value<float>(&max_met_tm_diff), "maximum Tm differences between methyl and unmethyl version.")
        ("max_cpg_in_primer", po::value<int>(&max_cpg_in_primer), "maximum number of CpG allowed within the primer.")
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
    	std::cerr << "During commandline options analysis, unhandled Exception reached the top of main: " 
              << e.what() << std::endl; 
    	return ERROR; 
	}

	std::ifstream f(inputfile);
	if(!f.is_open()){
		std::cerr << "No such file: " << inputfile << std::endl;
	}
	string seq;
	string line;
    // TODO: fix fasta parsing.
	while(f.good()){
		getline(f, line);
		if(line[0]=='>' || line[0]=='#')
			continue;
		seq += line;
	}



    std::streambuf * obuf;
    std::ofstream of;
    if(!outputfile.empty()) {
        of.open(outputfile);
        obuf = of.rdbuf();
    } else {
        obuf = std::cout.rdbuf();
    }
    std::ostream output(obuf);

    std::streambuf * lbuf;
    std::ofstream lf;
    if(!logfile.empty()) {
        lf.open(logfile);
        lbuf = lf.rdbuf();
    } else {
        lbuf = cnull.rdbuf();
    }
    std::ostream logging(lbuf);

	try{
		bisearch(seq.c_str(), output, logging, myPCRCondition, 
			product_len_min, product_len_max, primer_len_min, primer_len_max,
			max_tm_diff, max_met_tm_diff, max_cpg_in_primer, score_threshold, max_results);
	}
	catch(std::exception& e) 
	{ 
    	std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << std::endl; 
    	return ERROR; 
	}

	return SUCCESS;
};


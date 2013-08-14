#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

/*
sense
CG -> YG
anti-sense
CG -> CR


met-sense
CG -> CG
C* -> T*
# this==C and next==G -> C(this)
  this==C and next!=G -> T
  this!=C -> this

unmet-sense
CG -> TG
C* -> T*
# just convert every C to T

met-asense
CG -> CG
*G -> *A
# this=G and last==C -> G(this)
  this=G and last!=C -> A
  this!=G -> this

unmet-sense
CG -> CA
*G -> *A
# just convert every G to A

*/


template<char from, char to>
void conv(std::istream& in, std::ostream& out){
  std::string line;
  while(in && getline(in, line)){
    for(int i=0; i<line.size(); i++){
      char b = toupper(line[i]);
      out << ((b==from) ? to : b);
    }
    out << std::endl;
  }
}

void conv_met_pos(std::istream& in, std::ostream& out){
  std::string line;
  bool first=true;
  char b=0;
  while(in && getline(in, line)){
    for(int i=0; i<line.size(); i++){
      char next = toupper(line[i]);
      if(first){
        first = false;
        b = next;
        continue;
      }
      out << ((b=='C' && next!='G') ? 'T' : b);
      b = next;
    }
    out << std::endl;
  }
  if(b)
    out << b << std::endl;
}

void conv_met_neg(std::istream& in, std::ostream& out){
  std::string line;
  char last='N';
  while(in && getline(in, line)){
    for(int i=0; i<line.size(); i++){
      char b = toupper(line[i]);
      out << ((b=='G' && last!='C') ? 'A' : b);
      last = b;
    }
    out << std::endl;
  }
}


void print_usage(){
  std::cout << "convert [pn] [mu] filename" << std::endl;
}

int main(int argc, char **argv){
  if(argc!=4){
    print_usage();
    return -1;
  }
  bool sense = argv[1][0]=='p';
  bool methyl = argv[2][0]=='m';

  std::string name_addition = "-bs_";
  if(sense)
    name_addition += "pos_";
  else
    name_addition += "neg_";
  if(methyl)
    name_addition += "met";
  else
    name_addition += "unm";

  std::string filename(argv[3]);

  std::string record;
  std::ifstream ifs(filename);

  if(!ifs.is_open()){
    std::cerr << "cannot open file: " << filename << std::endl;
    return -1;
  }
  
  std::string line;
  getline(ifs,line);
  if(line[0]=='>'){
    record = line+name_addition;
  }
  else{
    record = "no-record";
    ifs.seekg(0, std::ifstream::beg);
  }
  std::cout << ">" << record << std::endl;

  if(!methyl){
    if(sense)
      conv<'C','T'>(ifs, std::cout);
    else
      conv<'G','A'>(ifs, std::cout);
  }
  else{
    if(sense)
      conv_met_pos(ifs, std::cout);
    else
      conv_met_neg(ifs, std::cout);
  }
  return 0;
}


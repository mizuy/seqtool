#ifndef NUCLEOTIDE_H
#define NUCLEOTIDE_H

#include <string>

namespace nucleotide{
  inline int complement_chr(char x){
    return ((x=='A') ? 'T' : ((x=='T') ? 'A': ((x=='G') ? 'C' : ((x=='C') ? 'G' : 'N'))));
  }

  inline std::string reverse_complement(const std::string& input){
    const int l = input.length();
    std::string ret;
    ret.resize(l);
    for(int i=0; i<l; i++){
      char v = input[l-1-i];
      ret[i] = complement_chr(v);
    }
    return ret;
  }


  template<bool methyl>
  inline void bisulfite(std::string& ret, const std::string& input){
    const int l = input.length();
    ret.resize(l);
    for(int i=0; i<l-1; i++){
        if(input[i]=='C' && !(methyl && input[i+1]=='G'))
          ret[i] = 'T';
        else
          ret[i] = input[i];
    }
    ret[l-1] = (input[l-1]=='C') ? 'T' : input[l-1];
  }

  const int N = -1;
  const int A = 0;
  const int T = 1;
  const int G = 2;
  const int C = 3;

  inline int base_int(char x){
      char c = toupper(x);
      return ((c=='A') ?  A  : ((c=='T') ?  T : ((c=='G') ?  G  : ((c=='C') ? C  :  N ))));
  }

  inline int complement(int x){
    return ((x==A) ?  T  : ((x==T) ?  A : ((x==G) ?  C  : ((x==C) ? G  :  N ))));
  }
};

#endif
#ifndef NUCLEOTIDE_H
#define NUCLEOTIDE_H

#include <string>

#define COMPLEMENT(x) (((x)=='A') ? 'T' : (((x)=='T') ? 'A': (((x)=='G') ? 'C' : (((x)=='C') ? 'G' : 'N'))))
inline std::string reverse_complement(const std::string& input){
  const int l = input.length();
  std::string ret;
  ret.resize(l);
  for(int i=0; i<l; i++){
    char v = input[l-1-i];
    ret[i] = COMPLEMENT(v);
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

#endif
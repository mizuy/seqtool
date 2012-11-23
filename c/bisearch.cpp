#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <cctype>
#include <cassert>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>

#include "bisearch.h"
#include "nucleotide.h"

const bool debug = false;

using namespace std;

const int A = 0;
const int T = 1;
const int G = 2;
const int C = 3;

inline int complement(int x){
  return ((x==A) ?  T  : ((x==T) ?  A : ((x==G) ?  C  : ((x==C) ? G  :  -1 ))));
}

const int anneals[4][4] = {
  {0, 2, 0, 0}, // AT
  {2, 0, 0, 0}, // TA
  {0, 0, 0, 4}, // GC
  {0, 0, 4, 0}, // CG
};
const int anneals_c[4][4] = {
  {2, 0, 0, 0},
  {0, 2, 0, 0},
  {0, 0, 4, 0},
  {0, 0, 0, 4},
};

struct Condition{
public:
  Condition(float a, float b, float c, float d):
    weight(a), optimal(b), minimum(c), maximum(d){}
  float weight;
  float optimal;
  float minimum;
  float maximum;
  inline bool within(float value){return (minimum <= value) && (value <= maximum);}
  inline float score(float value){return weight * abs(value-optimal);}
};

class PrimerCondition{
public:
  PrimerCondition(int min_len=20, int max_len=35, bool bisulfite=true):
    length(0.5, 23., min_len, max_len),
    gc(1.0, bisulfite ? 30 : 50, bisulfite ? 0 : 30, bisulfite ? 60 : 70),
    tm(1.0, 60., 45, 70),
    sa(0.1, 0., 0, 20),
    sea(0.2, 0., 0, 10),
    pa(0.1, 0., 0, 20),
    pea(0.2, 0., 0, 10){
  }
  Condition length;
  Condition gc;
  Condition tm;
  Condition sa;
  Condition sea;
  Condition pa;
  Condition pea;
};

  // ATGC
const float nnt_dh[4][4] = {
  {9.1, 8.6, 7.8, 6.5},   //AA, AT, AG, AC
  {6.0, 9.1, 5.8, 5.6},   //TA, TT, TG, TC
  {5.6, 6.5, 11.0, 11.1}, //GA, GT, GG, GC
  {5.8, 7.8, 11.9, 11.0}, //CA, CT, CG, CC
};
const float nnt_dg[4][4] = {
  {1.55, 1.25, 1.45, 1.40},
  {0.85, 1.55, 1.15, 1.15},
  {1.15, 1.40, 2.30, 2.70},
  {1.15, 1.45, 3.05, 2.30},
};

class MeltTemp{
public:
  MeltTemp(const PCRCondition& cond) : cond(cond) {
    c_salt = cond.na_conc + cond.k_conc + 4*pow(cond.mg_conc,0.5f);
    cc_primer = c_r*c_t0*log(cond.primer_conc); // RT0ln(c)
    cc_salt = 16.6 * log10(c_salt / (1+0.7*c_salt));

    // nnt_dh, nnt_dg validity. NNT of a sequence complemtely equals to that of reverse complement.
    if(debug){
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          int ii = complement(i);
          int jj = complement(j);
          //cout << i << " " << j << " " << ii << " " << jj << endl;
          assert(nnt_dh[i][j] == nnt_dh[jj][ii]);
          assert(nnt_dg[i][j] == nnt_dg[jj][ii]);
        }
      }
    }
  }
  float calc_tm(float sdh, float sdg){
    float dhp = -1000.0 * (2*c_dhe + sdh);
    float dgp = -1000.0 * (2*c_dge + c_dgi + sdg);
    return c_t0 * dhp / (dhp-dgp+cc_primer) + cc_salt - 269.3;
  }
  float seq_tm(const vector<char>& seq){
    float p=0; //enthalpy
    float q=0; //free energy
    for(int i=0; i<seq.size()-1; i++){
      p += nnt_dh[(int)seq[i]][(int)seq[i+1]];
      q += nnt_dg[(int)seq[i]][(int)seq[i+1]];
    }
    return calc_tm(p,q);
  }
private:
  float c_salt;
  float cc_primer;
  float cc_salt;
  const float c_r = 1.987f;
  const float c_t0 = 298.2f;
  const float c_dhe = 5.0; // delta He = 5 kcal/mol
  const float c_dge = 1.0; // delta Ge = 1 kcal/mol
  const float c_dgi = -2.2; // delta Gi = -2.2 kcal/mol
  PCRCondition cond;
};


class Seqint{
public:
  Seqint(const string& seq)
  {
    const int length = seq.length();
    seqint_.resize(length);
    for(int i=0; i<length; i++){
      char c = toupper(seq[i]);
      if(c=='A')
        seqint_[i] = A;
      else if(c=='T')
        seqint_[i] = T;
      else if(c=='G')
        seqint_[i] = G;
      else
        seqint_[i] = C;
    }
    tm_dh_m.resize(length-1);
    tm_dg_m.resize(length-1);
    tm_dh_u.resize(length-1);
    tm_dg_u.resize(length-1);
    // p     q
    // A T G C A T G C
    // |/|/|/|/|/|/|/
    // x x x x x x x
    // p   q
    for(int i=0; i<length-1; i++){
      int p = seqint_[i];
      int q = seqint_[i+1];
      tm_dh_m[i] = nnt_dh[p][q];
      tm_dg_m[i] = nnt_dg[p][q];
      if(p==C) p=T;
      if(q==C) q=T;
      tm_dh_u[i] = nnt_dh[p][q];
      tm_dg_u[i] = nnt_dg[p][q];
    }
    // Tm of seq[p:q] is calc_tm(sum(tm_dh_m[p:q-1]), sum(tm_dg_m[p:q-1]))
  }
  inline bool is_cpg(int i){
    return (i+1<seqint_.size()) ? (seqint_[i]==C and seqint_[i+1]==G) : false;
  }
  inline bool is_gc(int i){
    return seqint_[i]==C or seqint_[i]==G;
  }

  inline int s_c(int p, int q){
    return anneals[(int)seqint_[p]][(int)seqint_[q]];
  }
  inline int s_c_c(int p, int q){
    return anneals_c[(int)seqint_[p]][(int)seqint_[q]];
  }
  /*
  inline float tm_m(int i, int j, MeltTemp& mt){
    float e=0; //enthalpy
    float f=0; //free energy
    for(int i=0; i<j; i++){
      e += tm_dh_m[i];
      f += tm_dg_m[i];
    }
    return mt.calc_tm(e,f);
  }
  inline float tm_u(int i, int j, MeltTemp& mt){
    float e=0; //enthalpy
    float f=0; //free energy
    for(int i=0; i<j; i++){
      e += tm_dh_u[i];
      f += tm_dg_u[i];
    }
    return mt.calc_tm(e,f);
  }
*/
  vector<float> tm_dh_m;
  vector<float> tm_dg_m;
  vector<float> tm_dh_u;
  vector<float> tm_dg_u;
private:
  vector<char> seqint_;
};


class Value{
public:
  Value():len(-1), tm_m(-1), tm_u(-1), gc(-1), sa(-1), sa_k(-1), cpg(0), sea_pos(-1), 
          sea_pos_k(-1), sea_neg(-1), sea_neg_k(-1), valid(true), s_pos(0), s_neg(0){}
  string to_str(){
    std::stringstream ss;
    ss.setf(ios::fixed, ios::floatfield);
    ss << "Value: valid=" << valid 
       << setprecision(2)
       << " len=" << len
       << " Tm=" << tm_m << " Tmdiff=" << abs(tm_m-tm_u)
       << " gc=" << gc << " cpg=" <<cpg << " sa=" << sa <<" sa_k=" <<sa_k
       << " sea=" << sea_pos <<" sea_k=" <<sea_pos_k
       << " s_pos=" <<s_pos << " s_neg=" <<s_neg;
      //       << "SCORE: " << primer_cond.length.score(len) << " " << primer_cond.gc.score(gc) << " " << primer_cond.tm.score(tm) <<" " << primer_cond.sa.score(sa) <<" " << primer_cond.sea.score(sea_pos) << " " << primer_cond.sea.score(sea_neg);
    return ss.str();
  }
  int len; // for debug
  float tm_m;
  float tm_u;
  float gc;
  int sa;
  int sa_k;
  int cpg;
  int sea_pos;
  int sea_pos_k;
  int sea_neg;
  int sea_neg_k;
  bool valid;
  float s_pos;
  float s_neg;
};

template<typename T>
class array2{
public:
  array2(int n, int m, const T& value=T()):n_(n), m_(m), vec_(n*m, value){};
  inline T& get(int i,int j){return vec_[i*m_+j];}
private:
  int n_;
  int m_;
  vector<T> vec_;
};

class array_sa: public array2<int>{
public:
  array_sa(int n, int m, const int value):array2<int>(n,m,value){}
  inline int& get_se(int start, int end){return get(start, end-start+1);}
};

template<typename T>
class array3{
public:
  array3(int n, int m, int l, const T& value=T()):n_(n), m_(m), l_(l), lm_(l*m), vec_(n*m*l, value){};
  inline T& get(int i, int j, int k){return vec_[lm_*i + l_*j + k];}
private:
  int n_;
  int m_;
  int l_;
  int lm_;
  vector<T> vec_;
};


class PrimerPairResult{
public:
  PrimerPairResult(const string& input, float score, int i, int n, int j, int m, int pa, int pa_k, int pea, int pea_k)
    : score_(score), i_(i), n_(n), j_(j), m_(m), pa_(pa), pa_k_(pa_k), pea_(pea), pea_k_(pea_k) {
    fw_ = input.substr(i,n);
    rv_ = reverse_complement(input.substr(j,m));
    replace(fw_.begin(), fw_.end(), 'C', 'Y'); // Y = C or T
    replace(rv_.begin(), rv_.end(), 'G', 'R'); // R = G or A
  }
  bool operator<(const PrimerPairResult& rhs)const { return score_ < rhs.score_; }
  float score_;
  int i_, n_, j_, m_, pa_, pa_k_, pea_, pea_k_;
  string fw_, rv_;
};

class Results : public std::set<PrimerPairResult>{
public:
  Results(int size) : size_(size), lowest_(1000000){}

  bool add(const PrimerPairResult & ppr){

    if(set<PrimerPairResult>::size() < size_){
      set<PrimerPairResult>::insert(ppr);
      lowest_ = last()->score_;
      return true;
    }
    else if(ppr.score_ < lowest_){
      set<PrimerPairResult>::insert(ppr);
      set<PrimerPairResult>::erase(last());
      lowest_ = last()->score_;
      return true;
    }
    return false;
  }
  std::set<PrimerPairResult>::iterator last(){
    auto e = set<PrimerPairResult>::end();
    return --e;
  }
  float lowest(){return lowest_;}
private:
  int size_;
  float lowest_;
};

class AnnealStore{
public:
  AnnealStore():value(-1), index(0){}
  void store(float v, int i){
    if(v > value){
      value = v;
      index = i;
    }
  }
  float value;
  int index;
};
/*
usage
int pos=0, pos_k=0;
MAX_K(pos, pos_k, 0, length, seqint.something(k))

def max_k(function, range_):
    ret = None
    ret_k = None
    for k in range_:
        v = function(k)
        if (not ret) or ret<v:
            ret = v
            ret_k = k
    return ret,ret_k
*/
#define MAX_K(RET, RETK, FROM, TO, CALC) {for(int k=(FROM);k<=(TO); k++){int value=(CALC); if(value>RET){RET=value; RETK=k;}} }


void bisearch(const char* input, ostream& output,
                  int product_len_min, int product_len_max,
                  int primer_len_min, int primer_len_max,
                  PCRCondition cond,
                  float max_tm_diff,
                  float max_met_tm_diff,
                  float score_threshold,
                  int max_results)
{
  PrimerCondition primer_cond(primer_len_min, primer_len_max, true);
  Condition product_cond(0, 0, product_len_min, product_len_max);
  MeltTemp mt(cond);

  const string genomic(input);
  string seq;
  bisulfite<true>(seq, genomic);

  const int length = seq.length();
  Seqint seqint(seq);

  //output << "# bisearch length=" << length << " primer_len=" << primer_len_max << endl;

  // TODO: I dont need a field for primer_len_max<minimum. DROP IT OUT.
  array2<Value> cs(length, primer_len_max);
  
  //output << "# Tm, length, GC calculation" << endl;

  // TM of fw primer and rv primer is identical.
  // TODO. dynamic programming for TM calculation.
  for(int i=0; i<length; i++){
    int cpg=0;
    int gc=0;
    if(seqint.is_cpg(i)) cpg++;
    if(seqint.is_gc(i)) gc++;
    float dh = 0;
    float dg = 0; 
    float dh_u = 0;
    float dg_u = 0;

    int n=2;
    for(; n<=min(length-i, primer_len_min-1); n++){
      assert( i+n <= length );
      assert( ! primer_cond.length.within(n) );
      const int j = i+n-1;

      // calc state of the primer seq[i:j]
      if(seqint.is_cpg(j)) cpg++;
      if(seqint.is_gc(j)) gc++;
      dh += seqint.tm_dh_m[j-1];
      dg += seqint.tm_dg_m[j-1];
      dh_u += seqint.tm_dh_u[j-1];
      dg_u += seqint.tm_dg_u[j-1];

      Value& v = cs.get(i,n);
      v.valid = false;
    }
    for(; n<=min(length-i, primer_len_max); n++){
      assert( i+n <= length );
      assert( primer_cond.length.within(n) );
      const int j = i+n-1;

      // calc state of the primer seq[i:j]
      if(seqint.is_cpg(j)) cpg++;
      if(seqint.is_gc(j)) gc++;
      dh += seqint.tm_dh_m[j-1];
      dg += seqint.tm_dg_m[j-1];
      dh_u += seqint.tm_dh_u[j-1];
      dg_u += seqint.tm_dg_u[j-1];

      Value& v = cs.get(i,n);
      
      // num of CpG in primer
      v.gc = 100.0*gc/n;
      v.cpg = cpg;
      if(cpg>1){
        v.valid = false;
        continue;
      }
      // GC ratio;
      // TODO switch according to bsp or not.
      if( ! primer_cond.gc.within(v.gc) ){
        v.valid = false;
        continue;
      }

      // Tm(methylated) conditions
      float tm_m = mt.calc_tm(dh, dg);
      float tm_u = mt.calc_tm(dh_u, dg_u);
      v.tm_m = tm_m;
      v.tm_u = tm_u;
      if( ! (primer_cond.tm.within(tm_m) && primer_cond.tm.within(tm_u)) ){
        v.valid = false;
        continue;
      }
  
      // if no cpg, no tm difference.
      assert(cpg>0 || abs(tm_m-tm_u)<0.001);
      if(cpg>0 && (abs(tm_m-tm_u) > max_met_tm_diff)){
        v.valid = false;
        continue;
      }

      assert(v.valid);
      v.len = n;
    }
  }

  // self annealing calculation
  // sa(i...i+n-1) = max{k from -(n-1) to n-1} sa_k(i...i+n-1)
  // sa_k(i...i+n-1) = k<0 | sa_0(i.....i+n-1+k)
  //                   k=0 | sa_0(i.....i+n-1)
  //                   k>0 | sa_0(i+k...i+n-1)
  // sa first step, store all sa_0 value for sa
  // sa_0(i...i+n-1) = n=1 | score_c(i,i)
  //                   n=2 | 2*score_c(i,i+1)
  //                   n=n | sa_0(i+1,i+n-1-1) + 2*score_c(i,i+n-1)

  //output << "# SA first step" << endl;
  // SA first step
  array_sa sa_0(length, primer_len_max+1, -1);
  for(int end=0; end<length; end++){
    for(int n=1; n<=min(end+1, primer_len_max); n++){
      const int i = end-n+1;
      assert(i+n <= length);

      if(n==1){
        assert(i==end);
        sa_0.get(i,n) = seqint.s_c(i,i);
      }
      else if(n==2){
        sa_0.get(i,n) = 2 * seqint.s_c(i,end);
      }
      else{
        assert(sa_0.get_se(i+1, end-1) >=0);
        sa_0.get(i,n) = sa_0.get_se(i+1, end-1) + 2 * seqint.s_c(i,end);
      }
      assert(sa_0.get(i,n) >=0);
    }
  }
  // SA second step
  //output << "# SA 2nd step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_len_min; n<=min(length-i, primer_len_max); n++){
      assert( i+n <= length );
      Value& v = cs.get(i,n);
      // for debug, comment out
      //if(!v.valid) continue;

      int sa = -1;
      int sa_k = 0;
      // max_k
      MAX_K(sa, sa_k, -(n-1), n-1, (k<0) ? sa_0.get(i,n+k) : sa_0.get(i+k,n-k) );
      assert(sa >= 0);
      v.sa = sa;
      v.sa_k = sa_k;
      if(! primer_cond.sa.within(sa) )
        v.valid = false;
    }
  }

  // Self End Annealing calculation
  // sea(i...i+n-1) = max{k from 0 to n-1} sea_k(i...i+n-1)
  // sea_k(i...i+n-1) = sea_0(i+k...i+n-1)
  // sea_0(i...i+n-1) = n=1 | score_c(i,i)
  //                    n=2 | 2*score_c(i,i+1)
  //                    n=n | 

  // SEA first step
  //output << "# SEA 1st step" << endl;
  array_sa sea_0(length, primer_len_max+1, -1);
  array_sa sea_0_full(length, primer_len_max+1, 0);

  // TODO: combine SA and SEA loop.
  for(int end=0; end<length; end++){
    for(int n=1; n<=min(end+1, primer_len_max); n++){
      const int i = end-n+1;
      assert(0<= i && i<= end);
      assert(i+n <= length);
      // fill sea_0.get(i,n)
      if(n==1){
        assert(i==end);
        sea_0.get(i,n) = seqint.s_c(i,i);
        sea_0_full.get(i,n) = true;
      }
      else if(n==2){
        int s = seqint.s_c(i,end);
        sea_0.get(i,n) = 2 * s;
        sea_0_full.get(i,n) = int(s!=0);
      }
      else{
        int base = sea_0.get_se(i+1, end-1);
        assert(base>=0);
        int s = seqint.s_c(i,end);
        bool full = base>0 && s>0 && sea_0_full.get_se(i+1,end-1);
        sea_0_full.get(i,n) = int(full);
        if(full)
          sea_0.get(i,n) = base + 2*s;
        else if(s!=0)
          sea_0.get(i,n) = base + s;
        else
          sea_0.get(i,n) = 0;
      }
    }
  }

  // SEA second step
  //output << "# SEA 2nd step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_len_min; n<=min(length-i, primer_len_max); n++){
      assert( i+n <= length );
      Value& v = cs.get(i,n);
      if(!v.valid)
        continue;

      int pos = -1;
      int pos_k = 0;
      MAX_K(pos, pos_k, 0, n-1, sea_0.get(i+k,n-k) );
      int neg = -1;
      int neg_k = 0;
      MAX_K(neg, neg_k, -(n-1), 0, sea_0.get(i, n+k) );
      assert(neg >= 0 && pos>=0);
      
      v.sea_pos = pos;
      v.sea_pos_k = pos_k;
      v.sea_neg = neg;
      v.sea_neg_k = neg_k;
      
      if(!primer_cond.sea.within(pos) || !primer_cond.sea.within(neg)){
        v.valid = false;
        continue;
      }

      float c = primer_cond.length.score(n)
        + primer_cond.gc.score(v.gc)
        + primer_cond.sa.score(v.sa)
        + primer_cond.tm.score(v.tm_m);
      v.s_pos = c + primer_cond.sea.score(v.sea_pos);
      v.s_neg = c + primer_cond.sea.score(v.sea_neg);

      if( v.s_pos > score_threshold ||  v.s_neg > score_threshold )
        v.valid = false;
    }
  }


  if(debug){
    int c = 0;
    for(int i=0; i<length; i++){
      for(int n=primer_len_min; n<min(length-i, primer_len_max); n++){
        Value& v = cs.get(i,n);
        if(v.valid){
          output << seq.substr(i,n) << " " << v.to_str() << endl;
      
          c++;}
      }
    }
    //output << "# Valid Primers: " << c << endl;
  }

  // Pair Annealing
  //output << "# PA 1st step" << endl;

  array3<int> pa_0(length, length, primer_len_max+1, -1);
  array3<int> pea_0_l(length, length, primer_len_max+1, -1);
  array3<int> pea_0_r(length, length, primer_len_max+1, -1);
  array3<int> pea_0_l_full(length, length, primer_len_max+1, 1);

  /*
    pa0(i,j,k) means, anneal score of seq[i,k] and reverse_complement(seq[j,k])

    i         i+k-1          j       j+k-1
    |----------->------------<---------|
    <-----k----->            <----k---->

    seq[i,k] means seq.substr(i,k), [seq[x] for x in i <= x < i+k]

    k=1
     seq[j,1]     |
     seq[i,1]     |

    1 < k <= primer_len_max
     seq[j,k]   <----|-|
     seq[i,k]   |---->->

     seq[j,k-1] <----|
     seq[i,k-1] |---->

  
    pea0_l means, end score of 5' side of seq[i,i+k], termed left side.
    pea0_r means, end score of 5' side of reverse_complement(seq[j,j+k])), termed right side.

     seq[j,k]   <-------------|-|
     seq[i,k]   |------------->->
     seq[j,k-1] <-------------|
     seq[i,k-1] |------------->

     if s>0
      rt              |------|-|
     else
      rt                       |

     if full(i,j,k-1)
      lt       |-------------|-|
     else
      lt       |---------|      
  */

  // fill pa0
  for(int i=0; i<length; i++){
    for(int j=i; j<length; j++){
      const int primer_len = min(length-j+1, primer_len_max);
      
      // fill pa_0
      pa_0.get(i,j,1) = seqint.s_c_c(i,j);
      for(int k=2; k<=primer_len; k++){
          pa_0.get(i,j,k) = pa_0.get(i,j,k-1) + seqint.s_c_c(i+k-1, j+k-1);
      }

      // fill pea
      // k=1
      int s = seqint.s_c_c(i,j);
      pea_0_l_full.get(i,j,1) = static_cast<int>(s>0);
      pea_0_l.get(i,j,1) = s;
      pea_0_r.get(i,j,1) = s;
      assert(s>=0);

      for(int k=2; k<=primer_len; k++){
        int lt = pea_0_l.get(i,j,k-1);
        int rt = pea_0_r.get(i,j,k-1);
        
        //cout << lt << " " << rt << " " << i << " " << j << " " << k << " " << endl;;
        assert( lt>=0 && rt>=0 );

        bool full = pea_0_l_full.get(i,j,k-1);
        int s = seqint.s_c_c(i+k-1, j+k-1);

        // for left side, add s if full.
        pea_0_l_full.get(i,j,k) = static_cast<int>(full && (s>0));
        pea_0_l.get(i,j,k) = lt + (full ? s : 0);

        // for right side, reset if s==0
        pea_0_r.get(i,j,k) = (s==0)? 0 : rt+s;       
      }

    }
  }

  // Pair End Annealing

  /*
   pair annealing values for each pair 

   fw = seq.substr(i,n)
   rv = reverse_complement(seq.substr(j,m)

   n >= m

  (k=-(m-1)) nn=1
   rv <------|
   fw        |---------->
  
  (A) nn=m-k
   rv   <------|
   fw        |---------->
  
  (k=0) nn=m
   rv        <------|
   fw        |---------->

  (B) nn=m
   rv          <------|
   fw        |---------->

  (k=n-m) nn=m
   rv            <------|
   fw        |---------->

  (C) nn=n-k
   rv               <------|
   fw        |---------->

  (k=n-1) nn=1
   rv                   <------|
   fw        |---------->

  (A) -(m-1) < k < 0
   nn=m-k
   seq[j,m]  <------|
   seq[i,n]       |---------->
   -> pa(i, j-k, nn)

  (B) 0 < k < n-m
   nn=m
   seq[j,m]         <------|
   seq[i,n]       |---------->
   -> pa(i+k, j, nn)

  (C) n-m < k < n-1
   nn=n-k
   seq[j,m]              <------|
   seq[i,n]       |---------->
   -> pa(i+k, j, nn)


  n <= m

  (k=-(m-1)) nn=1
   rv <----------|
   fw            |------>
  
  (A) nn=m-k
   rv   <----------|
   fw            |------>
  
  (k=n-m) nn=n
   rv        <----------|
   fw            |------>

  (B) nn=n
   rv          <----------|
   fw            |------>

  (k=0) nn=n
   rv            <----------|
   fw            |------>

  (C) nn=n-k
   rv                <----------|
   fw            |------>

  (k=n-1) nn=1
   rv                   <----------|
   fw            |------>

  (A) -(m-1) < k < n-m
   nn=m-k
   seq[j,m]   <----------|
   seq[i,n]            |------>
   -> pa(i, j-k, nn)

  (B) n-m < k < 0
   nn=n
   seq[j,m]        <----------|
   seq[i,n]          |------>
   -> pa(i, j-k, nn)

  (C) 0 < k < n-1
   nn=n-k
   seq[j,m]            <----------|
   seq[i,n]       |------>
   -> pa(i+k, j, nn)
  */

  Results results(max_results);

  //output << "# PA 2st step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_len_min; n<=min(length-i, primer_len_max); n++){
      assert( i+n <= length );
      Value& fw = cs.get(i,n);
      if(!fw.valid) continue;

      for(int j=i+product_len_min; j<min(i+product_len_max,length-1); j++){
        for(int m=primer_len_min; m<=min(length-j,primer_len_max); m++){
          const int product_len = j+m-i+1;
          if( !product_cond.within(product_len) )
            continue;

          // rv primer is reverse_complement(seq[j:j+m])
          Value& rv = cs.get(j,m);

          if(!rv.valid) continue;
          if(abs(fw.tm_m-rv.tm_m) > max_tm_diff) continue;

          float score = fw.s_pos + rv.s_neg;
          if( score > score_threshold ) continue;
          if( score > results.lowest()) continue;
          // rank in check
                    
          AnnealStore pa;
          if(n>=m){
            int k = -(m-1);
            for(; k<0; k++)
              pa.store( pa_0.get(i, j-k, m-k), k);
            for(; k<n-m; k++)
              pa.store( pa_0.get(i+k, j, m), k);
            for(; k<n-1; k++ )
              pa.store( pa_0.get(i+k, j, n-k), k);
          }
          else{
            int k = -(m-1);
            for(; k<n-m; k++)
              pa.store( pa_0.get(i, j-k, m-k), k);
            for(; k<0; k++)
              pa.store( pa_0.get(i, j-k, n), k);
            for(; k<n-1; k++ )
              pa.store( pa_0.get(i+k, j, n-k), k);            
          }

          if(! primer_cond.pa.within(pa.value) )
            continue;
          score += primer_cond.pa.score(pa.value);

          AnnealStore pea;
          if(n>=m){
            int k=0;
            for(; k<n-m; k++)
              pea.store( pea_0_l.get(i+k, j, m), k);
            for(; k<n-1; k++ ){
              pea.store( pea_0_l.get(i+k, j, n-k), k);
              pea.store( pea_0_r.get(i+k, j, n-k), k);
            }
          }
          else{
            int k = n-m;
            for(; k<0; k++)
              pea.store( pea_0_r.get(i, j-k, n), k);
            for(; k<n-1; k++ ){
              pea.store( pea_0_r.get(i+k, j, n-k), k);            
              pea.store( pea_0_l.get(i+k, j, n-k), k);            
            }
          }

          if(! primer_cond.pea.within(pea.value) )
            continue;
          score += primer_cond.pea.score(pea.value);

          // cs.get(i,n), cs.get(j,m) for sa, sea
          results.add(PrimerPairResult(seq, score, i, n, j, m, pa.value, pa.index, pea.value, pea.index));
        }
      }
    }
  }
  

  output << "> general" << endl;
  output << "sequence: " << genomic << endl;
  output << endl;
  output << "> bs_pcr" << endl;

  int c = 0;
  for(auto i : results){
    c++;
    output << "// rank=" << c << " score=" << i.score_ << " pea=" << i.pea_ << " pea_k=" << i.pea_k_ << endl;
    output << "//   " << left << setfill(' ') << setw(35) << ("5'-"+i.fw_+"-3'") << ": " << cs.get(i.i_,i.n_).to_str() << endl;
    output << "//   " << left << setfill(' ') << setw(35) << ("5'-"+i.rv_+"-3'") << ": " << cs.get(i.j_,i.m_).to_str() << endl;
    output << "Bi-" << right << setfill('0') << setw(3) << c << ": " << i.fw_ << ", " << i.rv_ << endl;
  }
  output << "// Total: " << c << " Results." << endl;
}



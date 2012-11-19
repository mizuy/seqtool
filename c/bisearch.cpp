#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <cctype>
#include <cassert>
#include <sstream>

#include "bisearch.h"
#include "nucleotide.h"

using namespace std;

// ATGC
const float nnt_dh[4][4] = {
  {9.1, 8.6, 7.8, 6.5}, //AA, AT, AG, AC
  {6.0, 9.1, 5.8, 5.6},
  {5.6, 6.5, 11.0, 11.1},
  {5.8, 7.8, 11.9, 11.0},
};
const float nnt_dg[4][4] = {
  {1.55, 1.25, 1.45, 1.40},
  {0.85, 1.55, 1.15, 1.15},
  {1.15, 1.40, 2.30, 2.70},
  {1.15, 1.45, 3.05, 2.30},
};

const int A = 0;
const int T = 1;
const int G = 2;
const int C = 3;

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
  PrimerCondition():
    length(0.5, 23., 20, 35),
    gc(1.0, 50., 30, 70),
    gc_bsp(1.0, 30., 0, 60),
    tm(1.0, 60., 45, 70),
    sa(0.1, 0., 0, 20),
    sea(0.2, 0., 0, 10),
    pa(0.1, 0., 0, 20),
    pea(0.2, 0., 0, 10){
  }
  Condition length;
  Condition gc;
  Condition gc_bsp;
  Condition tm;
  Condition sa;
  Condition sea;
  Condition pa;
  Condition pea;
};
PrimerCondition primer_cond;

Condition product_cond(0, 0, 100, 400);

// setting values;
const float na_conc = 0.033;
const float mg_conc = 0.0065;
const float primer_conc = 0.0000005;

const float c_salt = na_conc + 4*pow(mg_conc,0.5f);

const float cc_salt = 16.6 * log10(c_salt / (1+0.7*c_salt)) - 269.3;
const float cc_primer = 10.2 - 0.001*1.987*298.2*log(primer_conc);

const int min_product = 100;
const int max_product = 600;
const float threshold = 200.0;

const float max_met_tm_diff = 2.5;
const float max_tm_diff = 8.0;

float calc_tm(float enthalpy, float free_energy){
  return cc_salt + 298.2 * (10.0+enthalpy) / (enthalpy-free_energy + cc_primer);
}
float seq_tm(const vector<char>& seq){
  float p=0; //enthalpy
  float q=0; //free energy
  for(int i=0; i<seq.size()-1; i++){
    p += nnt_dh[seq[i]][seq[i+1]];
    q += nnt_dg[seq[i]][seq[i+1]];
  }
  return calc_tm(p,q);
}

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
    tm_dh.resize(length-1);
    tm_dg.resize(length-1);
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
      tm_dh[i] = nnt_dh[p][q];
      tm_dg[i] = nnt_dg[p][q];
      if(p==C) p=T;
      if(q==C) q=T;
      tm_dh_u[i] = nnt_dh[p][q];
      tm_dg_u[i] = nnt_dg[p][q];
    }
    // Tm of seq[p:q] is calc_tm(sum(tm_dh[p:q-1]), sum(tm_dg[p:q-1]))
  }
  inline bool is_cpg(int i){
    return (i+1<seqint_.size()) ? (seqint_[i]==C and seqint_[i+1]==G) : false;
  }
  inline bool is_gc(int i){
    return seqint_[i]==C or seqint_[i]==G;
  }

  inline int s_c(int p, int q){
    return anneals[seqint_[p]][seqint_[q]];
  }
  inline int s_c_c(int p, int q){
    return anneals_c[seqint_[p]][seqint_[q]];
  }
  inline float tm_m(int i, int j){
    float e=0; //enthalpy
    float f=0; //free energy
    for(int i=0; i<j; i++){
      e += tm_dh[i];
      f += tm_dg[i];
    }
    return calc_tm(e,f);
  }
  inline float tm_u(int i, int j){
    float e=0; //enthalpy
    float f=0; //free energy
    for(int i=0; i<j; i++){
      e += tm_dh_u[i];
      f += tm_dg_u[i];
    }
    return calc_tm(e,f);
  }


  vector<float> tm_dh;
  vector<float> tm_dg;
  vector<float> tm_dh_u;
  vector<float> tm_dg_u;
private:
  vector<char> seqint_;
};


class Value{
public:
  Value():tm(-1), tm_u(-1), gc(-1), sa(-1), sa_k(-1), cpg(0), sea_pos(-1), 
          sea_pos_k(-1), sea_neg(-1), sea_neg_k(-1), valid(true), s_pos(0), s_neg(0){}
  string to_str(){
    std::stringstream ss;
    ss << "Value: valid=" << valid << " Tm=" << tm << " Tmdiff=" << abs(tm-tm_u)
       << " gc=" << gc << " cpg=" <<cpg << " sa=" << sa <<" sa_k=" <<sa_k
       << " sea=" << sea_pos <<" sea_k=" <<sea_pos_k
       << " s_pos=" <<s_pos << " s_neg=" <<s_neg;
    return ss.str();
  }
  float tm;
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
class array{
public:
  array(int n, int m, const T& value=T()):n_(n), m_(m), vec_(n*m, value){};
  inline T& get(int i,int j){return vec_[i*m_+j];}
private:
  int n_;
  int m_;
  vector<T> vec_;
};

class array_sa: public array<int>{
public:
  array_sa(int n, int m, const int value):array<int>(n,m,value){}
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
  }
  float score_;
  int i_, n_, j_, m_, pa_, pa_k_, pea_, pea_k_;
  string fw_, rv_;
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


void bisearch(const char* input){
  const string seq(input);
  const int length = seq.length();
  Seqint seqint(seq);

  const int primer_len = primer_cond.length.maximum;

  cout << "bisearch length=" << length << " primer_len=" << primer_len << endl;

  // TODO: I dont need a field for primer_len<minimum. DROP IT OUT.
  array<Value> cs(length, primer_len);
  
  cout << "Tm, length, GC calculation" << endl;

  // TODO. I need TM for reverse primer.
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
    for(; n<=min(length-i, (int)primer_cond.length.minimum-1); n++){
      assert( i+n <= length );
      assert( ! primer_cond.length.within(n) );
      const int j = i+n-1;

      // calc state of the primer seq[i:j]
      if(seqint.is_cpg(j)) cpg++;
      if(seqint.is_gc(j)) gc++;
      dh += seqint.tm_dh[j-1];
      dg += seqint.tm_dg[j-1];
      dh_u += seqint.tm_dh_u[j-1];
      dg_u += seqint.tm_dg_u[j-1];

      Value& v = cs.get(i,n);
      v.valid = false;
    }
    for(; n<=min(length-i, primer_len); n++){
      assert( i+n <= length );
      assert( primer_cond.length.within(n) );
      const int j = i+n-1;

      // calc state of the primer seq[i:j]
      if(seqint.is_cpg(j)) cpg++;
      if(seqint.is_gc(j)) gc++;
      dh += seqint.tm_dh[j-1];
      dg += seqint.tm_dg[j-1];
      dh_u += seqint.tm_dh_u[j-1];
      dg_u += seqint.tm_dg_u[j-1];

      Value& v = cs.get(i,n);
      //cout << i<<" "<<n<<" "<<seq.substr(i,n) << endl;
      
      // num of CpG in primer
      v.gc = 100.0*gc/n;
      v.cpg = cpg;
      if(cpg>1){
        v.valid = false;
        continue;
      }
      // GC ratio;
      // TODO switch according to bsp or not.
      if( ! primer_cond.gc_bsp.within(v.gc) ){
        //cout << "gc cond" << gc << "/" << n << endl;
        v.valid = false;
        continue;
      }

      // Tm(methylated) conditions
      float tm = calc_tm(dh, dg);
      float tm_u = calc_tm(dh_u, dg_u);
      v.tm = tm;
      v.tm_u = tm_u;
      //cout << "tm: " << tm << " tmu:" << tm_u << endl;
      if( ! (primer_cond.tm.within(tm) && primer_cond.tm.within(tm_u)) ){
        v.valid = false;
        //cout << "tm cond " << tm << " " << tm_u << endl;
        continue;
      }
  
      // if no cpg, no tm difference.
      assert(cpg>0 || abs(tm-tm_u)<0.001);
      if(cpg>0 && (abs(tm-tm_u) > max_met_tm_diff)){
        v.valid = false;
        //cout << "tm diff cond " << abs(tm-tm_u)<< " > " << max_met_tm_diff << endl;
        continue;
      }

      assert(v.valid);
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

  cout << "SA first step" << endl;
  // SA first step
  array_sa sa_0(length, primer_len+1, -1);
  for(int end=0; end<length; end++){
    for(int n=1; n<=min(end+1, primer_len); n++){
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
  cout << "SA 2nd step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_cond.length.minimum; n<=min(length-i, primer_len); n++){
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
  cout << "SEA 1st step" << endl;
  array_sa sea_0(length, primer_len+1, -1);
  array_sa sea_0_full(length, primer_len+1, 0);

  // TODO: combine SA and SEA loop.
  for(int end=0; end<length; end++){
    for(int n=1; n<=min(end+1, primer_len); n++){
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
  cout << "SEA 2nd step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_cond.length.minimum; n<=min(length-i, primer_len); n++){
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
        + primer_cond.tm.score(v.tm);
      v.s_pos = c + primer_cond.sea.score(v.sea_pos);
      v.s_neg = c + primer_cond.sea.score(v.sea_neg);

      if( v.s_pos > threshold ||  v.s_neg > threshold )
        v.valid = false;
    }
  }


  {
    int c = 0;
    for(int i=0; i<length; i++){
      for(int n=primer_cond.length.minimum; n<min(length-i, primer_len); n++){
        Value& v = cs.get(i,n);
        if(v.valid){
          cout << seq.substr(i,n) << " " << v.to_str() << endl;
      
          c++;}
      }
    }
    cout << "Valid Primers: " << c << endl;
  }

  // Pair Annealing
  // pair annealing values for each pair w=target[i:i+n], v=target[j:j+m]
  // pa(w,v) = max{k from -(n-1) to m-1} 
  //                 pa_0[1,1,n+k] if k<0, n+k<m
  cout << "PA 1st step" << endl;

  array3<int> pa_0(length, length, primer_len+1, -1);
  array3<int> pea_0_l(length, length, primer_len+1, -1);
  array3<int> pea_0_r(length, length, primer_len+1, -1);
  array3<int> pea_0_l_full(length, length, primer_len+1, 1);

  // fill pa0
  for(int i=0; i<length; i++){
    for(int j=0; j<length; j++){
      for(int k=0; k<primer_len+1; k++){ // TODO. fix.
        // fill pa0
        if(k==0)
          pa_0.get(i,j,k) = seqint.s_c_c(i,j);
        else
          pa_0.get(i,j,k) = pa_0.get(i,j,k-1) + seqint.s_c_c(i+k-1, j+k-1);
            
        // fill pea0
        if(k==0)
          continue;
        else if(k==1){
          int s = seqint.s_c_c(i,j);
          pea_0_l_full.get(i,j,k) = int(s!=0);
          pea_0_l.get(i,j,k) = s;
          pea_0_r.get(i,j,k) = s;
        }
        else{
          int lt = pea_0_l.get(i,j,k-1);
          int rt = pea_0_r.get(i,j,k-1);
          //cout << length << " " << primer_len+1 << " " << i << " " << j << " " << k << " lt:" << lt << " rt:" << rt << endl;
          //assert( lt>=0 && rt>=0 );
          // OUTRANGE
          int s = seqint.s_c_c(i+k-1, j+k-1);
          bool full = (s!=0) && pea_0_l_full.get(i,j,k-1);
          int new_lt = full ? (lt + s) : lt;
          int new_rt = (s!=0) ? (rt + s) : 0;
          pea_0_l_full.get(i,j,k) = int(full);
          //if(new_lt > 0 || new_rt > 0)
          //  cout << "s: " << s << " newlt: " << new_lt << " newrt: " << new_rt << endl;
          pea_0_l.get(i,j,k) = new_lt;
          pea_0_r.get(i,j,k) = new_rt;
        }
      }
    }
  }

  vector<PrimerPairResult> result;

  cout << "PA 2st step" << endl;
  for(int i=0; i<length; i++){
    for(int n=2; n<=min(length-i, primer_len); n++){

      assert( i+n <= length );
      Value& fw = cs.get(i,n);
      if(!fw.valid) continue;

      for(int j=i+product_cond.minimum-primer_len; j<=min(i+(int)product_cond.maximum-primer_len,length); j++){
        for(int m=2; m<=min(length-j,primer_len); m++){
          if( !product_cond.within(j+m-i+1) )
            continue;
          Value& rv = cs.get(j,m);

          if(!rv.valid)
            continue;

          //cout << " i " << i << " n " << n << " j " << j << " m " << m << endl;

          if(abs(fw.tm-rv.tm) > max_tm_diff){
            //cout << "max_tm_diff " << abs(fw.tm-rv.tm) << endl;
            continue;
          }
          float score = fw.s_pos + rv.s_neg;
          if( score > threshold ){
            //cout << "score " << score << endl;
            continue;
          }
                    
          int pa = 0;
          int pa_k = -1;
          int pea = 0;
          int pea_k = -1;
          // MAX_K
          for(int k=-(n-1); k<=m-1; k++){
            int pa_;
            int pea_;
            if(n<=m){
              if (k<=0)
                pa_ = pa_0.get(i+max(0,-k), j+max(0,k), n+k);
              else if(k<=m-n)
                pa_ = pa_0.get(i+max(0,-k), j+max(0,k), n);
              else
                pa_ = pa_0.get(i+max(0,-k), j+max(0,k), m-k);
              
              if(k<=0){
                int l = pea_0_l.get(i+max(0,-k), j+max(0,k), n+k);
                int r = pea_0_r.get(i+max(0,-k), j+max(0,k), n+k);
                pea_ = max(l,r);
              }
              else if(k<=m-n){
                pea_ =  pea_0_r.get(i+max(0,-k), j+max(0,k), n);
              }
            }
            else{
              if(k<=m-n)
                pa_ = pa_0.get(i+max(0,-k), j+max(0,k), n+k);
              else if(k<=0)
                pa_ = pa_0.get(i+max(0,-k), j+max(0,k), m);
              else
                pa_ = pa_0.get(i+max(0,-k), j+max(0,k), m-k);

              if(k<=m-n){
                int l = pea_0_l.get(i+max(0,-k), j+max(0,k), n+k);
                int r = pea_0_r.get(i+max(0,-k), j+max(0,k), n+k);
                pea_ = max(l,r);
              }
              else if(k<=0)
                pea_ = pea_0_l.get(i+max(0,-k), j+max(0,k), m);
              else
                pea_ = 0;
            }
            if(pa_ > pa){
              pa = pa_;
              pa_k = k;
            }
            if(pea_ > pea){
              pea = pea_;
              pea_k = k;
            }

          }

          if(! primer_cond.pa.within(pa) )
            continue;
          if(! primer_cond.pea.within(pea) )
            continue;
          score += primer_cond.pa.score(pa);
          score += primer_cond.pea.score(pea);

          // cs.get(i,n), cs.get(j,m) for sa, sea
          result.push_back(PrimerPairResult(seq, score, i, n, j, m, pa, pa_k, pea, pea_k));

          //ppr = PrimerPairResult(target, primerp, score,i,n,j,m,pa,pa_k,pea,pea_k,mode_bsp)
          //result[threadno].add(ppr)
          }
        }

    }
  }

  int c = 0;
  for(vector<PrimerPairResult>::iterator i=result.begin(); i!=result.end(); i++){
    cout << "primer !!: score=" << i->score_ << " fw=" << i->fw_ << " rv="<< i->rv_ << endl;
    cout << "   " << cs.get(i->i_,i->n_).to_str() << endl;
    cout << "   " << cs.get(i->j_,i->m_).to_str() << endl;
    c++;
  }
  cout << "Total: " << c << " Results." << endl;
}


const char * test = "TCCTCTTCTGGGAGTAGGCAGAAGACTCCCGGGAGGAGAGGCGAACAGCGGACGCCAATTCTTTTGAAAGCACTGTGTTCCTTAGCACCGCGGGTCGCTACGGGCCTCTTGCTGTCGCGGGATTTCGGTCCACCTTCCGATTGGGCCGCCGCATCCCGGATCAGATTTCGCGGGCGACCCACGGAACCCGCGGAGCCGGGACGTGAAAGGTTAGAAGGTTTCCCGTTCCCATCAAGCCCTAGGGCTCCTCGTGGCTGCTGGGAGTTGTAGTCTGAACGCTTCTATCTTGGCGAGAAGCGCCTACGCTCCCCCTACCGAGTCCCGCGGTAATTCTTAAAGCACCTGCACCGCCCCCCCGCCGCCTGCAGAGGGCGCAGCAGGTCTTGCACCTCTTCTGCATCTCA";
//const char * test = "TCCTCTTCTGGGAGTAGGCAGAAGACTCCCGGGAGGAGAGGCGAACAGCGGACGCCAATTCTTTTGAAAGCACTGTGTTCCTTAGCACCGCGGGTCGCTACGGGCCTCTTGCTGTCGCGGGATTTCGGTCCACCTTCCGATTGGGCCGCCGCATCCCGGATCAGATTTCGCGGGCGACCCACGGAACCCGCGGAGCCGGGACGTGAAAGGTTAGAAGGTTTCCCGTTCCCATCAAGCCCTAGGGCTCCTCGTGGCTGCTGGGAGTTGTAGTCTGAACGCTTCTATCTTGGCGAGAAGCGCCTACGCTCCCCCTACCGAGTCCCGCGGTAATTCTTAAAGCACCTGCACCGCCCCCCCGCCGCCTGCAGAGGGCGCAGCAGGTCTTGCACCTCTTCTGCATCTCATTCTCCAGGCTTCAGACCTGTCTCCCTCATTCAAAAAATATTTATTATCGAGCTCTTACTTGCTACCCAGCACTGATATAGGCACTCAGGAATACAACAATGAATAAGATAGTAGAAAAATTCTATATCCTCATAAGGCTTACGTTTCCATGTACTGAAAGCAATGAACAAATAAATCTTATCAGAGTGATAAGGGTTGTGAAGGAGATTAAATAAGATGGTGTGATATAAAGTATCTGGGAGAAAACGTTAGGGTGTGATATTACGGAAAGCCTTCCTAAAAAATGACATTTTAACTGATGAGAAGAAAGGATCCAGCTGAGAGCAAACGCAAAAGCTTTCTTCCTTCCACCCTTCATATTTGACACAATGCAGGATTCCTCCAAAATGATTTCCACCAATTCTGCCCTCACAGCTCTGGCTTGCAGAATTTTCCACCCCAAAATGTTAGTATCTACGGCACCAGGTCGGCGAGAATCCTGACTCTGCACCCTCCTCCCCAACTCCATTTCCTTTGCTTCCTCCGGCAGGCGGATTACTTGCCCTTACTTGTCATGGCGACTGTCCAGCTTTGTGCCAGGAGCCTCGCAGGGGTTGATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCG";

int main(){
  string seq(test);
  string b;
  bisulfite<true>(b, seq);
  cout << "bisearch length=" << seq.length() << endl;
  cout << "Sequence: " << seq << endl;
  cout << "Bisulfite treated sequence: " << b << endl;
  bisearch(b.c_str());
};

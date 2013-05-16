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

#include "melt_temp.h"
#include "container.h"
#include "bisearch.h"
#include "nucleotide.h"
#include "primer.h"

const bool debug = true;

using namespace std;
using namespace nucleotide;

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


class Seqint{
public:
  Seqint(const string& seq)
  {
    seqstr_ = seq;
    const int length = seq.length();
    seqint_.resize(length);
    for(int i=0; i<length; i++){
      seqint_[i] = base_int(seq[i]);
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
      tm_dh_m[i] = melt_temp::nnt_dh[p][q];
      tm_dg_m[i] = melt_temp::nnt_dg[p][q];
      if(p==C) p=T;
      if(q==C) q=T;
      tm_dh_u[i] = melt_temp::nnt_dh[p][q];
      tm_dg_u[i] = melt_temp::nnt_dg[p][q];
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

  vector<float> tm_dh_m;
  vector<float> tm_dg_m;
  vector<float> tm_dh_u;
  vector<float> tm_dg_u;
private:
  string seqstr_;
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

class array_sa: public array2<int>{
public:
  array_sa(int n, int m, int value):array2<int>(n,m,value){}
  inline int& get_se(int start, int end){return get(start, end-start+1);}
};

class PrimerPairResult{
public:
  PrimerPairResult(float score, int i, int n, int j, int m, int pa, int pa_k, int pea, int pea_k)
    : score_(score), i_(i), n_(n), j_(j), m_(m), pa_(pa), pa_k_(pa_k), pea_(pea), pea_k_(pea_k) {}
  bool operator<(const PrimerPairResult& rhs)const { return score() > rhs.score(); }
  bool operator>(const PrimerPairResult& rhs)const { return score() < rhs.score(); }
  float score() const{return score_;}

  float score_;
  int i_, n_, j_, m_, pa_, pa_k_, pea_, pea_k_;
};

class ResultChunk{
  typedef set<PrimerPairResult> container;
public:
  ResultChunk() : size_(0), lowest_(1000000){}
  ResultChunk(int size) : size_(size), lowest_(1000000){}

  void resize(int size){ size_ = size;}

  container::iterator begin(){return set_.begin();}
  container::iterator end(){return set_.end();}
  container::iterator last(){return --end();}
  float lowest(){return lowest_;}

  bool add(const PrimerPairResult & ppr){
    if(size_ <= 0){
      return false;
    }
    else if(set_.size() < size_){
      set_.insert(ppr);
      lowest_ = last()->score();
      return true;
    }
    else{
      auto l = last();
      if(ppr > *l){
        set_.erase(l);
        set_.insert(ppr);
        lowest_ = last()->score();
        return true;
      }
    }
    return false;
  }
private:
  container set_;
  int size_;
  float lowest_;
};

class Results{
public:
  typedef vector<ResultChunk> container;
  container::iterator begin(){return list_.begin();}
  container::iterator end(){return list_.end();}

  Results(int length, int chunk, int coverage){
    length_ = length;
    chunk_ = chunk;
    coverage_ = coverage;

    clength_ = length/chunk + ((length%chunk)?1:0);
    list_.resize(clength_);
    for(auto& i : list_){
      i.resize(coverage_);
    }
  }

  float lowest(int pos){
    int ii = pos/chunk_;
    assert(ii < clength_);
    return list_.at(ii).lowest();
  };
  void add(int pos, const PrimerPairResult & ppr){
    int ii = pos/chunk_;
    assert(ii < clength_);
    list_.at(ii).add(ppr);
  }
private:
  container list_;
  int length_;
  int clength_;
  int chunk_;
  int coverage_;
};


// Store highest score and k value.
class AnnealStore{
public:
  AnnealStore():value(-1), index(0){}
  inline void store(float v, int i){
    if(v > value){
      value = v;
      index = i;
    }
  }
  float value;
  int index;
};

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
  melt_temp::MeltTemp mt(cond);

  const string genomic(input);
  string seq;
  bisulfite<true>(seq, genomic);

  const int length = seq.length();
  Seqint seqint(seq);

  output << "// bisearch length=" << length << " primer_len=" << primer_len_max << endl;

  // TODO: I dont need a field for primer_len_max<minimum. DROP IT OUT.
  array2<Value> cs(length, primer_len_max);
  
  output << "// Tm, length, GC calculation" << endl;

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

  output << "// SA first step" << endl;
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
  output << "// SA 2nd step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_len_min; n<=min(length-i, primer_len_max); n++){
      assert( i+n <= length );
      Value& v = cs.get(i,n);
      // for debug, comment out
      //if(!v.valid) continue;

      AnnealStore sa;
      for(int k=-(n-1);k<=n-1; k++){
        int v = (k<0) ? sa_0.get(i,n+k) : sa_0.get(i+k,n-k);
        sa.store(v, k);
      }
      assert(sa.value >= 0);
      v.sa = sa.value;
      v.sa_k = sa.index;
      if(! primer_cond.sa.within(sa.value) )
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
  output << "// SEA 1st step" << endl;
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
  output << "// SEA 2nd step" << endl;
  for(int i=0; i<length; i++){
    for(int n=primer_len_min; n<=min(length-i, primer_len_max); n++){
      assert( i+n <= length );
      Value& v = cs.get(i,n);
      if(!v.valid)
        continue;

      AnnealStore pos;
      for(int k=0;k<=n-1; k++){
        pos.store(sea_0.get(i+k,n-k), k);
      }
      AnnealStore neg;
      for(int k=-(n-1);k<=0; k++){
        neg.store(sea_0.get(i, n+k), k);
      }
      assert(neg.value >= 0 && pos.value>=0);
      
      v.sea_pos = pos.value;
      v.sea_pos_k = pos.index;
      v.sea_neg = neg.value;
      v.sea_neg_k = neg.index;
      
      if(!primer_cond.sea.within(pos.value) || !primer_cond.sea.within(neg.value)){
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
    output << "// Valid Primers: " << c << endl;
  }

  // TODO: max_results
  Results results(length, length/100, 5);

  output << "// PA' 1st step" << endl;
  array2<int> pag(length, length);
  // fill pa0
  /*
      i     i+n
      |------>
      <---------|
      j        j+m
  
      seq[i,k] means seq.substr(i,k), [seq[x] for x in i <= x < i+k]
      */
  for(int i=0; i<length; i++){
    for(int g=0; g<length; g++){
        const int j = i+g;
        if(!(j<length))
            continue;

        pag.get(i,g) = seqint.s_c_c(i, j);
    }
  }

  output << "// PA' 2nd step" << endl;  
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

          const int position = (i+j+m)/2;
          if(results.lowest(position) < score ) continue;
          //if( score > results.lowest()) continue;
          // rank in check

          /*
            g_u 1
                11
            g_1 111
                111
                111
            g_0 111
                 11
            g_d   1
          */

          const int g_u = (j+m)-(i);   // j-i +m == g+m
          const int g_1 = (j+m)-(i+n);
          const int g_0 = j-i;
          const int g_d = (j)-(i+n-1); // j-i -(n+1) == g-n
          //const int g_min = min(g_0,g_1);
          //const int g_max = max(g_0,g_1);
  
          AnnealStore pa;
          AnnealStore pea;
          for(int gg=g_d; gg<g_u; gg++){
            const int ii_l = max(i,j-gg);
            const int ii_r = min(i+n, j+m-gg);
            assert(ii_l<ii_r);
            const int k = gg-g_0;
            
            int pav = 0;
            for(int ii=ii_l; ii<ii_r; ii++){
              pav += pag.get(ii,gg);
            }
            pa.store(pav, k);

            if(gg<=g_0){
              int peav_l = 0;
              for(int ii=ii_l; ii<ii_r; ii++){
                int v = pag.get(ii,gg);
                if(v==0) break;
                peav_l += v;
              }
              pea.store(peav_l, k);
            }
            if(gg<=g_1){
              int peav_r = 0;
              for(int ii=ii_r-1; ii_l<=ii; ii--){
                int v = pag.get(ii,gg);
                if(v==0) break;
                peav_r += v;
              }
              pea.store(peav_r, k);
            }
          }

          if(debug){
            string fw_str = seq.substr(i,n);
            string rv_str = reverse_complement(seq.substr(j,m));
            auto pa_ = annealing_score(fw_str,rv_str);
            auto pea_ = end_annealing_score(fw_str,rv_str);

            if(pa.value != get<0>(pa_)){
                cout << "//pa mismatch" << endl;
                cout << "//ref: " << get<0>(pa_) << endl;
                print_primer_bar(fw_str,rv_str, get<1>(pa_));
                cout << "//now: " << pa.value << endl;
                print_primer_bar(fw_str,rv_str, pa.index);
            }
            if(pea.value != get<0>(pea_)){
                cout << "//pea mismatch" << endl;
                cout << "//ref: " << get<0>(pea_) << endl;
                print_primer_bar(fw_str,rv_str, get<1>(pea_));
                cout << "//now: " << pea.value << endl;
                print_primer_bar(fw_str,rv_str, pea.index);
            }
          }

          if(! primer_cond.pa.within(pa.value) )
            continue;
          score += primer_cond.pa.score(pa.value);

          if(! primer_cond.pea.within(pea.value) )
            continue;
          score += primer_cond.pea.score(pea.value);

          // cs.get(i,n), cs.get(j,m) for sa, sea
          results.add(position, PrimerPairResult(score, i, n, j, m, pa.value, pa.index, pea.value, pea.index));

        }
      }
    }
  }

  output << "general:" << endl;
  output << "    sequence: " << genomic << endl;
  output << "    show_bisulfite: True" << endl;
  output << endl;
  output << "bs_pcr:" << endl;

  vector<PrimerPairResult> all;
  for(auto& rc: results){
    for(auto& ppr : rc){
      all.push_back(ppr);
    }
  }

  int count = 0;
  
  sort(all.rbegin(), all.rend());

  for(auto& ppr: all){
    count++;
    string fw_str = seq.substr(ppr.i_,ppr.n_);
    string rv_str = reverse_complement(seq.substr(ppr.j_,ppr.m_));
    replace(fw_str.begin(), fw_str.end(), 'C', 'Y'); // Y = C or T
    replace(rv_str.begin(), rv_str.end(), 'G', 'R'); // R = G or A

    output << "    Bi-" << right << setfill('0') << setw(3) << count << ": " << fw_str << ", " << rv_str << endl;
    output << "        // rank=" << count << " score=" << ppr.score_;
    output << " pea=" << ppr.pea_ << " pea_k=" << ppr.pea_k_ << endl;
    output << "        //   " << left << setfill(' ') << setw(35) << ("5'-"+fw_str+"-3'");
    output << ": " << cs.get(ppr.i_,ppr.n_).to_str() << endl;
    output << "        //   " << left << setfill(' ') << setw(35) << ("5'-"+rv_str+"-3'");
    output << ": " << cs.get(ppr.j_,ppr.m_).to_str() << endl;
  }
  output << "// Total: " << count << " Results." << endl;
}



#include "primer.h"
#include "nucleotide.h"

#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <numeric>
#include <cstdlib>

using namespace std;

const int anneals_[4][4] = {
  {0, 2, 0, 0}, // AT
  {2, 0, 0, 0}, // TA
  {0, 0, 0, 4}, // GC
  {0, 0, 4, 0}, // CG
};

inline int anneals(int l, int r){
    return anneals_[nucleotide::base_int(l)][nucleotide::base_int(r)];
}

inline int anneal_value(vector<int>& vv){
    return accumulate(vv.begin(), vv.end(), 0);
}
template<typename Iterator>
inline int end_anneal_value(Iterator begin, Iterator end){
    int ret = 0;
    for(auto i=begin; i!=end; i++){
        if(*i==0){
            return ret;
        }
        ret += *i;
    }
    return ret;
}
inline int end_anneal_value_l(vector<int>& vv){
    return end_anneal_value(vv.begin(), vv.end());
}
inline int end_anneal_value_r(vector<int>& vv){
    return end_anneal_value(vv.rbegin(), vv.rend());
}
inline int end_anneal_value_lr(vector<int>& vv){
    return max(end_anneal_value_l(vv),end_anneal_value_r(vv));
}

// Store highest score and k value.
class AnnealStore{
public:
  AnnealStore():value(-1), index(0){}
  inline void store(int v, int i){
    if(v > value){
      value = v;
      index = i;
    }
  }
  int value;
  int index;
};

tuple<int,int> annealing_score_(const string& fw, const string& rv, bool end_annealing){
    const string w = fw;
    const string v(rv.rbegin(), rv.rend());
    int n = w.length();
    int m = v.length();

    AnnealStore av;
    AnnealStore eav;

    if(n<=m){
        for(int k=-(n-1); k<=m-1; k++){
            if(k<=0){
                // 5'- w[0]....w[-k]....w[n-1] -3'
                //         3'- v[0].....v[n+k-1]....v[m-1] -5'
                vector<int> vv;
                for(int i=0; i<n+k; i++){
                    vv.push_back(anneals(w[-k+i],v[i]));
                }
                av.store(anneal_value(vv),k);
                eav.store(end_anneal_value_lr(vv),k);
            }
            else if(k<=m-n){
                //         w[0]....w[n-1]
                // v[0]....v[k]....v[k+n-1].....v[m-1]
                vector<int> vv;
                for(int i=0; i<n; i++){
                    vv.push_back(anneals(w[0+i],v[k+i]));
                }
                av.store(anneal_value(vv),k);
                eav.store(end_anneal_value_r(vv),k);
            }
            else{
                //        w[0]...w[m-k-1]....w[n-1]
                // v[0]...v[k]...v[m-1]
                vector<int> vv;
                for(int i=0; i<m-k; i++){
                    vv.push_back(anneals(w[i],v[k+i]));
                }
                av.store(anneal_value(vv),k);
            }
        }
    }
    else{
        for(int k=-(n-1); k<=m-1; k++){
            if(k<=m-n){
                // w[0]....w[-k]....w[n-1]
                //         v[0].....v[n+k-1]....v[m-1]
                vector<int> vv;
                for(int i=0; i<n+k; i++){
                    vv.push_back(anneals(w[-k+i],v[i]));
                }
                av.store(anneal_value(vv),k);
                eav.store(end_anneal_value_lr(vv),k);

            }
            else if(k<=0){
                // w[0]....w[k]....w[m-k-1].....w[n-1]
                //         v[0]....v[m-1]
                vector<int> vv;
                for(int i=0; i<m; i++){
                    vv.push_back(anneals(w[k+i],v[0+i]));
                }
                av.store(anneal_value(vv),k);
                eav.store(end_anneal_value_l(vv),k);

            }
            else{
                //        w[0]...w[m-k-1]....w[n-1]
                // v[0]...v[k]...v[m-1]
                vector<int> vv;
                for(int i=0; i<m-k; i++){
                    vv.push_back(anneals(w[i],v[k+i]));
                }
                av.store(anneal_value(vv),k);
            }
        }
    }

    if(end_annealing)
        return make_pair(eav.value,eav.index);
    else
        return make_pair(av.value,av.index);
}
tuple<int,int> annealing_score(const string& fw, const string& rv){
    return annealing_score_(fw,rv,false);
}
tuple<int,int> end_annealing_score(const string& fw, const string& rv){
    return annealing_score_(fw,rv,true);
}

void print_primer_bar(const std::string& fw, const std::string& rv, int index){
    int gap = abs(index);

    string rrv(rv.rbegin(), rv.rend());

    int n = fw.length();
    int m = rv.length();

    if(index>0){
        for(int i=0; i<gap; i++)
            cout << ' ';
        cout << "5'-";
        cout << fw;
        cout << "-3'" << endl;

        for(int i=0; i<gap; i++)
            cout << ' ';
        cout << "  <";
        for(int i=0; i<min(n,m-gap); i++){
            if(anneals(fw[i],rrv[gap+i])>0)
                cout << "|";
            else
                cout << " ";
        }
        cout << ">" << endl;

        cout << "3'-";
        cout << rrv;
        cout << "-5'" << endl;
    }
    else{
        cout << "5'-";
        cout << fw;
        cout << "-3'" << endl;

        for(int i=0; i<gap; i++)
            cout << ' ';
        cout << "  <";
        for(int i=0; i<min(n-gap,m); i++){
            if(anneals(fw[gap+i],rrv[i])>0)
                cout << "|";
            else
                cout << " ";
        }
        cout << ">" << endl;

        for(int i=0; i<gap; i++)
            cout << ' ';
        cout << "3'-";
        cout << rrv;
        cout << "-5'" << endl;
    }
}

void print_primer_pair(const string& fw, const string& rv, bool end_annealing){
    auto v = annealing_score_(fw,rv,end_annealing);
    int av = get<0>(v);
    int index = get<1>(v);

    cout << "value=" << av << ", index="  << index << endl;
    cout << "fw: " << fw << endl;
    cout << "rv: " << rv << endl;

    print_primer_bar(fw,rv,end_annealing);
}

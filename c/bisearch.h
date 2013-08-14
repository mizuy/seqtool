#ifndef BISEARCH_H
#define BISEARCH_H

#include <ostream>


struct PCRCondition{
  PCRCondition(float na_conc, 
               float k_conc,
               float mg_conc,
               float primer_conc)
  : na_conc(na_conc), k_conc(k_conc), mg_conc(mg_conc), primer_conc(primer_conc) {}
  float na_conc;
  float k_conc;
  float mg_conc;
  float primer_conc;
};


void bisearch(const char* input, std::ostream& output, std::ostream& logging,
                  const PCRCondition& cond,
                  int product_len_min=100, int product_len_max=600,
                  int primer_len_min=20, int primer_len_max=35,
                  float max_tm_diff=8.0f,
                  float max_met_tm_diff = 2.5f,
                  int max_cpg_in_primer=1,
                  float score_threshold=100.0f,
                  int max_results=200);

#endif

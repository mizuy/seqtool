#ifndef BISEARCH_H
#define BISEARCH_H

#include <ostream>


struct PCRCondition{
  PCRCondition(float na_conc=0, 
               float k_conc=0.05f,
               float mg_conc=0.0015,
               float primer_conc=0.000001)
    : na_conc(na_conc), k_conc(k_conc), mg_conc(mg_conc), primer_conc(primer_conc) {}
  float na_conc;
  float k_conc;
  float mg_conc;
  float primer_conc;
};
/*
// Na 33mM, K 50mM, Mg 6.5mM, Primer 0.5uM
const PCRCondition myPCRCondition(0.033, 0.050, 0.0065, 0.0000005);
*/

void bisearch(const char* input, std::ostream& output,
                  int product_len_min=100, int product_len_max=600,
                  int primer_len_min=20, int primer_len_max=35,
                  PCRCondition cond=PCRCondition(),
                  float max_tm_diff=8.0f,
                  float max_met_tm_diff = 2.5f,
                  int max_cpg_in_primer=1,
                  float score_threshold=100.0f,
                  int max_results=200);

#endif

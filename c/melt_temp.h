#ifndef MELT_TEMP_H
#define MELT_TEMP_H

#include <cmath>
#include <vector>
#include "bisearch.h"

namespace melt_temp{
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
	    cc_salt = 16.6 * log10(c_salt / (1.0+0.7*c_salt)) +  3.85;

	    // nnt_dh, nnt_dg validity. NNT of a sequence complemtely equals to that of reverse complement.
	    /*if(debug){
	      for(int i=0; i<4; i++){
	        for(int j=0; j<4; j++){
	          int ii = complement(i);
	          int jj = complement(j);
	          //cout << i << " " << j << " " << ii << " " << jj << endl;
	          assert(nnt_dh[i][j] == nnt_dh[jj][ii]);
	          assert(nnt_dg[i][j] == nnt_dg[jj][ii]);
	        }
	      }
	    }*/
	  }
	  float calc_tm(float sdh, float sdg){
	    float dhp = -1000.0f * (2*c_dhe + sdh);
	    float dgp = -1000.0f * (2*c_dge + c_dgi + sdg);
	    return c_t0 * dhp / (dhp-dgp+cc_primer) + cc_salt - c_k;
	  }
	  float seq_tm(const std::vector<char>& seq){
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
      const float c_k =  273.15f;
	  const float c_t0 = 298.2f;
	  const float c_dhe = 5.0; // delta He = 5 kcal/mol
	  const float c_dge = 1.0; // delta Ge = 1 kcal/mol
	  const float c_dgi = -2.2; // delta Gi = -2.2 kcal/mol
	  PCRCondition cond;
	};

};

#endif

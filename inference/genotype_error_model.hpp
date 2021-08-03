//
//  genotype_error_model.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 17/03/2021.
//
#include <vector>

#ifndef genotype_error_model_hpp
#define genotype_error_model_hpp


extern "C"
{
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}

#include <string>
class GenotypeErrorModel{
    unsigned int states;
    std::string name;
    pll_state_t undefinedState;
    long double seqErrorRate;
    long double ADOErrorRate;
public:
    GenotypeErrorModel(std::string name, long double seqErrorRate,long double ADOErrorRate, int states );
   long double getSeqErrorRate() const;
   long double getADOErrorRate() const;
   GenotypeErrorModel &operator=(const GenotypeErrorModel &original);
   void computeStateErrorProbPT19(pll_state_t state,
                              std::vector<double>::iterator &clvp) const;


   void computeStateErrorProbPT17(pll_state_t state,
                                                      std::vector<double>::iterator &clvp) const;
    void computeStateErrorProbPT20(pll_state_t state,
                                                       std::vector<double>::iterator &clvp) const;
    
  
};

#endif /* genotype_error_model_hpp */

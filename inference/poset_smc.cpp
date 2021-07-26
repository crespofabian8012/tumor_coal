//
//  poset_smc.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#include "poset_smc.hpp"

#include <iostream>

std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    
    std::shared_ptr<State> result = std::make_shared<State>(0.0);
                                                            
    result->initForest( params.sampleSize,params.msa, params.positions, params.getProgramOptions());
     
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    std::shared_ptr<State> result;// = std::make_shared<State>(<#_Args &&__args...#>) ;
    //
    //choose 2 random trees to merge
    
    return result;
}
   //double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<double> > &genealogy, const SVModelParams &params);
void PosetSMC::generate_data(gsl_rng *random, size_t T, SMCOptions &params, std::vector<double> &latent, std::vector<double> &obs){
    
    
    
    
}
PosetSMC::~PosetSMC()
{
    
}

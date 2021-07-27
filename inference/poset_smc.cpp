//
//  poset_smc.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#include "poset_smc.hpp"
#include "state.hpp"

#include "poset_smc_params.hpp"


#include <iostream>

std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    
    std::shared_ptr<State> result = std::make_shared<State>(params);
                                                            
    result->initForest( params.sampleSize,params.msa, params.positions, params.getProgramOptions());
     
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    std::shared_ptr<State> result;
    //
    //choose 2 random trees to merge
    
    return result;
}
   //double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<double> > &genealogy, const SVModelParams &params);
void PosetSMC::generate_data(gsl_rng *random, size_t T, SMCOptions &params, std::vector<double> &latent, std::vector<double> &obs){
    
    
    
    
}

double PosetSMC::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<State> > &genealogy, const PosetSMCParams &p){
    
    double result=0.0;
    return result;
    
}
   
   
void PosetSMC::set_particle_population(const std::vector<shared_ptr<State> > &particles) {
    
    
    
    
    
}
PosetSMC::~PosetSMC()
{
    
}

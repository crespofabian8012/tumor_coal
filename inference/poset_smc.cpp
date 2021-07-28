//
//  poset_smc.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#include "poset_smc.hpp"
#include "state.hpp"
#include "random.h"
#include "poset_smc_params.hpp"


#include <iostream>

std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    
    std::shared_ptr<State> result = std::make_shared<State>(params);
                                                            
    result->initForest( params.sampleSize,params.msa, params.positions, params.getProgramOptions());
     
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    std::shared_ptr<State> result;
    
    if (curr.root_count()>1){
        // select 2 nodes to coalesce
        Population *pop;
        long double minProposedTime = 100000000 ;
        int idx_pop;
        for(size_t i=0; i <  curr.getPopulationSet()->getPopulations().size(); i++){
            
            pop= curr.getPopulationSet()->getPopulationbyIndex(i);
            long double time= pop->proposeTimeNextCoalEvent(random, pop->numActiveGametes, params.getProgramOptions().K);
            
            if (time<minProposedTime){
                minProposedTime=time;
                idx_pop=i;
                
            }
        }
        //choose 2 active lineages to coalesce inside population i
        
        
        
    }
    
  
    
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

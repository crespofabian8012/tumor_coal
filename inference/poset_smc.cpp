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

PosetSMC::PosetSMC(size_t numClones, size_t num_iter)
  {
        this->numClones=numClones;
      this->num_iter = num_iter;
         std::cout<<"creating posetSMC"<< std::endl;
}


std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    
    std::shared_ptr<State> result = std::make_shared<State>(params, random );
     
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    
    //make a  copy of curr
    std::shared_ptr<State> result = std::make_shared<State>(curr);
    
    if (curr.root_count()>1){
        // select 2 nodes to coalesce
        Population *pop;
        long double minProposedTime = 100000000 ;
        int idx_pop;
        for(size_t i=0; i <  curr.getNumberPopulations(); i++){
            
            pop= curr.getPopulationByIndex( i);
            
        
            
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
//void PosetSMC::generate_data(gsl_rng *random, size_t T, SMCOptions &params, std::vector<double> &latent, std::vector<double> &obs){
//    

//}
unsigned long PosetSMC::num_iterations(){
    
    return  num_iter;
    
}
double PosetSMC::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<State> > &genealogy, const PosetSMCParams &p){
    
    double result=0.0;
     int a=3;
      a++;
    return result;
    
}
   
   
void PosetSMC::set_particle_population(const std::vector<shared_ptr<State> > &particles) {
    
  
      int a=3;
      a++;
    
    
}


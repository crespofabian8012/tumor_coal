//
//  poset_smc.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#include "utils.hpp"
#include "poset_smc.hpp"
#include "state.hpp"
#include "random.h"
#include "poset_smc_params.hpp"


#include <iostream>

PosetSMC::PosetSMC(size_t numClones, size_t num_iter)
{
    this->numClones=numClones;
    this->num_iter = num_iter;
    
}
std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    log_w = 0;
    std::shared_ptr<State> result = std::make_shared<State>(params, random );
    log_w = result->getInitialLogWeight();
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    
    //make a  copy of curr
    std::shared_ptr<State> result(make_shared<State>(curr));
    log_w=0;
    double logLikNewHeight = 0.0;
    double logWeight = 0.0;
    double newHeight = 0.0;
    if (result->root_count()>1){
        // select 2 nodes to coalesce
        Population *leftNodePop, *rightNodePop;
        int  idxleftNodePop = 0;
        int  idxRightNodePop = 0;
        double timeNextCoalEvent = 0.0;
        timeNextCoalEvent = result->getNextCoalTime(random, idxleftNodePop, idxRightNodePop, logLikNewHeight, params.getProgramOptions().K);
        
        assert(timeNextCoalEvent > 0);
        newHeight = timeNextCoalEvent* result->getTheta();
        std::cout << " new height " << newHeight << std::endl;
        std::cout << " loglik new height " << logLikNewHeight << std::endl;
    
        result->insertNewCoalTime( newHeight);
       // result->setHeightScaledByTheta(timeNextCoalEvent* result->getTheta());
         
        leftNodePop = result->getPopulationByIndex(idxleftNodePop);
        rightNodePop = result->getPopulationByIndex(idxRightNodePop);
       
        if (params.doPriorPost){
            
               logWeight = result->proposalPriorPost( random, leftNodePop,rightNodePop, newHeight, logLikNewHeight);
            }
            else{//PriorPrior
                logWeight =  result->proposalPriorPrior(random, leftNodePop,rightNodePop, newHeight, logLikNewHeight);
            }
            log_w += logWeight ;
        std::cout << " loglik new particle " << log_w << std::endl;
     }
       
    return result;
}
unsigned long PosetSMC::num_iterations(){
    
    return  num_iter;
    
}
double PosetSMC::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<State> > &genealogy, const PosetSMCParams &p){
    
    double result=0.0;
    return result;
    
}



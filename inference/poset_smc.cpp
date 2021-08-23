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
    
}
std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    log_w = 0;
    std::shared_ptr<State> result = std::make_shared<State>(params, random );
    log_w = result->getInitialLogWeight();
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    
    //make a  copy of curr
    // to avoid 2 allocations like in the next line
    // std::shared_ptr<State> result(new State(curr));
    //std::cout << " particle with delta : "<< curr.getPopulationByIndex(0)->delta <<std:: endl;
    std::shared_ptr<State> result(make_shared<State>(curr));
    log_w=0;

    if (result->root_count()>1){
        // select 2 nodes to coalesce
        Population *pop,*chosenPop, *oldestPop;
        long double minTimeNextEvent = DOUBLE_INF ;
        
        long double timeNextEvent;

        for(size_t i=0; i <  result->getNumberPopulations(); i++){
            
            pop= result->getPopulationByIndex( i);
            
            oldestPop = result->getPopulationByIndex( result->getNumberPopulations()-1);
            //std::cout << "num active gametes pop " << pop->numActiveGametes<< std::endl;
            
            timeNextEvent = params.getPopulationEvent(i, result->getIdNextCoalEventForPopulation(i));
            
            if (timeNextEvent<minTimeNextEvent){
                minTimeNextEvent=timeNextEvent;
                chosenPop=pop;
            }
        }
        assert(timeNextEvent > 0);
        
        result->moveNextIdEventForPopulation(chosenPop->index);
        if (chosenPop->numActiveGametes >1){
            //next event a coalescence
            //choose 2 active lineages to coalesce inside population chosenPop
            int idxFirst = 0,  idxSecond = 0;
            int firstInd = 0,  secondInd = 0;
            int choosePairIndividuals = YES;
            
            //std::cout << "num active gametes chosenPop " << chosenPop->numActiveGametes<< std::endl;
        
            if (result->getNumberPopulations()==1)
                assert(chosenPop->numActiveGametes == result->root_count());
            
            chosenPop->ChooseRandomIndividual(&idxFirst, numClones,   &idxSecond, random, choosePairIndividuals);
            
            firstInd = chosenPop->idsActiveGametes[idxFirst];
            secondInd = chosenPop->idsActiveGametes[idxSecond];
            
            assert(firstInd != secondInd);
        
            
            unsigned int chosep_pop_idx = chosenPop->index;
            
            assert(minTimeNextEvent* result->getTheta()>0);
            result->setHeightScaledByTheta(minTimeNextEvent* result->getTheta());
    
            
            std::shared_ptr<PartialTreeNode> node = result->connect(firstInd, secondInd, chosep_pop_idx);
        
            
            result->updateIndexesActiveGametes( idxFirst, idxSecond,  chosep_pop_idx, chosep_pop_idx,  node->index, node->index_population);
       
            
            chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
         
            double weight = result->likelihood_factor(node);
            
            // log_w += weight;
            assert(!isnan(weight) && !isinf(weight));
            
            for (size_t i =0; i < result->root_count(); i ++)
                log_w += result->getRootAt(i)->ln_likelihood;
        
            
        }
        else{//next event is an origin
            
            result->setHeightScaledByTheta(minTimeNextEvent* result->getTheta());
            
        }
        
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



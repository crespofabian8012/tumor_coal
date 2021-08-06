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
    State * newState= new State(curr);
    std::shared_ptr<State> result(newState);
    log_w=0;
    
    if (result->root_count()>1){
        // select 2 nodes to coalesce
        Population *pop,*chosenPop, *oldestPop;
        long double minTimeNextEvent = 100000000 ;
        int idx_pop=1;
        long double waitingTime, timeNextEvent;
        bool isCoalescentEvent=false;
        for(size_t i=0; i <  result->getNumberPopulations(); i++){
            
            pop= result->getPopulationByIndex( i);
            oldestPop = result->getPopulationByIndex( result->getNumberPopulations()-1);
        
            if (pop->numActiveGametes >1){
                waitingTime= pop->proposeTimeNextCoalEvent(random, pop->numActiveGametes, params.getProgramOptions().K);
                
                timeNextEvent = result->getHeight()+ waitingTime;
                isCoalescentEvent=true;
                //do I need to convert this proposed time?
            }
            else{//already reached the MRCA
                
                if (result->getNumberPopulations()==1 || pop == oldestPop){
                    timeNextEvent = pop->timeOriginSTD;
                }
                else{
                    
                    timeNextEvent= pop->timeOriginSTD * pop->x / oldestPop->x;
                }
                isCoalescentEvent=false;
            }
            if (timeNextEvent<minTimeNextEvent){
                minTimeNextEvent=timeNextEvent;
                idx_pop=i;
            }
        }
        
        chosenPop= result->getPopulationByIndex( idx_pop);
        if (chosenPop->numActiveGametes >1){
            //next event a coalescence
            //choose 2 active lineages to coalesce inside population idx_pop
            
            int firstInd,  secondInd=0;
            int choosePairIndividuals = YES;
               
            chosenPop->ChooseRandomIndividual(&firstInd, numClones,   &secondInd, random, choosePairIndividuals);
            
            std::cout<< "first idx "<< firstInd << " second " << secondInd<< std::endl;
            std::cout<< "timeNextEvent "<<  minTimeNextEvent<<  std::endl;
            
            //TODO: the indexes firstInd and secondInd are respect to the chosePop
            //for more than one clon we have to fix this
            
            unsigned int index_population = idx_pop;
            std::shared_ptr<PartialTreeNode> node = result->connect(firstInd, secondInd, minTimeNextEvent-result->getHeight(), index_population);
            
            
            
            long double logLikNextCoal = chosenPop->logLikelihoodNextCoalescent(minTimeNextEvent, result->getHeight(), chosenPop->numActiveGametes, params.getProgramOptions().K);
            
            node->ln_likelihood = node->ln_likelihood +logLikNextCoal;
            
    
            long double logLikNoCoal;
            //do we  need to add loglik of no coal events in  other populations?
          
            for(size_t i=0; i <  result->getNumberPopulations(); i++){
                    pop= result->getPopulationByIndex( i);
                // do we need to include this in the node->ln_likelihood?
                    log_w +=pop->LogDensityTimeSTDFrom(pop->timeOriginSTD, 0);
                    
                    if (pop!=chosenPop){
                        logLikNoCoal=pop->LogProbNoCoalescentEventBetweenTimes(result->getHeight(), minTimeNextEvent,  pop->numActiveGametes, pop->timeOriginSTD, pop->delta, params.getProgramOptions().K);
                       
                        node->ln_likelihood = node->ln_likelihood +logLikNoCoal;
                        
                    }
                }
            
            
            double weight = result->likelihood_factor(node);
            
            assert(!isnan(weight) && !isinf(weight));
                   
            log_w = weight;
            
        }
        else{//next event an origin
            
            result->setHeight(minTimeNextEvent);
            //the new log_weight is the log of probability of no events
            // in any of the other populations
//             long double logLikNoCoal=pop->LogProbNoCoalescentEventBetweenTimes(curr.getHeight(), minTimeNextEvent,  pop->numActiveGametes, pop->timeOriginSTD, pop->delta, params.getProgramOptions().K);
            
        }
       
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
   
    return result;
    
}
   
   
void PosetSMC::set_particle_population(const std::vector<shared_ptr<State> > &particles) {
    
    
    
}


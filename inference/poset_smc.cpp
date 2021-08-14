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
        long double minTimeNextEvent = DOUBLE_INF ;
        int idx_pop=1;
        long double waitingTime, timeNextEvent;
        bool isCoalescentEvent=false;
        double currentTime=0.0;
        double currentTimeKingman=0.0;
        for(size_t i=0; i <  result->getNumberPopulations(); i++){
            
            pop= result->getPopulationByIndex( i);
            oldestPop = result->getPopulationByIndex( result->getNumberPopulations()-1);
        
            if (pop->numActiveGametes >1){
                
                currentTime = result->getHeightModelTime() ;
                currentTimeKingman = Population::FmodelTstandard (currentTime , pop->timeOriginSTD, pop->delta,   params.getProgramOptions().K);//this is current time in Kingman coal time
                
                waitingTime= pop->proposeTimeNextCoalEvent(random, pop->numActiveGametes, params.getProgramOptions().K);
                
                currentTimeKingman = currentTimeKingman+waitingTime;
                
                timeNextEvent =   Population::GstandardTmodel(currentTimeKingman, pop->timeOriginSTD, pop->delta, params.getProgramOptions().K);;
                isCoalescentEvent=true;
                
            }
            else{//last event
                
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
        // std::cout<< "chosen pop "<< idx_pop <<  " T " << chosenPop->timeOriginSTD << std::endl;
        if (chosenPop->numActiveGametes >1){
            //next event a coalescence
            //choose 2 active lineages to coalesce inside population idx_pop
            
            int firstInd,  secondInd=0;
            int choosePairIndividuals = YES;
               
            chosenPop->ChooseRandomIndividual(&firstInd, numClones,   &secondInd, random, choosePairIndividuals);
            
            //std::cout<< "first idx "<< firstInd <<  " second " << secondInd<< std::endl;
           // std::cout<< "first label "<< result->getRootAt(firstInd)->label <<  " second " << result->getRootAt(secondInd)->label<< std::endl;
           // std::cout<< "timeNextEvent "<<  minTimeNextEvent<<  " T "<< chosenPop->timeOriginSTD  << std::endl;
            
            //TODO: the indexes firstInd and secondInd are respect to the chosePop
            //for more than one clon we have to fix this
            
            unsigned int index_population = idx_pop;
            result->setHeightScaledByTheta(minTimeNextEvent* result->getTheta());
            
           // std::cout << "new  height scaled by theta"<<  " is "<< result->getHeightScaledByTheta()<< std::endl;
            
            std::shared_ptr<PartialTreeNode> node = result->connect(firstInd, secondInd, minTimeNextEvent* result->getTheta()-result->getHeightScaledByTheta(), index_population);
            
            
            long double logLikNextCoal = 0.0;
            
       
            
            logLikNextCoal = chosenPop->logLikelihoodNextCoalescent(minTimeNextEvent, result->getHeightModelTime(), chosenPop->numActiveGametes, params.getProgramOptions().K);
            
            node->ln_likelihood = node->ln_likelihood +logLikNextCoal;
            
    
            long double logLikNoCoal;
            //do we  need to add loglik of no coal events in  other populations?
          
            for(size_t i=0; i <  result->getNumberPopulations(); i++){
                    pop= result->getPopulationByIndex( i);
                // do we need to include this in the node->ln_likelihood?
                    log_w +=pop->LogDensityTimeSTDFrom(pop->timeOriginSTD, 0);
                    
                    if (pop!=chosenPop){
                        logLikNoCoal=pop->LogProbNoCoalescentEventBetweenTimes(result->getHeightModelTime(), minTimeNextEvent,  pop->numActiveGametes, pop->timeOriginSTD, pop->delta, params.getProgramOptions().K);
                       
                        node->ln_likelihood = node->ln_likelihood +logLikNoCoal;
                        
                    }
                }
            
         //    std::cout << "log lik of parent node after the coal part "<<  " is "<< node->ln_likelihood<< std::endl;
            double weight = result->likelihood_factor(node);
            
             
           assert(!isnan(weight) && !isinf(weight));
                   
            log_w = weight;
            
            //std::cout << "log weight "<<  " is "<<  log_w<< std::endl;
            
        }
        else{//next event is an origin
            
            result->setHeightScaledByTheta(minTimeNextEvent* result->getTheta());
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


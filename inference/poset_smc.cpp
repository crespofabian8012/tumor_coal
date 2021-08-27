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
    // to avoid 2 allocations like in the next line
    // std::shared_ptr<State> result(new State(curr));
    // I used the following
    std::shared_ptr<State> result(make_shared<State>(curr));
    log_w=0;
    double weight = 0.0;
    
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
            int idxFirstRoot = 0,  idxSecondRoot = 0;
            int choosePairIndividuals = YES;
            unsigned int chosep_pop_idx = chosenPop->index;
            //std::cout << "num active gametes chosenPop " << chosenPop->numActiveGametes<< std::endl;
            
            assert(minTimeNextEvent* result->getTheta()>0);
            result->setHeightScaledByTheta(minTimeNextEvent* result->getTheta());
                           
            
            if (params.doPriorPost){
                
                std::vector<std::pair<int, int>> allCoalPairs = Utils::allCombinations(chosenPop->numActiveGametes, 2);
                int numberProposals = allCoalPairs.size();
                if (t>1)
                    numberProposals = 2* numberProposals;
              
                std::vector<double> logWeights(numberProposals, 0.0);
                std::vector<double> normWeights(numberProposals, 0.0);
                
                std::vector<std::shared_ptr<PartialTreeNode>> nodeProposals(numberProposals);
                
                for(size_t i=0; i <  allCoalPairs.size(); i++){
                    idxFirst = allCoalPairs[i]. first;
                    idxSecond = allCoalPairs[i]. second;
                    
                    idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
                    idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
                    
                   // std::cout << "idx first" << idxFirst << " idx second "<<idxSecond << std::endl;
                    assert(idxFirstRoot != idxSecondRoot);
                    if (t==1){
                       nodeProposals[i] = result->proposeNewNode( idxFirstRoot,  idxSecondRoot, chosep_pop_idx );
                        logWeights[i] = nodeProposals[i]->ln_likelihood;
                    }
                    else{
                        nodeProposals[2*i] = result->proposeNewNode( idxFirstRoot,  idxSecondRoot, chosep_pop_idx );
                        logWeights[2*i] = nodeProposals[2*i]->ln_likelihood;
                        nodeProposals[2*i+1] = result->proposeNewNode(idxSecondRoot, idxFirstRoot,   chosep_pop_idx );
                        logWeights[2*i+1] = nodeProposals[2*i+1] ->ln_likelihood;
                    }
         
                }
                unsigned int pos;
                weight = Utils::normalize(logWeights, normWeights);
                pos = Random::randomDiscreteFromProbabilityVector(random, &normWeights[0], normWeights.size());
              

                idxFirst = allCoalPairs[pos]. first;//first idx inside the list of active gametes of the chosen pop
                idxSecond = allCoalPairs[pos]. second;//second idx inside the list of active gametes of the chosen pop
                
                std::cout<< "first " << idxFirst << "second " << idxSecond << " weight " << weight << " norm weight " << normWeights[pos] << std::endl;
                
                idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
                idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
                
                assert(idxFirstRoot != idxSecondRoot);
                
                int posFirstRoot= result->getNodeIdxById(idxFirstRoot);
                int posSecondRoot = result->getNodeIdxById(idxSecondRoot);
                
           
                result->remove_roots(posFirstRoot, posSecondRoot);
                result->addRoot(nodeProposals[pos]);
                
              
                result->increaseNextAvailable();
                result->updateIndexesActiveGametes( idxFirst, idxSecond,  chosep_pop_idx, chosep_pop_idx,  nodeProposals[pos]->index, nodeProposals[pos]->index_population);
                               
                chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
            }
            else{//PriorPrior
                
               
                if (result->getNumberPopulations()==1)
                    assert(chosenPop->numActiveGametes == result->root_count());
                
        
                chosenPop->ChooseRandomIndividual(&idxFirst, numClones,   &idxSecond, random, choosePairIndividuals);
                
                idxFirstRoot = chosenPop->idsActiveGametes[idxFirst];
                idxSecondRoot = chosenPop->idsActiveGametes[idxSecond];
                
                assert(idxFirstRoot != idxSecondRoot);
                
                std::shared_ptr<PartialTreeNode> node = result->connect(idxFirstRoot, idxSecondRoot, chosep_pop_idx);
                
                result->updateIndexesActiveGametes( idxFirst, idxSecond,  chosep_pop_idx, chosep_pop_idx,  node->index, node->index_population);
                
                chosenPop->numActiveGametes= chosenPop->numActiveGametes-1;
                
                weight = result->likelihood_factor(node);
                assert(!isnan(weight) && !isinf(weight));
                
            }
            
            log_w += weight;
        
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



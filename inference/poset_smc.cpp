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
#include "gnu_plotter.hpp"
#include <Eigen/Core>
#include <cmath>
#include <boost/tuple/tuple.hpp>




#include <iostream>
using ListListDouble = std::vector<std::vector<double>>;
using PairListListDouble = std::pair<ListListDouble,ListListDouble>;
std::unordered_map<size_t , std::set<pairs > > PosetSMC::sizeCombinationMap = std::unordered_map<size_t , std::set<pairs > >() ;
PosetSMC::PosetSMC(size_t numClones, size_t num_iter, bool doPlots)
{
    this->numClones=numClones;
    this->numIter = num_iter;
    for( size_t i = 0; i< num_iter; ++i ){
        
        if (doPlots)
         doPlotPerIteration.push_back(true);
        else
            doPlotPerIteration.push_back(false);
    }
    
}
std::shared_ptr<State> PosetSMC::propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params){
    log_w = 0;
    std::shared_ptr<State> result = std::make_shared<State>(params, random );
    log_w = result->getLogWeight();
    return result;
}
std::shared_ptr<State> PosetSMC::propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params){
    
    //make a  copy of curr
    std::shared_ptr<State> result(make_shared<State>(curr));
    log_w=0;
    double logLikNewHeight = 0.0;
    double logWeightDiff = 0.0;
    double currentLogWeight = result->getLogWeight();
   // std::cout<< "current log weight " << currentLogWeight << std::endl;
    double newHeight = 0.0;
    double K = params.getProgramOptions().K;
    std::set<std::pair<int, int>> setPairs;
 
    if (result->root_count()>1){
        // select 2 nodes to coalesce
        Population *receiverPop, *incomingPop;
        int  idxReceiverPop = 0;
        int  idxIncomingPop = 0;
        double timeNextCoalEvent = 0.0;
        if (doPlotPerIteration[t]){
            Population *currPop = result->getCurrentPopulation();
                 //double from = currPop->currentModelTime * currPop->x * result->getTheta()+0.001;
                 double from = 0.0001;
                 double to = 0.98* currPop->timeOriginSTD *currPop->x * result->getTheta() ;
                 int n= 29;
                 Eigen::ArrayXf xseries = Eigen::ArrayXf::LinSpaced(n, from , to);
                 std::vector<std::string> labels;
                 GNUPlotter plotter;
                 //std::vector<double> logLikGridPerIncrement = result->evalLogLikRatioGridPerIncrement( from,  to,  n, labels);
                 //plotter.plot2dSerie(from, to,  n, logLikGridPerIncrement, labels );
                 PairListListDouble pairIncrementsLogLiks = result->evalLogLikRatioGrid(from, to, n, labels);
                 std::string fileName = to_string(currPop->sampleSize)+"_iteration_"+to_string(t)+".png";
                 
                 //std::filesystem::path cwd = std::filesystem::current_path() / fileName;
                 //std::string path = ""
                 std::string filePath = to_string(currPop->sampleSize)+"_iteration_"+to_string(t)+".png";
                 double trueCurrenCoalTime = params.trueCoalTimes[t-1];
                 double trueNextCurrenCoalTime = params.trueCoalTimes[t];
                 plotter.plot2dSeveralSeries(pairIncrementsLogLiks.first, pairIncrementsLogLiks.second, labels,  true, filePath, trueCurrenCoalTime, trueNextCurrenCoalTime );
                doPlotPerIteration[t] = false;
        }
        
        if(kernelType == PRIORPRIOR) {//PriorPrior: samples increment from the prior and pair from the prior
            timeNextCoalEvent = result->getNextCoalTime(random, idxReceiverPop, idxIncomingPop, logLikNewHeight, K);
            
            assert(timeNextCoalEvent > 0);
            newHeight = timeNextCoalEvent* result->getTheta();
            result->insertNewCoalTime( newHeight);
            receiverPop = result->getPopulationByIndex(idxReceiverPop);
            incomingPop = result->getPopulationByIndex(idxIncomingPop);
            logWeightDiff =  result->proposalPriorPrior(random, receiverPop, incomingPop, newHeight, logLikNewHeight);
            result->acceptNextCoalTimeProposal(timeNextCoalEvent, random, idxReceiverPop, idxIncomingPop, K);
        }
        else{
            
            
            if (kernelType == PRIORPOST){//PriorPrior: samples increment from the prior and pair from the posterior
                
                
                timeNextCoalEvent = result->getNextCoalTime(random,  idxReceiverPop, idxIncomingPop, logLikNewHeight, K);
                assert(timeNextCoalEvent > 0);
                
                
                newHeight = timeNextCoalEvent* result->getTheta();
               // std::cout<< "new height " << newHeight << " logLikNewHeight "<< logLikNewHeight<< std::endl;
                result->insertNewCoalTime( newHeight);
                receiverPop = result->getPopulationByIndex(idxReceiverPop);
                incomingPop = result->getPopulationByIndex(idxIncomingPop);
                logWeightDiff = result->proposalPriorPost( random, receiverPop,incomingPop, newHeight, logLikNewHeight, t);
                
                result->acceptNextCoalTimeProposal(timeNextCoalEvent, random, idxReceiverPop, idxIncomingPop,K);
            }
            else if (kernelType == POSTPOST1){//POSTPOST: samples increment from the posterior and then samples pair from the posterior conditional on  increment
                          
                logWeightDiff =  result->proposalPostPost1(random, timeNextCoalEvent, logLikNewHeight, numIncrementsPOSTPOST, K);
        
                      }
                      
            else if (kernelType == POSTPOST2){//POSTPOST: samples increment  and pair from the posterior
                
                
                logWeightDiff =  result->proposalPostPost(random, timeNextCoalEvent, logLikNewHeight, numIncrementsPOSTPOST, K);
                
            }
            
        }
        
        assert(timeNextCoalEvent > 0);
        newHeight = timeNextCoalEvent* result->getTheta();
        if (params.verbose>1){
            std::cout << " new height " << newHeight << std::endl;
            std::cout << " loglik new height " << logLikNewHeight << std::endl;
        }
        unsigned int numNonTrivialTrees = result->numberNonTrivialTrees();
        assert(numNonTrivialTrees>0);
        //log_w += logWeight -logLikNewHeight-log(result->getTheta())- log(numNonTrivialTrees) ;
        if (result->getNumberPopulations() > 1.0)//nonClockTrees
            logWeightDiff += logWeightDiff- log(numNonTrivialTrees) ;
        else
            logWeightDiff += logWeightDiff;
        if (params.verbose>1)
            std::cout << " loglik new particle " << log_w << std::endl;
    }
    //std::cout<< "new log weight " << currentLogWeight +logWeightDiff << std::endl;
    result->setLogWeight(currentLogWeight +logWeightDiff );
    log_w = currentLogWeight +logWeightDiff  ;
    return result;
}
unsigned long PosetSMC::num_iterations(){
    
    return  numIter;
    
}
double PosetSMC::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<State> > &genealogy, const PosetSMCParams &p){
    
    double result=0.0;
    return result;
    
}


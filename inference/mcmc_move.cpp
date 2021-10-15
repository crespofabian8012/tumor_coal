/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/

/*
 * mcmc move  class
 */
#include "mcmc_move.hpp"
#include "mcmc_chain.hpp"
#include "data_utils.hpp"
#include "utils.hpp"

#include <algorithm>

#include "random.h"
//#include <algorithm>

//using namespace std;

MCMCmove::MCMCmove(int chainId,  std::string nameMove, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences)
{
    this->chainId = chainId;
    this->nameMove = nameMove;
    this->mcmcParamKernels= mcmcParKernel;
    this->affectedParameters = affectedParameters;
    this->vectorParam = &vectorParam;
    this->numberAccept=0;
    this->numberReject=0;
    this->numberAttemps=0;
    this->isInvalidMove=false;
    this->newLogConditionalLikelihoodTree = currentlogConditionalLikelihoodTree;
    this->newLogConditionalLikelihoodSequences = currentlogConditionalLikelihoodSequences;
    this->chainLogConditionalLikelihoodTree = currentlogConditionalLikelihoodTree;
    this->chainLogConditionalLikelihoodSequences= currentlogConditionalLikelihoodSequences;
}
double MCMCmove::numberAccepted()
{
    return this->numberAccept;
}
std::string MCMCmove::name() {
    return this->nameMove;
}
double MCMCmove::numberRejected()
{
    return this->numberReject;
}
void  MCMCmove::resetNumberAccepted(){
   this->numberAccept=0;
}
void  MCMCmove::resetNumberRejected(){
    this->numberReject=0;
}
int  MCMCmove::getChainId()
{
    return chainId;
}
void  MCMCmove::setCurrentLogLikelihoods(long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences)
{
    this->chainLogConditionalLikelihoodTree= currentlogConditionalLikelihoodTree;
    this->chainLogConditionalLikelihoodSequences= currentlogConditionalLikelihoodSequences;
}
long double  MCMCmove::getCurrentLogLikelihoodTree() const
{
    return this->chainLogConditionalLikelihoodTree;
}
long double  MCMCmove::getCurrentLogLikelihoodSequences() const
{
    return this->chainLogConditionalLikelihoodSequences;
}

NewProportionsVectorMove::NewProportionsVectorMove(int chainId, std::string nameMove,std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences ):MCMCmove(chainId, nameMove, mcmcParKernel, affectedParameters, vectorParam,
                                                                currentlogConditionalLikelihoodTree,
                                                                 currentlogConditionalLikelihoodSequences)
{
    
}
void MCMCmove::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const   gsl_rng *randomGsl, boost::mt19937* rngBoost, Chain* chain)
{
    saveCurrentValue(programOptions,mcmcOptions, chain );
    numberAttemps=0;
    makeProposal(programOptions, mcmcOptions, randomGsl, chain);
    
    if (isInvalidMove){
        rollbackMove(programOptions,mcmcOptions, chain);
        this->numberReject++;
        isInvalidMove=false;
        if (mcmcOptions.verbose>=2){
            std::cout << "\n Rejected new move because invalid " << nameMove.c_str() << std::endl;
        }
        return;
    }
   long  double logAcceptanceRate =computeLogAcceptanceProb(programOptions, mcmcOptions, chain);
  
    if (mcmcOptions.verbose>=1)
      std::cout << "\n log Acceptance rate " << logAcceptanceRate <<std::endl;
    
    long double randomNumber;
    
    if(mcmcOptions.useGSLRandomGenerator)
        randomNumber= Random::randomUniformFromGsl2(randomGsl);
    else
        randomNumber= Random::randomUniformBoost(  rngBoost);
    
    long double logRandom=log(randomNumber);
    if (mcmcOptions.verbose>=1)
      std::cout << "\n logRandom " << logRandom <<std::endl;
    if (logRandom < logAcceptanceRate )
    {//accept the move
        if (mcmcOptions.verbose>=2)
           std::cout<< "\n Accepted new move" << nameMove.c_str() << std::endl;
        this->numberAccept++;
        
        //Chain *chain=getChain();
        if ( (!mcmcOptions.noData && chain->currentlogConditionalLikelihoodTree !=newLogConditionalLikelihoodTree))
        {
            chainLogConditionalLikelihoodTree = newLogConditionalLikelihoodTree;
            chain->currentlogConditionalLikelihoodTree =chainLogConditionalLikelihoodTree;
        }
        if ((!mcmcOptions.noData && mcmcOptions.useSequencesLikelihood ==1 && chainLogConditionalLikelihoodSequences!=newLogConditionalLikelihoodSequences)){
            chainLogConditionalLikelihoodSequences = newLogConditionalLikelihoodSequences;
            chain->currentlogConditionalLikelihoodSequences =chainLogConditionalLikelihoodSequences;
        }
            
    }
    else
    {
        rollbackMove(programOptions,mcmcOptions, chain);
        this->numberReject++;
        if (mcmcOptions.verbose>=2)
         std::cout<< "\n Rejected new move " << nameMove.c_str()<< std::endl;
    }
}
//////////////////////////////////////////////
void NewProportionsVectorMove::saveCurrentValue(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
    int i=0;
    Population *popI;
    chain->oldproportionsVector= chain->proportionsVector;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->oldx= popI->x;

        popI->oldFatherPop = popI->FatherPop;
        popI->oldDeathRate = popI->deathRate;
        popI->oldGrowthRate = popI->growthRate;
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewProportionsVectorMove::makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,const  gsl_rng *randomGenerator,Chain* chain)
{
   // Chain *chain=getChain();
    int i=0;
    Population *popI;
    //long double *proportionsVectorArray;
    std::vector<long double> effectPopSizesvector(chain->numClones);
    bool allPopulationPopSizesSet=false;
    double totalSampleSize= chain->totalSampleSize();
    for( i = 0 ; i < chain->numClones; i++)
    {
        effectPopSizesvector.at(i) = chain->oldproportionsVector.at(i) * totalSampleSize;
        //fprintf (stderr, "\n old effect population size  at %d: %lf \n",i,effectPopSizesvector.at(i) );
    }
    do {
        allPopulationPopSizesSet=true;
        numberAttemps++;
        Random::randomDirichletFromVector (effectPopSizesvector, chain->proportionsVector, true, randomGenerator, NULL);
        //        randomDirichletFromGsl(chain->numClones, proportionsVectorArray, &(chain->proportionsVector[0]));
        chain->updateEffectPopSizesCurrentProportionsVector();
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=chain->populations[i];
            popI->x= chain->proportionsVector.at(i);
            
            popI->growthRate =popI->delta  / popI->x;
            
            popI->deltaT =popI->growthRate * totalSampleSize;
            popI->theta = chain->theta * popI->x;
            //popI->popSize=popI->x * popI->birthRate;
            popI->deathRate= popI->birthRate - popI->growthRate;
            
            popI->numCompletedCoalescences=0;
            popI->CoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime.clear();
            popI->numIncomingMigrations=0;

        }
        if (!allPopulationPopSizesSet)
            break;
        chain->InitListPossibleMigrations();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        allPopulationPopSizesSet = chain->checkMigrationsOrder();
        if (numberAttemps >mcmcOptions.maxNumberProposalAttempts)
            break;
    }
    while(!allPopulationPopSizesSet);
    if (allPopulationPopSizesSet)
    {
        //        for( i = 0 ; i < chain->numClones; i++)
        //         fprintf (stderr, "\n new proportions vector at %d: %lf \n",i,chain->proportionsVector.at(i) );
    }
}
void NewProportionsVectorMove::rollbackMove(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    chain->proportionsVector= chain->oldproportionsVector ;
    int i=0;
    Population *popI;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->x= popI->oldx;
        popI->FatherPop = popI->oldFatherPop;
        popI->deathRate = popI->oldDeathRate;
        popI->growthRate = popI->oldGrowthRate;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
        popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
    }
}
long double NewProportionsVectorMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
    // ListClonesAccordingTimeToOrigin(chain->populations);
   
    std::vector<long double> oldConcentrationParameter(chain->numClones);
    std::vector<long double> proposedConcentrationParameter(chain->numClones);
    double totalSampleSize= chain->totalSampleSize();
    for(unsigned i = 0 ; i < chain->numClones; i++){
        oldConcentrationParameter.at(i)  = chain->oldproportionsVector.at(i) * totalSampleSize;
        proposedConcentrationParameter.at(i)  = chain->proportionsVector.at(i) * totalSampleSize;
    }
    
   long double priorDensityNewProportionsVector =  chain->DirichletDensity( chain->proportionsVector, oldConcentrationParameter , chain->numClones);
    
    long double priorDensityCurrentProportionsVector = chain->DirichletDensity(chain->oldproportionsVector, proposedConcentrationParameter, chain->numClones);
    
    long double sumLogNumerators= priorDensityNewProportionsVector;
     long double sumLogDenominators=priorDensityCurrentProportionsVector;
    if(!mcmcOptions.noData)
    {
         newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions,mcmcOptions);
        sumLogNumerators+= newLogConditionalLikelihoodTree ;
    
     sumLogDenominators+=chain->currentlogConditionalLikelihoodTree ;
    }
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>=0.0)?0.0: difference;
    
    return(LogAcceptanceRate);
}
NewGrowthRateMoveForPopulation::NewGrowthRateMoveForPopulation(int chainId, std::string nameMove,Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam,long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences):MCMCmove(chainId, nameMove, mcmcParKernel, affectedParameters, vectorParam, currentlogConditionalLikelihoodTree,
    currentlogConditionalLikelihoodSequences)
{
    this->pop = pop;
}
void NewGrowthRateMoveForPopulation::saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
   pop->olddelta= pop->delta;

    Population *popI;
    
    for( int i = 0 ; i < chain->numClones; i++)
    {
        popI= chain->populations[i];
        popI->oldFatherPop = popI->FatherPop;
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewGrowthRateMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator, Chain* chain)
{
    //Chain *chain=getChain();
    
    double randomDelta = Random::RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, mcmcOptions.useGSLRandomGenerator, randomGenerator, NULL);
  
    pop->delta =  randomDelta;
    pop->growthRate =pop->delta  / pop->x;
    pop->deathRate= pop->birthRate - pop->growthRate;
    
    chain->InitListPossibleMigrations();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
void NewGrowthRateMoveForPopulation::rollbackMove(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
    Population *popI;
    pop->delta= pop->olddelta;
    for(int i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->FatherPop = popI->oldFatherPop;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
        popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
    }
}
long double NewGrowthRateMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
    
    //fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total growth rate for population %d,  of the chain %d is = %lf  \n", pop->index, chain->chainNumber, newLogConditionalLikelihoodTree );
    
    long double priorDensityNewScaledGrowthRate =  Distributions::LogUniformDensity(pop->delta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    long double priorDensitScaledGrowthRate = Distributions::LogUniformDensity(pop->olddelta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    long double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewScaledGrowthRate;
    
    long double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensitScaledGrowthRate;
    
    
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>=0.0)?0.0: difference;
    
    return LogAcceptanceRate;
}


void NewProportionsVectorMove::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,const  gsl_rng *rngGsl,boost::mt19937* rngBoost,Chain* chain)
{
    MCMCmove::move(programOptions,mcmcOptions, rngGsl, rngBoost, chain);
}
void NewGrowthRateMoveForPopulation::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain)
{
    MCMCmove::move(programOptions,mcmcOptions, rngGsl, rngBoost, chain);
}

void NewTimeOriginOnTreeforPopulationMove::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,const  gsl_rng *rngGsl,boost::mt19937* rngBoost,Chain* chain)
{
    MCMCmove::move(programOptions,mcmcOptions, rngGsl, rngBoost, chain);
}
void NewTimeOriginOnEdgeforPopulationMove::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain)
{
    MCMCmove::move(programOptions,mcmcOptions, rngGsl, rngBoost, chain);
}
void NewGlobalScaledGrowthRateForPopulationMove::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const  gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain)
{
    MCMCmove::move(programOptions,mcmcOptions, rngGsl, rngBoost, chain);
}
NewGlobalScaledGrowthRateForPopulationMove::NewGlobalScaledGrowthRateForPopulationMove(int chainId,std::string nameMove,  Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel, std::vector<IMCMCParameter<long double> *> &affectedParameters,MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences ):MCMCmove( chainId, nameMove, mcmcParKernel, affectedParameters, vectorParam,  currentlogConditionalLikelihoodTree,
 currentlogConditionalLikelihoodSequences)
{
    this->pop =pop;
}
NewTimeOriginOnTreeforPopulationMove::NewTimeOriginOnTreeforPopulationMove(int chainId,std::string nameMove,  Population *pop,std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences):MCMCmove(chainId, nameMove, mcmcParKernel, affectedParameters, vectorParam,  currentlogConditionalLikelihoodTree,
    currentlogConditionalLikelihoodSequences)
{
    this->pop =pop;
}
NewTimeOriginOnEdgeforPopulationMove::NewTimeOriginOnEdgeforPopulationMove(int chainId,std::string nameMove,  Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam,long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences):MCMCmove(chainId, nameMove, mcmcParKernel, affectedParameters, vectorParam,  currentlogConditionalLikelihoodTree,
    currentlogConditionalLikelihoodSequences)
{
    this->pop =pop;
}
void NewTimeOriginOnTreeforPopulationMove::saveCurrentValue(ProgramOptions &programOptions,MCMCOptions &mcmcOptions,Chain* chain)
{
    //Chain *chain=getChain();
    //save information for the population
    bool isOldestPop= chain->isOldestPopulation(pop, chain->rMRCAPopulation);
    
    pop->oldrMRCA =  pop->rMRCA ;//this only changes to population pop
    
    pop->oldTimeOriginInput = pop->timeOriginInput;
    pop->oldScaledTimeOriginInput = pop->scaledtimeOriginInput;
    pop->oldTimeOriginSTD = pop->timeOriginSTD;
    pop->TimeOriginSTD->saveCurrentValue();
   
    //save information for the other populations
   
         Population *popI;
        for(int i = 0 ; i < chain->numClones; i++)
        {
            popI=chain->populations[i];
      
            popI->oldSampleSize = popI->sampleSize;
           // popI->oldeffectPopSize = popI->effectPopSize;
            popI->oldTimeOriginSTD = popI->timeOriginSTD;
            pop->TimeOriginSTD->saveCurrentValue();
            popI->oldFatherPop =  popI->FatherPop;
            
            if (!isOldestPop)
            {
              popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
              popI->CoalescentEventTimes.clear();
              popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
              popI->immigrantsPopOrderedByModelTime.clear();
            }
        
        }
}
void NewTimeOriginOnTreeforPopulationMove::rollbackMove(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
    int i=0;
    Population *popI;
     bool isOldestPop= chain->isOldestPopulation(pop, chain->rMRCAPopulation);
    //save information for the population
   
    pop->rMRCA = pop->oldrMRCA ;
    
    pop->timeOriginInput = pop->oldTimeOriginInput;
   

    pop->scaledtimeOriginInput = pop->oldScaledTimeOriginInput ;
    pop->timeOriginSTD = pop->oldTimeOriginSTD ;
    
    pop->TimeOriginSTD->rollbackValue();
    //save information for the other populations
    if(programOptions.numClones >1 )
    {
        for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->sampleSize = popI->oldSampleSize;
        //  popI->effectPopSize = popI->oldeffectPopSize;
        if (popI != pop){
          popI->timeOriginSTD = popI->oldTimeOriginSTD;
            popI->TimeOriginSTD->rollbackValue();
        }
        
        popI->FatherPop = popI->oldFatherPop ; //the father population can change also for other populations
        if (!isOldestPop){
            
            popI->CoalescentEventTimes = popI->oldCoalescentEventTimes ;
            popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
            popI->oldCoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime = popI->oldimmigrantsPopOrderedByModelTime;
            popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
            popI->oldimmigrantsPopOrderedByModelTime.clear();
        }
      }
    }
}
void NewTimeOriginOnTreeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator,Chain* chain)
{
    //Chain *chain=getChain();
    //first compute the list of 3 adjacent edges
    //build a cummulative list with the 3 adjacent edges
    //then select at random number between 0 and the sum of the  3 adjacent edges
    //recompute the sample sizes
    // recompute the coalescent and migration events
    //bool existsZeroSampleSizePop=false;
    //double alpha[chain->numClones];
    //int numberPoints = chain->numClones -1;
    std::map<pll_rnode_t*, std::vector<Population*>>  newrmrcaOfPopulation;
    std::map<pll_rnode_t*, std::vector<Population*>>::iterator it;
    
    bool isOldestPop = chain->isOldestPopulation(pop, chain->rMRCAPopulation) ;
    
    if (!isOldestPop)
    {
        chain->proposedrMRCAPopulation=  chain->chooseAvailableEdgeOnRootedTreeForPopulation(pop, chain->rMRCAPopulation, programOptions.healthyTipLabel,randomGenerator );
    }
    else
    {
        chain->proposedrMRCAPopulation=chain->chooseNewTimeofOriginOnEdge(pop,mcmcOptions, randomGenerator, NULL);
        if (pop->timeOriginInput < pop->lowerBoundTimeOriginInput)
        {
            // fprintf (stderr, "\n invalid move \n");
            isInvalidMove = true;
            //exit(-1);
            
        }
    }
    bool stillTheOldestPop=  chain->isOldestPopulation(pop, chain->proposedrMRCAPopulation) ;
    
    if (!stillTheOldestPop)//not the oldest population any more
    {
        chain->currentrMRCAPopulation.clear();
        chain->currentrMRCAPopulation.insert(chain->rMRCAPopulation.begin(), chain->rMRCAPopulation.end());
        
        chain->rMRCAPopulation.clear();
        chain->rMRCAPopulation.insert(chain->proposedrMRCAPopulation.begin(), chain->proposedrMRCAPopulation.end());
        
        chain->initPopulationsSampleSizes( chain->proposedrMRCAPopulation, programOptions.healthyTipLabel);
        
        chain->initTimeOriginSTDYoungerPopulations(mcmcOptions);
        chain->InitListPossibleMigrations();//after setting the timeSTD
        
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->proposedrMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        
    }
    else{//still the oldest population
        //save current dictionary MRCA-population
        chain->currentrMRCAPopulation.clear();
        chain->currentrMRCAPopulation.insert(chain->rMRCAPopulation.begin(), chain->rMRCAPopulation.end());
        
        chain->rMRCAPopulation.clear();
        chain->rMRCAPopulation.insert(chain->proposedrMRCAPopulation.begin(), chain->proposedrMRCAPopulation.end());
        for (unsigned int i = 0; i < chain->numClones; ++i)
        {
            auto popI =  chain->populations[i];
            popI->sampleSize= popI->oldSampleSize ;
            if (popI != pop){
                popI->timeOriginSTD = popI->oldTimeOriginSTD;
                popI->TimeOriginSTD->setParameterValue(popI->timeOriginSTD);
                
            }
            
            popI->FatherPop = popI->oldFatherPop ;
            if(!isOldestPop)
            {
                popI->restoreOldCoalescentTimes();
            }
           
            if (popI != pop){//the other younger populations
                popI->restoreOldImmigrationTimes();
            }
            else
            {//the oldest population itself
               
            }
            
        }
    }

}
long double  NewTimeOriginOnTreeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
   
    std::vector<std::pair<double, pll_tree_edge_t *> > currentAvailableEdges;
    std::vector<std::pair<double, pll_tree_edge_t *> > proposedAvailableEdges;
    
    
    
    bool isOldestPop = chain->isOldestPopulation(pop, chain->rMRCAPopulation);
    long double sumLogNumerators=0.0 ;
    
   long double sumLogDenominators=0.0;
    
  
      
    
    if (!isOldestPop)
    {
         fprintf (stderr, "\n not the oldest pop out of %d populations\n", programOptions.numClones);
        long double currentSumAvailBranchLengths = chain->computeAdjacentEdges(currentAvailableEdges, chain->currentrMRCAPopulation, programOptions.healthyTipLabel, pop, pop->oldrMRCA);
        
        long double newSumAvailBranchLengths =chain->computeAdjacentEdges(proposedAvailableEdges, chain->proposedrMRCAPopulation, programOptions.healthyTipLabel, pop, pop->rMRCA);
        
         double currentProportionsArray[chain->numClones];
        unsigned int currentSampleSizesArray[chain->numClones];
        unsigned int proposedSampleSizesArray[chain->numClones];
        for (unsigned int i = 0; i < chain->numClones; ++i)
        {
            auto pop =  chain->populations[i];
            proposedSampleSizesArray[i]= pop->sampleSize ;
            currentSampleSizesArray[i]= pop->oldSampleSize ;
            currentProportionsArray[i]= pop->x;
            
        }
        long double currentlogMultinomialProb= Distributions::logMultinomialProbability(chain->numClones, currentProportionsArray, currentSampleSizesArray);
        
         long double proposedlogMultinomialProb=Distributions::logMultinomialProbability(chain->numClones, currentProportionsArray, proposedSampleSizesArray);
        
        long double numeratorQ=log(pop->oldrMRCA->length) - log(newSumAvailBranchLengths);
        //if (mcmcOptions.priorsType==0)   priorDensityNewTotalEffectivePopulationSize = LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
        long double denominatorQ= log(pop->rMRCA->length) - log(currentSumAvailBranchLengths);
        
        //  if (mcmcOptions.priorsType==0)   //priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
        
        sumLogNumerators= newLogConditionalLikelihoodTree +numeratorQ+proposedlogMultinomialProb ;
        
       sumLogDenominators=chain->currentlogConditionalLikelihoodTree + denominatorQ+currentlogMultinomialProb;
        
        
    }
    else//is the oldest pop
    {
        if(!mcmcOptions.noData)
        {
            newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
            sumLogNumerators= newLogConditionalLikelihoodTree ;
             sumLogDenominators=chain->currentlogConditionalLikelihoodTree ;
            if(mcmcOptions.kernelType ==0)//multiplier move
            {
                if (mcmcOptions.doMCMCMoveTimeOriginInputOldestPop)
                {
                    
                    sumLogNumerators += log((pop->timeOriginInput - pop->lowerBoundTimeOriginInput)/ ( pop->oldTimeOriginInput- pop->lowerBoundTimeOriginInput));
                    sumLogNumerators +=log(( mcmcOptions.upperBoundTimeOriginInputOldestPop- pop->timeOriginInput )/ ( mcmcOptions.upperBoundTimeOriginInputOldestPop-pop->oldTimeOriginInput));
                }
                else
                    sumLogNumerators += log((pop->timeOriginSTD)/ ( pop->oldTimeOriginSTD));
            }
        }
        else
        {
            if (pop->timeOriginInput >=pop->lowerBoundTimeOriginInput)
            {
            // sumLogNumerators= pop->LogDensityTime2( pop->timeOriginSTD);
            // sumLogDenominators=pop->LogDensityTime2( pop->oldTimeOriginSTD);
            sumLogNumerators=pop->LogDensityTimeSTDFrom(pop->timeOriginSTD, 0);
            sumLogDenominators=pop->LogDensityTimeSTDFrom( pop->oldTimeOriginSTD,0);
                
             if(mcmcOptions.kernelType ==0)//multiplier move
             {
                 if (mcmcOptions.doMCMCMoveTimeOriginInputOldestPop)
                 {
                    
                     sumLogNumerators += log((pop->timeOriginInput - pop->lowerBoundTimeOriginInput)/ ( pop->oldTimeOriginInput- pop->lowerBoundTimeOriginInput));
                    sumLogNumerators +=log(( mcmcOptions.upperBoundTimeOriginInputOldestPop- pop->timeOriginInput )/ ( mcmcOptions.upperBoundTimeOriginInputOldestPop-pop->oldTimeOriginInput));
                 }
                  else
                    sumLogNumerators += log((pop->timeOriginSTD)/ ( pop->oldTimeOriginSTD));
              }
            }
            else
            {
                 fprintf (stderr, "\n The move is invalid you should not compute acceptance rate \n");
                  sumLogNumerators= log(0);
                  sumLogDenominators= log(0);
                  exit(-1);
            }
        }
    }
    long double difference=sumLogNumerators - sumLogDenominators;
    
    long double LogAcceptanceRate = (difference>=0.0)? 0.0: difference;
    if (isinf(LogAcceptanceRate) || isnan(LogAcceptanceRate))
    {
        if (mcmcOptions.verbose>=1)
          fprintf (stderr, "\n isNan LogAcceptanceRate \n");
        
    }
    return LogAcceptanceRate;
}
void NewTimeOriginOnEdgeforPopulationMove::saveCurrentValue(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
  
    //save information for the population
    //pop->oldFatherPop =  pop->FatherPop;
    //pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
    //pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
    //pop->oldSampleSize = pop->sampleSize;
    //    for(unsigned i=0; i< pop->FatherPop->immigrantsPopOrderedByModelTime.size();i++){
    //
    //         fprintf (stderr, "\n inmigrant %d, time %lf to father population order %d of population of order  %d \n", i , pop->FatherPop->immigrantsPopOrderedByModelTime.at(i).first, pop->FatherPop->order, pop->order);
    //
    //    }
    
    //    for(unsigned i=0; i< pop->CoalescentEventTimes.size();i++){
    //
    //        fprintf (stderr, "\n coalescent %d, time %lf  of population order %d and sample size %d  \n", i , pop->CoalescentEventTimes.at(i), pop->order, pop->sampleSize);
    //
    //    }
    
    pop->FatherPop->oldimmigrantsPopOrderedByModelTime =  pop->FatherPop->immigrantsPopOrderedByModelTime;
    // pop->FatherPop->immigrantsPopOrderedByModelTime.clear();
    
    pop->oldTimeOriginInput = pop->timeOriginInput;
    pop->oldTimeOriginSTD = pop->timeOriginSTD;
    //save information for the other populations
    //    for( i = 0 ; i < chain->numClones; i++)
    //    {
    //        popI=chain->populations[i];
    //       // popI->oldFatherPop =  pop->FatherPop;
    //        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
    //        popI->CoalescentEventTimes.clear();
    //        popI->numCompletedCoalescences=0;
    //        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
    //        popI->immigrantsPopOrderedByModelTime.clear();
    //        popI->numIncomingMigrations=0;
    //    }
}
void NewTimeOriginOnEdgeforPopulationMove::rollbackMove(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    pop->FatherPop->immigrantsPopOrderedByModelTime = pop->FatherPop->oldimmigrantsPopOrderedByModelTime;
    
    //pop->FatherPop->oldimmigrantsPopOrderedByModelTime.clear();
    pop->timeOriginInput = pop->oldTimeOriginInput;
    pop->timeOriginSTD = pop->oldTimeOriginSTD;
    pop->TimeOriginSTD->rollbackValue();
    
    chain->rMRCAPopulation.clear();
    chain->rMRCAPopulation.insert(chain->currentrMRCAPopulation.begin(), chain->currentrMRCAPopulation.end());
    
    chain->proposedrMRCAPopulation.clear();
    //save information for the other populations
    //    for( i = 0 ; i < chain->numClones; i++)
    //    {
    //        popI=chain->populations[i];
    //        //pop->FatherPop = pop->oldFatherPop ;
    //        pop->CoalescentEventTimes = pop->oldCoalescentEventTimes ;
    //        pop->oldCoalescentEventTimes.clear();
    //        pop->immigrantsPopOrderedByModelTime = pop->oldimmigrantsPopOrderedByModelTime;
    //        pop->oldimmigrantsPopOrderedByModelTime.clear();
    //    }
}
void NewTimeOriginOnEdgeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator, Chain* chain)
{
   // Chain *chain=getChain();
    
    std::map<pll_rnode_t*, Population*>  newrmrcaOfPopulation;
    if (pop->timeOriginInput <=0)
        fprintf(stderr, "the time origin input is 0");
    
    chain->chooseNewTimeofOriginOnEdge(pop, mcmcOptions, randomGenerator, NULL);
    //chain->initTimeOriginSTD();
    //    chain->initPopulationMigration();//after setting the timeSTD
    //
    //    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    //    chain->filterSortPopulationsCoalescentEvents();
    sort(pop->FatherPop->immigrantsPopOrderedByModelTime.begin(), pop->FatherPop->immigrantsPopOrderedByModelTime.end(), Population::comparePopulationsPairByTimeOrigin);
    
    //    for (int i = 0; i < pop->FatherPop->immigrantsPopOrderedByModelTime.size(); ++i)
    //        printf("\n ordered migrations: time(father pop units) : %lf, pop order: %d, time of origin input%lf \n", pop->FatherPop->immigrantsPopOrderedByModelTime[i].first,  pop->FatherPop->immigrantsPopOrderedByModelTime[i].second->order , pop->FatherPop->immigrantsPopOrderedByModelTime[i].second->timeOriginInput);
}
long double  NewTimeOriginOnEdgeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
   // Chain *chain=getChain();
    
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions,mcmcOptions);
    //long double currentSumAvailBranchLengths;// chain->sumAvailableBranchLengths(chain->rMRCAPopulation);
    //long double newSumAvailBranchLengths ;// chain->sumAvailableBranchLengths(chain->proposedrMRCAPopulation);
    // double numeratorQ=log(pop->oldrMRCA->length) + log(newSumAvailBranchLengths);
    long double numeratorQ =0;
    //    priorDensityNewTotalEffectivePopulationSize = LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    //double denominatorQ= log(pop->rMRCA->length) + log(currentSumAvailBranchLengths);
    long double denominatorQ= 0;
    //    priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
   long  double sumLogNumerators= newLogConditionalLikelihoodTree +numeratorQ;
    
    long double sumLogDenominators=chain->currentlogConditionalLikelihoodTree + denominatorQ;
    
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>=0.0)?0.0: difference;
   
    return LogAcceptanceRate;
}
//////////////////////////////////////////
NewGlobalScaledMutationRateMove::NewGlobalScaledMutationRateMove(int chainId,std::string nameMove,  Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
long double currentlogConditionalLikelihoodSequences):MCMCmove(chainId, nameMove, mcmcParKernel, affectedParameters, vectorParam, currentlogConditionalLikelihoodTree,
currentlogConditionalLikelihoodSequences)
{
    
    this->pop =pop;
}
void NewGlobalScaledMutationRateMove::move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain)
{
    MCMCmove::move(programOptions,mcmcOptions, rngGsl,rngBoost, chain);
}
void NewGlobalScaledMutationRateMove::saveCurrentValue(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
   // Chain *chain=getChain();
   
    chain->oldtheta = chain->theta;
    chain->thetaPar->saveCurrentValue();
    
    if (!mcmcOptions.noData)
       chain->saveTreeNodeCurrentTimePUnits();
    
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
        
        popI->oldTheta =popI->theta;
        popI->Theta->saveCurrentValue();
        
        if (!mcmcOptions.splitThetaDeltaTmoves)
        {
            
            popI->olddeltaT= popI->deltaT;
            popI->olddelta= popI->delta;
            popI->DeltaT->saveCurrentValue();
        }
      
        popI->oldScaledTimeOriginInput= popI->scaledtimeOriginInput ;
        popI->oldTimeOriginSTD = popI->timeOriginSTD;
        
        popI->TimeOriginSTD->saveCurrentValue();
      
       if (!mcmcOptions.noData)
        {
           
            popI->oldFatherPop = popI->FatherPop;
            popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
            popI->CoalescentEventTimes.clear();
            popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
            popI->immigrantsPopOrderedByModelTime.clear();
        }
        
    }
}
void NewGlobalScaledMutationRateMove::makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,const  gsl_rng *randomGenerator, Chain* chain)
{
    long double kernelPar= chain->thetaPar->getKernelParameter();
    
    long double newValue = chain->thetaPar->stepKernel(randomGenerator);
    
    chain->theta =chain->thetaPar->getParameter()->getCurrentValue();
    assert(chain->theta==newValue);
    if (mcmcOptions.verbose>=2){
      std::cout << "\n proposed new theta " << newValue << " old " <<  chain->oldtheta << " The logPrior is  " << this->mcmcParamKernels[0]->getParameter()->getLogPriorCurrentValue() << std::endl ;
    }
    
    assert(chain->theta >0);
    if (!mcmcOptions.noData)
        chain->updateNodeScaledTimeForRootedTree(chain->theta);
    else{
        if (chain->theta == 0)
          fprintf (stderr, "\n the new theta cannot be 0\n");
      }

    
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
     
        popI->theta = chain->theta * popI->x;
       
        if (!mcmcOptions.splitThetaDeltaTmoves)
        {
            
               kernelPar=  popI->DeltaT->getKernelParameter();
               
               newValue = popI->DeltaT->stepKernel(randomGenerator);
               
               popI->deltaT =popI->DeltaT->getParameter()->getCurrentValue();
            
        
              if (mcmcOptions.verbose>= 3)
             {
                std::cout << "\n proposed deltaT " << popI->deltaT << ", old " << popI->olddeltaT  << std::endl;
              }
            popI->delta = popI->deltaT * popI->x;
            if (mcmcOptions.verbose>=3)
                std::cout << "\n proposed delta " << popI->delta << " old " << popI->olddelta << std::endl;
             }
        
      //  if (mcmcOptions.doInferenceWithoutData==0)
      //  {
            popI->scaledtimeOriginInput = popI->timeOriginInput / chain->theta;
            popI->timeOriginSTD = popI->scaledtimeOriginInput / popI->x;
        
            popI->TimeOriginSTD->setParameterValue(popI->timeOriginSTD);
            
            if (mcmcOptions.verbose>=2)
                std::cout << "\n time origin input " << popI->timeOriginInput << ", scaled by theta " << popI->scaledtimeOriginInput << "  std "<< popI->timeOriginSTD << std::endl;
           
        if (!mcmcOptions.noData)
        {
            popI->numCompletedCoalescences=0;
            popI->CoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime.clear();
            popI->numIncomingMigrations=0;
        }
            //pop->oldFatherPop = pop->FatherPop;
    }
    if (!mcmcOptions.noData)
    {
        chain->InitListPossibleMigrations();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        bool allPopulationPopSizesSet = chain->checkMigrationsOrder();
        if (!allPopulationPopSizesSet)
            fprintf (stderr, "\n The immigrants order is incorrect"  );
    }
}
void NewGlobalScaledMutationRateMove::rollbackMove(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
    assert(chain->oldtheta>0);
    
    chain->theta = chain->oldtheta;
    
    chain->thetaPar->rollbackValue();
 
    // chain->rollbackTreeNodeCurrentTimePUnits();
    if ( !mcmcOptions.noData )
        chain->updateNodeScaledTimeForRootedTree(chain->theta);
    
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
        popI->theta =popI->oldTheta;
        popI->scaledtimeOriginInput = popI->oldScaledTimeOriginInput;
        popI->timeOriginSTD = popI->oldTimeOriginSTD;
        popI->TimeOriginSTD->rollbackValue();
        
        if (!mcmcOptions.splitThetaDeltaTmoves)
        {   popI->deltaT= popI->olddeltaT;
            popI->delta= popI->olddelta;
            popI->DeltaT->rollbackValue();
        }
       
        
        if (!mcmcOptions.noData )
        {
            popI->FatherPop = popI->oldFatherPop;
            popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
            popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
            popI->oldCoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
            popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
            popI->oldimmigrantsPopOrderedByModelTime.clear();
        }
    }
}
long double NewGlobalScaledMutationRateMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
    //Chain *chain=getChain();
    long double temp;
    
    long double from =0.0;
   // fprintf (stderr, "\n>> log conditional Likelihood tree for the new pair  theta and  r %d(old %d),  of the chain %d is = %lf  \n",chain->totalEffectPopSize, chain->oldTotalEffectPopSize,chain->chainNumber,newLogConditionalLikelihoodTree );

    long double priorDensityNewMutationRate=0.0;
    long double priorDensityCurrentMutationRate=0.0;
    
    if(mcmcOptions.priorsType==0)
    { priorDensityNewMutationRate= Distributions::LogUniformDensity(chain->theta, mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
    }
    if(mcmcOptions.priorsType==2)
    {
        priorDensityNewMutationRate =
    Distributions::LogPowerLawDistibutionDensity(chain->theta, mcmcOptions.parameterPowerLawDistributionMutationRate,from);
    }
    if(mcmcOptions.priorsType==1)
    {   priorDensityNewMutationRate =Distributions::LogExponentialDensity(chain->theta, from, mcmcOptions.lambdaExponentialPriorMutationRate);
        
    }
    
    if (isinf(priorDensityNewMutationRate) ||  isnan(priorDensityNewMutationRate))
    {
        if (mcmcOptions.verbose>=1)
            fprintf (stderr, "\n isNan priorDensityNewMutationRate: new mut rate :%.40Lf \n", chain->theta);
        
    }
    
    assert( this->mcmcParamKernels[0]->getParameter()->getLogPriorCurrentValue()==priorDensityNewMutationRate);
    
    if(mcmcOptions.priorsType==0)
    { priorDensityCurrentMutationRate= Distributions::LogUniformDensity(chain->oldtheta, mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
    }
    
    if(mcmcOptions.priorsType==2)
    {
        priorDensityCurrentMutationRate =
    Distributions::LogPowerLawDistibutionDensity(chain->oldtheta,mcmcOptions.parameterPowerLawDistributionMutationRate,from);
    }
    if(mcmcOptions.priorsType==1)
    {  priorDensityCurrentMutationRate =Distributions::LogExponentialDensity(chain->oldtheta,from, mcmcOptions.lambdaExponentialPriorMutationRate);
    }
    
    assert( this->mcmcParamKernels[0]->getParameter()->getOldLogPriorValue()==priorDensityCurrentMutationRate);
    long double priorDensityNewGrowthRate= 0;
    long double priorDensityCurrentGrowthRate= 0;
    
    long double logKernelGrowthRate =0.0;
    
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
        if (!mcmcOptions.splitThetaDeltaTmoves)
        {
            if(mcmcOptions.priorsType==0)
            {
                priorDensityNewGrowthRate +=Distributions::LogUniformDensity(popI->deltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
            }
            if(mcmcOptions.priorsType==2)
            {
                priorDensityNewGrowthRate +=Distributions::LogPowerLawDistibutionDensity(popI->deltaT,mcmcOptions.parameterPowerLawDistributionGrowthRate,from);
            }
            
            if(mcmcOptions.priorsType==1) {
                priorDensityNewGrowthRate +=Distributions::LogExponentialDensity(popI->deltaT,from, mcmcOptions.lambdaExponentialPriorGrowthRate);
            }
            
            if (mcmcOptions.kernelType==0)
            {
                  temp= mcmcParamKernels[0]->getKernel()->logKernel(pop->olddeltaT, pop->deltaT);
                logKernelGrowthRate += log(popI->deltaT / popI->olddeltaT);
                // logKernelGrowthRate += log(popI->olddeltaT / popI->deltaT);
            }
            
            if (isinf(priorDensityNewGrowthRate) ||  isnan(priorDensityNewGrowthRate))
            {
                if (mcmcOptions.verbose>=1)
                    fprintf (stderr, "\n isNan DensityNewGrowthRate: new growth rate :%.40Lf \n", popI->growthRate);
            }
            if(mcmcOptions.priorsType==0)
                priorDensityCurrentGrowthRate += Distributions::LogUniformDensity(popI->olddeltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
            
            if(mcmcOptions.priorsType==2)
                priorDensityCurrentGrowthRate +=Distributions::LogPowerLawDistibutionDensity(popI->olddeltaT,mcmcOptions.parameterPowerLawDistributionGrowthRate,from);
            if(mcmcOptions.priorsType==1)
                priorDensityCurrentGrowthRate +=Distributions::LogExponentialDensity(popI->olddeltaT,from, mcmcOptions.lambdaExponentialPriorGrowthRate);
            
            
        }
    
    }
    long double logkernelTheta=0;
    
    if (mcmcOptions.kernelType==0)
    {
            temp= mcmcParamKernels[0]->getKernel()->logKernel(chain->oldtheta, chain->theta);
        logkernelTheta=log(chain->theta/chain->oldtheta);//log(multiplier factor)
    }
    
    //double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize+ priorDensityNewMutationRate+ logkernelTotalEffectPopulationSize + logkernelMutationRate + priorDensityNewGrowthRate;
    
    //double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize+priorDensityCurrentNewMutationRate+ priorDensityCurrentGrowthRate;
    
   long  double sumLogNumerators=
    priorDensityNewMutationRate+logkernelTheta;
    
    if (!mcmcOptions.splitThetaDeltaTmoves)
       sumLogNumerators+=priorDensityNewGrowthRate+ logKernelGrowthRate;

    long double sumLogDenominators= priorDensityCurrentMutationRate;

    if (!mcmcOptions.splitThetaDeltaTmoves)
        sumLogDenominators+=priorDensityCurrentGrowthRate;
    
   if (!mcmcOptions.noData)
   {  //fprintf (stderr, "\n Log likelihood tree included!\n");
        newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
        sumLogNumerators+= newLogConditionalLikelihoodTree;
        sumLogDenominators+=chain->currentlogConditionalLikelihoodTree;
   }
   
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>0.0)?0.0: difference;
    
    if (isinf(LogAcceptanceRate) ||  isnan(LogAcceptanceRate))
    {
        if (mcmcOptions.verbose>=1)
           fprintf (stderr, "\n isNan  Log acceptance rate\n");
    }
    return LogAcceptanceRate;
}


long double NewGlobalScaledGrowthRateForPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)
{
  //  Chain *chain=getChain();
    
    long double priorDensityNewGrowthRate= 0.0;
    long double priorDensityCurrentGrowthRate= 0.0;
    long double from =0.0;
    //long double temp= mcmcParamKernels[0]->getKernel()->logKernel(pop->olddeltaT, pop->deltaT);
    long double logKernelGrowthRate =0.0;

        if(mcmcOptions.priorsType==0)
        {
            priorDensityNewGrowthRate +=Distributions::LogUniformDensity(pop->deltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
        }
        if(mcmcOptions.priorsType==2)
        {
            priorDensityNewGrowthRate +=Distributions::LogPowerLawDistibutionDensity(pop->deltaT,mcmcOptions.parameterPowerLawDistributionGrowthRate,from);
        }
        
        if(mcmcOptions.priorsType==1) {
            priorDensityNewGrowthRate +=Distributions::LogExponentialDensity(pop->deltaT,from, mcmcOptions.lambdaExponentialPriorGrowthRate);
            //pop->DeltaT->getParameter()->getLogPriorCurrentValue()
        }
        
        if (mcmcOptions.kernelType==0)
        {
            logKernelGrowthRate += log(pop->deltaT / pop->olddeltaT);
            // logKernelGrowthRate += log(popI->olddeltaT / popI->deltaT);
        }
    
        if (isinf(priorDensityNewGrowthRate) ||  isnan(priorDensityNewGrowthRate))
        {
            if (mcmcOptions.verbose>=1)
                fprintf (stderr, "\n isNan DensityNewGrowthRate: new growth rate :%.40Lf \n", pop->growthRate);
        }
        if(mcmcOptions.priorsType==0)
            priorDensityCurrentGrowthRate += Distributions::LogUniformDensity(pop->olddeltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
        
    if(mcmcOptions.priorsType==2){
         priorDensityCurrentGrowthRate +=Distributions::LogPowerLawDistibutionDensity(pop->olddeltaT, mcmcOptions.parameterPowerLawDistributionGrowthRate,from);
    }
    
    if(mcmcOptions.priorsType==1){
          //priorDensityCurrentGrowthRate +=Distributions::LogExponentialDensity(pop->olddeltaT,0, mcmcOptions.lambdaExponentialPriorGrowthRate);
        priorDensityCurrentGrowthRate +=Distributions::LogExponentialDensity(pop->olddeltaT,from, mcmcOptions.lambdaExponentialPriorGrowthRate);
        //pop->DeltaT->getParameter()->getOldLogPriorValue()
    }
    
     long  double sumLogNumerators = priorDensityNewGrowthRate+ logKernelGrowthRate ;
    
    long double sumLogDenominators =  priorDensityCurrentGrowthRate;
   
    if (!mcmcOptions.noData)
    {
        //fprintf (stderr, "\n Log likelihood tree included!\n");
        newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
        sumLogNumerators+= newLogConditionalLikelihoodTree;
        sumLogDenominators+=chain->currentlogConditionalLikelihoodTree;
    }
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>0.0)?0.0: difference;
    
    if (isinf(LogAcceptanceRate) ||  isnan(LogAcceptanceRate))
    {
        if (mcmcOptions.verbose>=1)
            fprintf (stderr, "\n isNan  Log acceptance rate\n");
    }
    return LogAcceptanceRate;
}
void NewGlobalScaledGrowthRateForPopulationMove::rollbackMove(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
   
        pop->deltaT= pop->olddeltaT;
        pop->delta= pop->olddelta;
    
    for(unsigned int i = 0 ; i < mcmcParamKernels.size(); i++){
        mcmcParamKernels[i]->rollbackValue();
        
    }
    for(unsigned int i = 0 ; i < affectedParameters.size(); i++){
        affectedParameters[i]->rollbackValue();
        
    }
        if (!mcmcOptions.noData )
        {
           
            //pop->delta= pop->olddelta;
        }
  
}
void NewGlobalScaledGrowthRateForPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator, Chain* chain)
{

       pop->DeltaT->stepKernel(randomGenerator);
        
        pop->deltaT = this->mcmcParamKernels[0]->getParameter()->getCurrentValue();
    
        if (mcmcOptions.verbose>= 2)
        {fprintf (stderr, "\n proposed deltaT  %.10Lf, old %.10Lf \n",pop->deltaT, pop->olddeltaT );
        }
    
        pop->delta = pop->deltaT * pop->x;
       
        if (mcmcOptions.verbose>=2)
            fprintf (stderr, "\n proposed delta %.10Lf, old %.10Lf \n", pop->delta, pop->olddelta );
      //}
}
void NewGlobalScaledGrowthRateForPopulationMove::saveCurrentValue(ProgramOptions &programOptions,MCMCOptions &mcmcOptions, Chain* chain)
{
  
    //for (unsigned int i = 0; i < chain->numClones; ++i)
   // {
        //auto popI =  chain->populations[i];
       pop->olddeltaT= pop->deltaT;
       pop->olddelta= pop->delta;
    
    for(unsigned int i = 0 ; i < mcmcParamKernels.size(); i++){
        mcmcParamKernels[i]->saveCurrentValue();
    }
    for(unsigned int i = 0 ; i < affectedParameters.size(); i++){
        affectedParameters[i]->saveCurrentValue();
        
    }
        if (!mcmcOptions.noData)
        {
            //pop->olddelta= pop->delta;

            //pop->oldFatherPop = pop->FatherPop;
            
            //pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
            //pop->CoalescentEventTimes.clear();
            //pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
            //pop->immigrantsPopOrderedByModelTime.clear();
        }
    //}
}
void NewTimeOriginOnTreeforPopulationMove::increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{
    //Chain *chain=getChain();check if it is the oldest pop
    //mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop += mcmcOptions.updateLengthMultiplierMCMCMove;
    long double newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() +mcmcOptions.updateLengthMultiplierMCMCMove;
    this->mcmcParamKernels[0]->modifyKernelParameter( newKernelParameter);
}
void NewTimeOriginOnTreeforPopulationMove::decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions )
{
    //Chain *chain=getChain();check if it is the oldest pop
    //mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop -= mcmcOptions.updateLengthMultiplierMCMCMove;
    long double newKernelParameter;
    if (this->mcmcParamKernels[0]->getKernelParameter() >=mcmcOptions.updateLengthMultiplierMCMCMove){
              
              newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() - mcmcOptions.updateLengthMultiplierMCMCMove;
          }
          else{
              
                newKernelParameter= this->mcmcParamKernels[0]->getKernelParameter() - mcmcOptions.updateLengthMultiplierMCMCMove;
               newKernelParameter= newKernelParameter <=0 ? 0.1 :  newKernelParameter;
            //newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() / 2.0;
          }
    //long double newKernelParameter= newParameterValue <=0 ? 0.1 :  newParameterValue;
    this->mcmcParamKernels[0]->modifyKernelParameter( newKernelParameter);
}
void NewGlobalScaledMutationRateMove::increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{
    //mcmcOptions.lengthIntervalMultiplierTheta += mcmcOptions.updateLengthMultiplierMCMCMove;
    
    long double newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() +mcmcOptions.updateLengthMultiplierMCMCMove;
    this->mcmcParamKernels[0]->modifyKernelParameter( newKernelParameter);
    
    if (!mcmcOptions.splitThetaDeltaTmoves){
        newKernelParameter=this->mcmcParamKernels[1]->getKernelParameter() +mcmcOptions.updateLengthMultiplierMCMCMove;
        this->mcmcParamKernels[1]->modifyKernelParameter( newKernelParameter);
      //mcmcOptions.lengthIntervalMultiplierDeltaT += mcmcOptions.updateLengthMultiplierMCMCMove;
    }
}
void NewGlobalScaledMutationRateMove::decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{
    //mcmcOptions.lengthIntervalMultiplierTheta -= mcmcOptions.updateLengthMultiplierMCMCMove;
    long double newKernelParameter;
    if (this->mcmcParamKernels[0]->getKernelParameter() >=mcmcOptions.updateLengthMultiplierMCMCMove){
        
        newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() - mcmcOptions.updateLengthMultiplierMCMCMove;
    }
    else{
        
      newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() / 2.0;
    }
    this->mcmcParamKernels[0]->modifyKernelParameter( newKernelParameter);
    
    if (!mcmcOptions.splitThetaDeltaTmoves){
        
        if (this->mcmcParamKernels[1]->getKernelParameter() >=mcmcOptions.updateLengthMultiplierMCMCMove){
               
               newKernelParameter=this->mcmcParamKernels[1]->getKernelParameter() - mcmcOptions.updateLengthMultiplierMCMCMove;
           }
           else{
             newKernelParameter=this->mcmcParamKernels[1]->getKernelParameter() / 2.0;
           }
        this->mcmcParamKernels[1]->modifyKernelParameter( newKernelParameter);
       //mcmcOptions.lengthIntervalMultiplierDeltaT -= mcmcOptions.updateLengthMultiplierMCMCMove;
    }
}
void NewGlobalScaledGrowthRateForPopulationMove::decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{

    long double newKernelParameter;
    if (this->mcmcParamKernels[0]->getKernelParameter() >=mcmcOptions.updateLengthMultiplierMCMCMove){
           
           newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() - mcmcOptions.updateLengthMultiplierMCMCMove;
       }
       else{
           newKernelParameter= this->mcmcParamKernels[0]->getKernelParameter() - mcmcOptions.updateLengthMultiplierMCMCMove;
           newKernelParameter= newKernelParameter <=0 ? 0.1 :  newKernelParameter;
        // newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() / 2.0;
       }
    this->mcmcParamKernels[0]->modifyKernelParameter( newKernelParameter);
}
void NewGlobalScaledGrowthRateForPopulationMove::increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{
   // mcmcOptions.lengthIntervalMultiplierTheta += mcmcOptions.updateLengthMultiplierMCMCMove;
    long double newKernelParameter=this->mcmcParamKernels[0]->getKernelParameter() +mcmcOptions.updateLengthMultiplierMCMCMove;
    this->mcmcParamKernels[0]->modifyKernelParameter( newKernelParameter);
}
void NewProportionsVectorMove::decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{
    //mcmcOptions.lengthIntervalMultiplierTheta -= mcmcOptions.updateLengthMultiplierMCMCMove;
}
void NewProportionsVectorMove::increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)
{
    //mcmcOptions.lengthIntervalMultiplierTheta += mcmcOptions.updateLengthMultiplierMCMCMove;
}

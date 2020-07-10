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
#include "data_utils.hpp"
#include "utils.hpp"
#include "random.h"
#include <algorithm>

using namespace std;

MCMCmove::MCMCmove(Chain *chain,  std::string nameMove)
{
    this->chain = chain;
    this->nameMove = nameMove;
    this->numberAccept=0;
    this->numberReject=0;
    this->numberAttemps=0;
    this->isInvalidMove=false;
    this->newLogConditionalLikelihoodTree =chain->currentlogConditionalLikelihoodTree;
    this->newLogConditionalLikelihoodSequences=chain->currentlogConditionalLikelihoodSequences;
    
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
Chain * MCMCmove::getChain()
{
    return chain;
}
NewTotalEffectPopSizeMove::NewTotalEffectPopSizeMove(Chain * chain, std::string nameMove ):MCMCmove(chain, nameMove)
{
    
}
NewProportionsVectorMove::NewProportionsVectorMove(Chain * chain, std::string nameMove ):MCMCmove(chain, nameMove)
{
    
}
NewEffectPopSizeMoveForPopulation::NewEffectPopSizeMoveForPopulation(Chain *chain, std::string nameMove, Population *pop):MCMCmove(chain, nameMove)
{
    
}
void NewTotalEffectPopSizeMove::safeCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    Population *popI;
    chain->oldTotalEffectPopSize = chain->totalEffectPopSize;
    for(unsigned int i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->oldeffectPopSize=popI->effectPopSize;
        popI->oldPopSize = popI->popSize;
        popI->oldDeathRate = popI->deathRate;
        popI->oldFatherPop = popI->FatherPop;
        popI->oldGrowthRate = popI->growthRate;
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewTotalEffectPopSizeMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    int i=0;
    Population *popI;
    Chain *chain=getChain();
    double newTotalPopulationSize;
    mcmcOptions.slidingWindowSizeTotalEffectPopSize= 2 * 0.01* (chain->totalEffectPopSize - (chain->totalSampleSize() / mcmcOptions.fixedLambda));
    bool allPopulationPopSizesSet=false;
    
    do {
        numberAttemps++;
        allPopulationPopSizesSet=true;
        newTotalPopulationSize= chain->proposalSlidingWindow(chain->oldTotalEffectPopSize,  mcmcOptions.slidingWindowSizeTotalEffectPopSize);
        
        chain->totalEffectPopSize= newTotalPopulationSize;
        //fprintf (stderr, "\n proposed total effect pop size %lu \n",chain->totalEffectPopSize );
        chain->updateEffectPopSizesCurrentProportionsVector();
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=chain->populations[i];
            popI->growthRate =popI->delta  / popI->effectPopSize;
            popI->popSize=popI->effectPopSize * popI->birthRate;
            popI->deathRate= popI->birthRate - popI->growthRate;
            popI->numCompletedCoalescences=0;
            popI->CoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime.clear();
            popI->numIncomingMigrations=0;
            
            if (popI->popSize < popI->sampleSize)
            {
                allPopulationPopSizesSet=false;
                break;
            }
        }
        if (!allPopulationPopSizesSet)
            fprintf (stderr, "\n The  popSize < sampleSize\n");
        if (!allPopulationPopSizesSet)
            break;
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        allPopulationPopSizesSet = chain->checkMigrationsOrder();
        if (!allPopulationPopSizesSet)
            fprintf (stderr, "\n The order of migrations if not correct, attempts %d\n", numberAttemps);
        if (numberAttemps >mcmcOptions.maxNumberProposalAttempts)
            break;
    }
    while(!allPopulationPopSizesSet);
    
    if (allPopulationPopSizesSet)
        fprintf (stderr, "\n new total effect pop size %lu \n",chain->totalEffectPopSize );
}
void NewTotalEffectPopSizeMove::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    chain->totalEffectPopSize= chain->oldTotalEffectPopSize ;
    int i=0;
    
    Population *popI;
    
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->effectPopSize=popI->oldeffectPopSize;
        popI->popSize = popI->oldPopSize;
        popI->deathRate = popI->oldDeathRate;
        popI->growthRate = popI->oldGrowthRate;
        popI->FatherPop = popI->oldFatherPop;
        popI->CoalescentEventTimes = popI->oldCoalescentEventTimes;
        popI->numCompletedCoalescences=popI->oldCoalescentEventTimes.size();
        popI->oldCoalescentEventTimes.clear();
        popI->immigrantsPopOrderedByModelTime =  popI->oldimmigrantsPopOrderedByModelTime;
        popI->numIncomingMigrations=popI->oldimmigrantsPopOrderedByModelTime.size();
        popI->oldimmigrantsPopOrderedByModelTime.clear();
    }
    
}
long double NewTotalEffectPopSizeMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
    //fprintf (stderr, "\n>> log conditional Likelihood tree for the new total effective population size %d(old %d),  of the chain %d is = %lf  \n",chain->totalEffectPopSize, chain->oldTotalEffectPopSize,chain->chainNumber,newLogConditionalLikelihoodTree );
   // long double priorDensityNewTotalEffectivePopulationSize= LogUniformDensity(chain->totalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    long double priorDensityNewTotalEffectivePopulationSize=LogExponentialDensity(mcmcOptions.lambdaExponentialPriorTotalEffectivePopSize, chain->totalEffectPopSize,mcmcOptions.numberTumorCells / mcmcOptions.fixedLambda);
    
   // long double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
     long double priorDensityCurrentTotalEffectivePopulationSize=LogExponentialDensity(mcmcOptions.lambdaExponentialPriorTotalEffectivePopSize, chain->oldTotalEffectPopSize,mcmcOptions.numberTumorCells / mcmcOptions.fixedLambda);
    
   long double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize;
   long  double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>=0.0)?0.0: difference;
   
    return LogAcceptanceRate;
}

void MCMCmove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    safeCurrentValue(programOptions,mcmcOptions );
    numberAttemps=0;
    makeProposal(programOptions, mcmcOptions);
    
    if (isInvalidMove){
        rollbackMove(programOptions,mcmcOptions);
        this->numberReject++;
        if (mcmcOptions.verbose>=2)
            printf("\n Rejected new move %s \n", nameMove.c_str());
        return;
    }
    
   long  double logAcceptanceRate =computeLogAcceptanceProb(programOptions, mcmcOptions);
  
    if (mcmcOptions.verbose>=1)
      printf("\n logAcceptancerate %.20Lf \n", logAcceptanceRate);
    
    long double randomNumber;
    
    if(mcmcOptions.useGSLRandomGenerator)
        randomNumber= randomUniformFromGsl();
    else
        randomNumber= randomUniformBoost();
    
    long double logRandom=log(randomNumber);
    if (mcmcOptions.verbose>=1)
      printf("\n logRandom %Lf \n", logRandom);
    if (logRandom < logAcceptanceRate )
    {//accept the move
        if (mcmcOptions.verbose>=3)
           printf("\n Accepted new move %s \n", nameMove.c_str());
        this->numberAccept++;
        
        Chain *chain=getChain();
        if ( (!mcmcOptions.noData && chain->currentlogConditionalLikelihoodTree !=newLogConditionalLikelihoodTree))
        {
            chain->currentlogConditionalLikelihoodTree = newLogConditionalLikelihoodTree;
        }
        if ((!mcmcOptions.noData && mcmcOptions.useSequencesLikelihood ==1 && chain->currentlogConditionalLikelihoodSequences !=newLogConditionalLikelihoodSequences))
            chain->currentlogConditionalLikelihoodSequences = newLogConditionalLikelihoodSequences;
    }
    else
    {
        rollbackMove(programOptions,mcmcOptions);
        this->numberReject++;
        if (mcmcOptions.verbose>=2)
          printf("\n Rejected new move %s \n", nameMove.c_str());
    }
}
//////////////////////////////////////////////
void NewProportionsVectorMove::safeCurrentValue(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    chain->oldproportionsVector= chain->proportionsVector;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->oldx= popI->x;
        popI->oldeffectPopSize=popI->effectPopSize;
        popI->oldPopSize = popI->popSize;
        popI->oldFatherPop = popI->FatherPop;
        popI->oldDeathRate = popI->deathRate;
        popI->oldGrowthRate = popI->growthRate;
        popI->oldCoalescentEventTimes = popI->CoalescentEventTimes;
        popI->CoalescentEventTimes.clear();
        popI->oldimmigrantsPopOrderedByModelTime= popI->immigrantsPopOrderedByModelTime;
        popI->immigrantsPopOrderedByModelTime.clear();
    }
}
void NewProportionsVectorMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
    long double *proportionsVectorArray;
    std::vector<long double> effectPopSizesvector(chain->numClones);
    bool allPopulationPopSizesSet=false;
    //double totalSampleSize= chain->totalSampleSize();
    for( i = 0 ; i < chain->numClones; i++)
    {
        effectPopSizesvector.at(i) = chain->oldproportionsVector.at(i) * chain->totalEffectPopSize;
        //fprintf (stderr, "\n old effect population size  at %d: %lf \n",i,effectPopSizesvector.at(i) );
    }
    do {
        allPopulationPopSizesSet=true;
        numberAttemps++;
        randomDirichletFromVector (effectPopSizesvector, chain->proportionsVector);
        //        randomDirichletFromGsl(chain->numClones, proportionsVectorArray, &(chain->proportionsVector[0]));
        chain->updateEffectPopSizesCurrentProportionsVector();
        for( i = 0 ; i < chain->numClones; i++)
        {
            popI=chain->populations[i];
            popI->x= chain->proportionsVector.at(i);
            
            popI->growthRate =popI->delta  / popI->effectPopSize;
            
            popI->deltaT =popI->growthRate * chain->totalEffectPopSize;
            popI->theta = chain->theta * popI->x;
            popI->popSize=popI->effectPopSize * popI->birthRate;
            popI->deathRate= popI->birthRate - popI->growthRate;
            
            popI->numCompletedCoalescences=0;
            popI->CoalescentEventTimes.clear();
            popI->immigrantsPopOrderedByModelTime.clear();
            popI->numIncomingMigrations=0;
            if (popI->popSize < popI->sampleSize)
            {
                allPopulationPopSizesSet=false;
            }
        }
        if (!allPopulationPopSizesSet)
            break;
        chain->initPopulationMigration();//after setting the timeSTD
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
void NewProportionsVectorMove::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    chain->proportionsVector= chain->oldproportionsVector ;
    int i=0;
    Population *popI;
    for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->x= popI->oldx;
        popI->effectPopSize=popI->oldeffectPopSize;
        popI->popSize = popI->oldPopSize;
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
long double NewProportionsVectorMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    // ListClonesAccordingTimeToOrigin(chain->populations);
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions,mcmcOptions);
    std::vector<long double> oldConcentrationParameter(chain->numClones);
    std::vector<long double> proposedConcentrationParameter(chain->numClones);
    for(unsigned i = 0 ; i < chain->numClones; i++){
        oldConcentrationParameter.at(i)  = chain->oldproportionsVector.at(i) * chain->totalEffectPopSize;
        proposedConcentrationParameter.at(i)  = chain->proportionsVector.at(i) * chain->totalEffectPopSize;
    }
    
   long double priorDensityNewProportionsVector =  chain->DirichletDensity( chain->proportionsVector, oldConcentrationParameter , chain->numClones);
    
    long double priorDensityCurrentProportionsVector = chain->DirichletDensity(chain->oldproportionsVector, proposedConcentrationParameter, chain->numClones);
    
    double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewProportionsVector;
    
    double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentProportionsVector;
    
    double logAcceptanceRate = std::min(sumLogNumerators - sumLogDenominators,0.0);
    
    return(logAcceptanceRate);
}
NewGrowthRateMoveForPopulation::NewGrowthRateMoveForPopulation(Chain *chain,std::string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
void NewGrowthRateMoveForPopulation::safeCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
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
void NewGrowthRateMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    double randomDelta = RandomLogUniform(mcmcOptions.Deltafrom, mcmcOptions.Deltato, mcmcOptions.useGSLRandomGenerator);
  
    pop->delta =  randomDelta;
    pop->growthRate =pop->delta  / pop->effectPopSize;
    pop->deathRate= pop->birthRate - pop->growthRate;
    
    chain->initPopulationMigration();//after setting the timeSTD
    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    chain->filterSortPopulationsCoalescentEvents();
}
void NewGrowthRateMoveForPopulation::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
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
long double NewGrowthRateMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
    
    //fprintf (stderr, "\n>> New log conditional Likelihood tree after the new total growth rate for population %d,  of the chain %d is = %lf  \n", pop->index, chain->chainNumber, newLogConditionalLikelihoodTree );
    
    long double priorDensityNewScaledGrowthRate =  LogUniformDensity(pop->delta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    long double priorDensitScaledGrowthRate =  LogUniformDensity(pop->olddelta, mcmcOptions.Deltafrom, mcmcOptions.Deltato);
    
    long double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewScaledGrowthRate;
    
    long double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensitScaledGrowthRate;
    
    
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>=0.0)?0.0: difference;
    
    return LogAcceptanceRate;
}
void NewEffectPopSizeMoveForPopulation::safeCurrentValue(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{    Chain *chain=getChain();
    pop->oldeffectPopSize= pop->effectPopSize;
    chain->oldTotalEffectPopSize =chain->totalEffectPopSize;
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
void NewEffectPopSizeMoveForPopulation::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{  Chain *chain=getChain();
    pop->effectPopSize= pop->oldeffectPopSize;
    chain->totalEffectPopSize =chain->oldTotalEffectPopSize;
    Population *popI;
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
void NewEffectPopSizeMoveForPopulation::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{

    Chain *chain=getChain();
    
    bool populationPopSizesSet=false;
    double newEffectPopulationSize;
    //int totalSampleSize = chain->totalSampleSize();
    
    double slidingWindowSize= 2* 0.1* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
    do {
        newEffectPopulationSize= chain->proposalSlidingWindow(pop->oldeffectPopSize,  slidingWindowSize);
        
        pop->effectPopSize= newEffectPopulationSize;
        
        chain->totalEffectPopSize =  chain->totalEffectPopSize + pop->effectPopSize - pop->oldeffectPopSize;
        
        chain->updateEffectPopSizesCurrentProportionsVector();
        
        pop->growthRate =pop->delta  / pop->effectPopSize;
        pop->popSize=pop->effectPopSize * pop->birthRate;
        pop->deathRate= pop->birthRate - pop->growthRate;
        if (pop->popSize < pop->sampleSize)
        {
            populationPopSizesSet=false;
        }
        if (!populationPopSizesSet)
            break;
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        populationPopSizesSet = chain->checkMigrationsOrder();
    }
    while(!populationPopSizesSet);
    
}
long double NewEffectPopSizeMoveForPopulation::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions){
    
    Chain *chain=getChain();
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions,mcmcOptions);
    
    // double slidingWindowSize= 2* (pop->oldeffectPopSize - pop->birthRate * pop->sampleSize);
    
    long double priorDensityNewEffectivePopulationSize= LogUniformDensity(pop->effectPopSize,  mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    long double priorDensityCurrentTotalEffectivePopulationSize= LogUniformDensity(chain->oldTotalEffectPopSize, mcmcOptions.totalEffectPopSizefrom, mcmcOptions.totalEffectPopSizeto);
    
    long double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewEffectivePopulationSize;
    
   long  double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize;
    
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>=0.0)?0.0: difference;
    
    return LogAcceptanceRate;
}
void NewTotalEffectPopSizeMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewProportionsVectorMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewGrowthRateMoveForPopulation::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewEffectPopSizeMoveForPopulation::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewTimeOriginOnTreeforPopulationMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewTimeOriginOnEdgeforPopulationMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewGlobalScaledGrowthRateForPopulationMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
NewGlobalScaledGrowthRateForPopulationMove::NewGlobalScaledGrowthRateForPopulationMove(Chain *chain,std::string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
NewTimeOriginOnTreeforPopulationMove::NewTimeOriginOnTreeforPopulationMove(Chain *chain,std::string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
NewTimeOriginOnEdgeforPopulationMove::NewTimeOriginOnEdgeforPopulationMove(Chain *chain,std::string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
void NewTimeOriginOnTreeforPopulationMove::safeCurrentValue(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    //save information for the population
    bool isOldestPop= chain->isOldestPopulation(pop, chain->rMRCAPopulation);
    
    pop->oldrMRCA =  pop->rMRCA ;//this only changes to population pop
    
    pop->oldTimeOriginInput = pop->timeOriginInput;
    pop->oldScaledTimeOriginInput = pop->scaledtimeOriginInput;
    pop->oldTimeOriginSTD = pop->timeOriginSTD;
    //save information for the other populations
   
         Population *popI;
        for(int i = 0 ; i < chain->numClones; i++)
        {
            popI=chain->populations[i];
      
            popI->oldSampleSize = popI->sampleSize;
           // popI->oldeffectPopSize = popI->effectPopSize;
            popI->oldTimeOriginSTD = popI->timeOriginSTD;
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
void NewTimeOriginOnTreeforPopulationMove::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
     bool isOldestPop= chain->isOldestPopulation(pop, chain->rMRCAPopulation);
    //save information for the population
   
    pop->rMRCA = pop->oldrMRCA ;
    
    pop->timeOriginInput = pop->oldTimeOriginInput;

    pop->scaledtimeOriginInput = pop->oldScaledTimeOriginInput ;
    pop->timeOriginSTD = pop->oldTimeOriginSTD ;
    //save information for the other populations
    if(programOptions.numClones >1 )
    {
        for( i = 0 ; i < chain->numClones; i++)
    {
        popI=chain->populations[i];
        popI->sampleSize = popI->oldSampleSize;
        //  popI->effectPopSize = popI->oldeffectPopSize;
        if (popI != pop)
          popI->timeOriginSTD = popI->oldTimeOriginSTD;
        
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
void NewTimeOriginOnTreeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    //first compute the list of 3 adjacent edges
    //build a cummulative list with the 3 adjacent edges
    //then select at random number between 0 and the sum of the  3 adjacent edges
    //recompute the sample sizes
    // recompute the coalescent and migration events
    //bool existsZeroSampleSizePop=false;
    //double alpha[chain->numClones];
    //int numberPoints = chain->numClones -1;
    std::map<pll_rnode_t*, std::vector<Population*>>  newrmrcaOfPopulation;
    map<pll_rnode_t*, std::vector<Population*>>::iterator it;
    
    bool isOldestPop = chain->isOldestPopulation(pop, chain->rMRCAPopulation) ;
    
    if (!isOldestPop)
    {
        chain->proposedrMRCAPopulation=  chain->chooseAvailableEdgeOnRootedTreeForPopulation(pop, chain->rMRCAPopulation, programOptions.healthyTipLabel);
    }
    else
    {
        chain->proposedrMRCAPopulation=chain->chooseNewTimeofOriginOnEdge(pop,mcmcOptions);
         if (pop->timeOriginInput <pop->lowerBoundTimeOriginInput)
             isInvalidMove = true;
    }
    bool stillTheOldestPop=  chain->isOldestPopulation(pop, chain->proposedrMRCAPopulation) ;
    
    if (!stillTheOldestPop)//not the oldest population any more
    {
        chain->currentrMRCAPopulation.clear();
        chain->currentrMRCAPopulation.insert(chain->rMRCAPopulation.begin(), chain->rMRCAPopulation.end());
        
        chain->rMRCAPopulation.clear();
        chain->rMRCAPopulation.insert(chain->proposedrMRCAPopulation.begin(), chain->proposedrMRCAPopulation.end());
        
        chain->initPopulationsSampleSizes( chain->proposedrMRCAPopulation, programOptions.healthyTipLabel);
        
        chain->initTimeOriginSTD(mcmcOptions);
        chain->initPopulationMigration();//after setting the timeSTD
        
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->proposedrMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        
    }
    else{//still the oldest population
        chain->currentrMRCAPopulation.clear();
        chain->currentrMRCAPopulation.insert(chain->rMRCAPopulation.begin(), chain->rMRCAPopulation.end());
        
        chain->rMRCAPopulation.clear();
        chain->rMRCAPopulation.insert(chain->proposedrMRCAPopulation.begin(), chain->proposedrMRCAPopulation.end());
        for (unsigned int i = 0; i < chain->numClones; ++i)
        {
            auto popI =  chain->populations[i];
            popI->sampleSize= popI->oldSampleSize ;
            if (popI != pop)
               popI->timeOriginSTD = popI->oldTimeOriginSTD;
            
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
long double  NewTimeOriginOnTreeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    TreeNode *u, *v;
    std::vector<pair<double, pll_tree_edge_t *> > currentAvailableEdges;
    std::vector<pair<double, pll_tree_edge_t *> > proposedAvailableEdges;
    
    
    
    bool isOldestPop = chain->isOldestPopulation(pop, chain->rMRCAPopulation);
    long double sumLogNumerators=0.0 ;
    
   long double sumLogDenominators=0.0;
    
    if (!isOldestPop)
    {
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
        long double currentlogMultinomialProb= logMultinomialProbability(chain->numClones, currentProportionsArray, currentSampleSizesArray);
        
         long double proposedlogMultinomialProb=logMultinomialProbability(chain->numClones, currentProportionsArray, proposedSampleSizesArray);
        
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
        }
        else
        {
            if (pop->timeOriginInput >=pop->lowerBoundTimeOriginInput)
            {
             sumLogNumerators= pop->LogDensityTimeSTDFrom( pop->timeOriginSTD, 0);
             sumLogDenominators=pop->LogDensityTimeSTDFrom( pop->oldTimeOriginSTD, 0);
                
             if(mcmcOptions.kernelType ==0)//multiplier move
             {
                    sumLogNumerators += log(( mcmcOptions.upperBoundTimeOriginInputOldestPop-pop->timeOriginInput)/ ( mcmcOptions.upperBoundTimeOriginInputOldestPop-pop->oldTimeOriginInput));
                    sumLogNumerators += log((pop->timeOriginInput - pop->lowerBoundTimeOriginInput)/ (pop->oldTimeOriginInput-pop->lowerBoundTimeOriginInput));
              }
            }
            else
            {
                  sumLogNumerators= log(0);
                  sumLogDenominators= log(0);
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
void NewTimeOriginOnEdgeforPopulationMove::safeCurrentValue(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    int i=0;
    Population *popI;
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
void NewTimeOriginOnEdgeforPopulationMove::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
 
    Population *popI;
    //save information for the population
    
    pop->FatherPop->immigrantsPopOrderedByModelTime = pop->FatherPop->oldimmigrantsPopOrderedByModelTime;
    
    //pop->FatherPop->oldimmigrantsPopOrderedByModelTime.clear();
    pop->timeOriginInput = pop->oldTimeOriginInput;
    pop->timeOriginSTD = pop->oldTimeOriginSTD;
    
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
void NewTimeOriginOnEdgeforPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    std::map<pll_rnode_t*, Population*>  newrmrcaOfPopulation;
    if (pop->timeOriginInput <=0)
        fprintf(stderr, "the time origin input is 0");
    
    chain->chooseNewTimeofOriginOnEdge(pop, mcmcOptions);
    //chain->initTimeOriginSTD();
    //    chain->initPopulationMigration();//after setting the timeSTD
    //
    //    chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
    //    chain->filterSortPopulationsCoalescentEvents();
    sort(pop->FatherPop->immigrantsPopOrderedByModelTime.begin(), pop->FatherPop->immigrantsPopOrderedByModelTime.end(), Population::comparePopulationsPairByTimeOrigin);
    
    //    for (int i = 0; i < pop->FatherPop->immigrantsPopOrderedByModelTime.size(); ++i)
    //        printf("\n ordered migrations: time(father pop units) : %lf, pop order: %d, time of origin input%lf \n", pop->FatherPop->immigrantsPopOrderedByModelTime[i].first,  pop->FatherPop->immigrantsPopOrderedByModelTime[i].second->order , pop->FatherPop->immigrantsPopOrderedByModelTime[i].second->timeOriginInput);
}
long double  NewTimeOriginOnEdgeforPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions,mcmcOptions);
    long double currentSumAvailBranchLengths;// chain->sumAvailableBranchLengths(chain->rMRCAPopulation);
    long double newSumAvailBranchLengths ;// chain->sumAvailableBranchLengths(chain->proposedrMRCAPopulation);
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
NewGlobalScaledMutationRateMove::NewGlobalScaledMutationRateMove(Chain *chain,std::string nameMove,  Population *pop):MCMCmove(chain, nameMove)
{
    this->pop =pop;
}
void NewGlobalScaledMutationRateMove::move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    MCMCmove::move(programOptions,mcmcOptions);
}
void NewGlobalScaledMutationRateMove::safeCurrentValue(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
   
 
    chain->oldtheta = chain->theta;
    if (!mcmcOptions.noData)
      chain->safeTreeNodeCurrentTimePUnits();
    
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
        popI->olddeltaT= popI->deltaT;
        popI->oldTheta =popI->theta;
        popI->olddelta= popI->delta;
        popI->oldScaledTimeOriginInput= popI->scaledtimeOriginInput ;
        popI->oldTimeOriginSTD = popI->timeOriginSTD;
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
void NewGlobalScaledMutationRateMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    //mcmcOptions.paramMultiplierEffectPopSize = parameterMultiplierMCMCmove (mcmcOptions.lengthIntervalMultiplier);
   
    //generate a uniform variable between 0 and 1
    double randomNumber;
    if(mcmcOptions.useGSLRandomGenerator)
       randomNumber= randomUniformFromGsl();
    else
       randomNumber= randomUniformBoost();
    
    //mcmcOptions.paramMultiplierTheta =parameterMultiplierMCMCmove (mcmcOptions.lengthIntervalMultiplier);
    long double m2 = exp(2 * log(mcmcOptions.paramMultiplierTheta) *(randomNumber -0.5));
    
    if (mcmcOptions.kernelType==0)
    {chain->theta = m2 * chain->theta;}
    
    if (mcmcOptions.kernelType==1)
    {
        chain->theta = randomNormalGreaterThan(chain->mutationRate, mcmcOptions.sigmaNormalKernelMutationRate,0);
    }
    if (mcmcOptions.verbose>=2)
      fprintf (stderr, "\n proposed new theta %.10Lf, old %.10Lf \n", chain->theta, chain->oldtheta );
    
    if (chain->theta != 0 &&  (!mcmcOptions.noData))
        chain->updateNodeScaledTimeForRootedTree(chain->theta);
    else{
        if (chain->theta == 0)
          fprintf (stderr, "\n the new theta cannot be 0\n");
      }
    long double m3;
    //mcmcOptions.paramMultiplierGrowthRate = parameterMultiplierMCMCmove(mcmcOptions.lengthIntervalMultiplier);
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
     
            popI->theta = chain->theta * popI->x;
       
        if (!mcmcOptions.splitThetaDeltaTmoves)
        {
           if(mcmcOptions.useGSLRandomGenerator)
              randomNumber= randomUniformFromGsl();
           else
             randomNumber= randomUniformBoost();
        
            m3 = exp(2 * log(mcmcOptions.paramMultiplierGrowthRate) *(randomNumber -0.5));
            //  popI->delta =  chain->proposalSlidingWindow(popI->Delta, mcmcOptions.slidingWindowSizeGrowtRate);
            //popI->delta =  chain->proposalSlidingWindow(popI->Delta, 0.1 * popI->oldGrowthRate);
            if (mcmcOptions.kernelType==1)
            {
               popI->deltaT = randomNormalGreaterThan(popI->deltaT, mcmcOptions.sigmaNormalKernelGrowthRate,0);
            }
            if (mcmcOptions.kernelType==0)
            {
                popI->deltaT = m3 * popI->deltaT ;
             }
            if (mcmcOptions.verbose>= 3)
             {
               fprintf (stderr, "\n proposed deltaT  %.10Lf, old %.10Lf \n",popI->deltaT, popI->olddeltaT );
             }
             popI->delta = popI->deltaT * popI->x;
             if (mcmcOptions.verbose>=3)
               fprintf (stderr, "\n proposed delta %.10Lf, old %.10Lf \n", popI->delta, popI->olddelta );
        }
      //  if (mcmcOptions.doInferenceWithoutData==0)
      //  {
            popI->scaledtimeOriginInput = popI->timeOriginInput / popI->theta;
            popI->timeOriginSTD = popI->scaledtimeOriginInput / popI->x;
            
            if (mcmcOptions.verbose>=2)
                fprintf (stderr, "\n time origin input   %.10Lf  , scaled by theta %.10Lf, std %.10Lf \n",  popI->timeOriginInput,  popI->scaledtimeOriginInput, popI->timeOriginSTD  );
           
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
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        bool allPopulationPopSizesSet = chain->checkMigrationsOrder();
    }
}
void NewGlobalScaledMutationRateMove::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    chain->theta = chain->oldtheta;
 
    // chain->rollbackTreeNodeCurrentTimePUnits();
    if (chain->theta != 0 && !mcmcOptions.noData )
        chain->updateNodeScaledTimeForRootedTree(chain->theta);
    
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
        popI->theta =popI->oldTheta;
        popI->scaledtimeOriginInput = popI->oldScaledTimeOriginInput;
        popI->timeOriginSTD = popI->oldTimeOriginSTD;
  
         if (!mcmcOptions.splitThetaDeltaTmoves)
            popI->deltaT= popI->olddeltaT;
        if (!mcmcOptions.splitThetaDeltaTmoves)
            popI->delta= popI->olddelta;
        
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
long double NewGlobalScaledMutationRateMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
   // fprintf (stderr, "\n>> log conditional Likelihood tree for the new pair  theta and  r %d(old %d),  of the chain %d is = %lf  \n",chain->totalEffectPopSize, chain->oldTotalEffectPopSize,chain->chainNumber,newLogConditionalLikelihoodTree );

    long double priorDensityNewMutationRate=0.0;
    long double priorDensityCurrentMutationRate=0.0;
    
    if(mcmcOptions.priorsType==0)
    { priorDensityNewMutationRate= LogUniformDensity(chain->theta, mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
    }
    if(mcmcOptions.priorsType==2)
    {
        priorDensityNewMutationRate =
    LogPowerLawDistibutionDensity(mcmcOptions.parameterPowerLawDistributionMutationRate, chain->theta,0);
    }
    if(mcmcOptions.priorsType==1)
    {   priorDensityNewMutationRate =LogExponentialDensity(mcmcOptions.lambdaExponentialPriorMutationRate,chain->theta,0);
        
    }
    
    if (isinf(priorDensityNewMutationRate) ||  isnan(priorDensityNewMutationRate))
    {
        if (mcmcOptions.verbose>=1)
            fprintf (stderr, "\n isNan priorDensityNewMutationRate: new mut rate :%.40Lf \n", chain->theta);
        
    }
    
    if(mcmcOptions.priorsType==0)
    { priorDensityCurrentMutationRate= LogUniformDensity(chain->oldtheta, mcmcOptions.MutRatefrom, mcmcOptions.MutRateto);
    }
    
    if(mcmcOptions.priorsType==2)
    {
        priorDensityCurrentMutationRate =
    LogPowerLawDistibutionDensity(mcmcOptions.parameterPowerLawDistributionMutationRate, chain->oldtheta,0);
    }
    if(mcmcOptions.priorsType==1)
    {  priorDensityCurrentMutationRate =LogExponentialDensity(mcmcOptions.lambdaExponentialPriorMutationRate,chain->oldtheta,0);
    }
    
    long double priorDensityNewGrowthRate= 0;
    long double priorDensityCurrentGrowthRate= 0;
    
    long double logKernelGrowthRate =0.0;
    for (unsigned int i = 0; i < chain->numClones; ++i)
    {
        auto popI =  chain->populations[i];
        if(mcmcOptions.priorsType==0)
        {
            priorDensityNewGrowthRate +=LogUniformDensity(popI->deltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
            
        }
        if(mcmcOptions.priorsType==2)
        {
            priorDensityNewGrowthRate +=
            LogPowerLawDistibutionDensity(mcmcOptions.parameterPowerLawDistributionGrowthRate, popI->deltaT,0);
            
        }
        
        if(mcmcOptions.priorsType==1) {
            priorDensityNewGrowthRate +=LogExponentialDensity(mcmcOptions.lambdaExponentialPriorGrowthRate,popI->deltaT,0);
        }
        
        if (mcmcOptions.kernelType==0)
           {
               logKernelGrowthRate += log(popI->deltaT / popI->olddeltaT);
              // logKernelGrowthRate += log(popI->olddeltaT / popI->deltaT);
         }
        
        if (isinf(priorDensityNewGrowthRate) ||  isnan(priorDensityNewGrowthRate))
        {
           if (mcmcOptions.verbose>=1)
             fprintf (stderr, "\n isNan DensityNewGrowthRate: new growth rate :%.40Lf \n", popI->growthRate);
            
        }
        
        if(mcmcOptions.priorsType==0)
            priorDensityCurrentGrowthRate += LogUniformDensity(popI->olddeltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
        
        if(mcmcOptions.priorsType==2)
            priorDensityCurrentGrowthRate +=LogPowerLawDistibutionDensity(mcmcOptions.parameterPowerLawDistributionGrowthRate, popI->olddeltaT,0);
        if(mcmcOptions.priorsType==1)
            priorDensityCurrentGrowthRate +=LogExponentialDensity(mcmcOptions.lambdaExponentialPriorGrowthRate,popI->olddeltaT,0);
    }
    

    long double logkernelTheta=0;
    
    if (mcmcOptions.kernelType==0)
    {
        logkernelTheta=log(chain->theta/chain->oldtheta);//log(multiplier factor)
    }
    
    //double sumLogNumerators= newLogConditionalLikelihoodTree +priorDensityNewTotalEffectivePopulationSize+ priorDensityNewMutationRate+ logkernelTotalEffectPopulationSize + logkernelMutationRate + priorDensityNewGrowthRate;
    
    //double sumLogDenominators=chain->currentlogConditionalLikelihoodTree +priorDensityCurrentTotalEffectivePopulationSize+priorDensityCurrentNewMutationRate+ priorDensityCurrentGrowthRate;
    
   long  double sumLogNumerators=
    priorDensityNewMutationRate+logkernelTheta;
    
    sumLogNumerators+=priorDensityNewGrowthRate+ logKernelGrowthRate;

    long double sumLogDenominators= priorDensityCurrentMutationRate;
    sumLogDenominators+=priorDensityCurrentGrowthRate;
    
   if (!mcmcOptions.noData)
   {  fprintf (stderr, "\n Log likelihood tree included!\n");
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


long double NewGlobalScaledGrowthRateForPopulationMove::computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    
    
    
    long double priorDensityNewGrowthRate= 0;
    long double priorDensityCurrentGrowthRate= 0;
    
    long double logKernelGrowthRate =0.0;
   // for (unsigned int i = 0; i < chain->numClones; ++i)
   // {
       // auto popI =  chain->populations[i];
        if(mcmcOptions.priorsType==0)
        {
            priorDensityNewGrowthRate +=LogUniformDensity(pop->deltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
            
        }
        if(mcmcOptions.priorsType==2)
        {
            priorDensityNewGrowthRate +=
            LogPowerLawDistibutionDensity(mcmcOptions.parameterPowerLawDistributionGrowthRate, pop->deltaT,0);
            
        }
        
        if(mcmcOptions.priorsType==1) {
            priorDensityNewGrowthRate +=LogExponentialDensity(mcmcOptions.lambdaExponentialPriorGrowthRate,pop->deltaT,0);
            
            
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
            priorDensityCurrentGrowthRate += LogUniformDensity(pop->olddeltaT, mcmcOptions.GrowthRatefrom, mcmcOptions.GrowthRateto);
        
        if(mcmcOptions.priorsType==2)
            priorDensityCurrentGrowthRate +=LogPowerLawDistibutionDensity(mcmcOptions.parameterPowerLawDistributionGrowthRate, pop->olddeltaT,0);
        if(mcmcOptions.priorsType==1)
            priorDensityCurrentGrowthRate +=LogExponentialDensity(mcmcOptions.lambdaExponentialPriorGrowthRate,pop->olddeltaT,0);
   // }
    
     long  double sumLogNumerators = priorDensityNewGrowthRate+ logKernelGrowthRate ;
    
    long double sumLogDenominators =  priorDensityCurrentGrowthRate;
   
   // if (mcmcOptions.doInferenceWithoutData ==0)
   // {
        //fprintf (stderr, "\n Log likelihood tree included!\n");
        newLogConditionalLikelihoodTree= chain->LogConditionalLikelihoodTree(programOptions, mcmcOptions);
        sumLogNumerators+= newLogConditionalLikelihoodTree;
        sumLogDenominators+=chain->currentlogConditionalLikelihoodTree;
   // }
    long double difference=sumLogNumerators - sumLogDenominators;
    long double LogAcceptanceRate = (difference>0.0)?0.0: difference;
    
    if (isinf(LogAcceptanceRate) ||  isnan(LogAcceptanceRate))
    {
        if (mcmcOptions.verbose>=1)
            fprintf (stderr, "\n isNan  Log acceptance rate\n");
    }
    return LogAcceptanceRate;
}
void NewGlobalScaledGrowthRateForPopulationMove::rollbackMove(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();

   // for (unsigned int i = 0; i < chain->numClones; ++i)
  //  {
   //     auto popI =  chain->populations[i];
        // popI->growthRate =popI->oldGrowthRate;
        // popI->effectPopSize = popI->oldeffectPopSize;
        // popI->popSize = popI->oldPopSize;
        // popI->deathRate = popI->oldDeathRate;
        pop->deltaT= pop->olddeltaT;
    
        if (!mcmcOptions.noData )
        {
            pop->theta =pop->oldTheta;
            pop->delta= pop->olddelta;
            pop->FatherPop = pop->oldFatherPop;
            pop->scaledtimeOriginInput = pop->oldScaledTimeOriginInput;
            pop->timeOriginSTD = pop->oldTimeOriginSTD;
            pop->CoalescentEventTimes = pop->oldCoalescentEventTimes;
            pop->numCompletedCoalescences=pop->oldCoalescentEventTimes.size();
            pop->oldCoalescentEventTimes.clear();
            pop->immigrantsPopOrderedByModelTime =  pop->oldimmigrantsPopOrderedByModelTime;
            pop->numIncomingMigrations=pop->oldimmigrantsPopOrderedByModelTime.size();
            pop->oldimmigrantsPopOrderedByModelTime.clear();
        }
  //  }
}
void NewGlobalScaledGrowthRateForPopulationMove::makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    double randomNumber;
    long double m3;
    mcmcOptions.paramMultiplierGrowthRate = parameterMultiplierMCMCmove(mcmcOptions.lengthIntervalMultiplier);
    
    //for (unsigned int i = 0; i < chain->numClones; ++i)
    //{
    //    auto popI =  chain->populations[i];
        
        if (!mcmcOptions.noData)
        {
            pop->theta = chain->theta * pop->x;
        }
        if(mcmcOptions.useGSLRandomGenerator)
            randomNumber= randomUniformFromGsl();
        else
            randomNumber= randomUniformBoost();
        m3 = exp(2 * log(mcmcOptions.paramMultiplierGrowthRate) *(randomNumber -0.5));
        //  popI->delta =  chain->proposalSlidingWindow(popI->Delta, mcmcOptions.slidingWindowSizeGrowtRate);
        //popI->delta =  chain->proposalSlidingWindow(popI->Delta, 0.1 * popI->oldGrowthRate);
        if (mcmcOptions.kernelType==1)
        {
            pop->deltaT = randomNormalGreaterThan(pop->deltaT, mcmcOptions.sigmaNormalKernelGrowthRate,0);
        }
        if (mcmcOptions.kernelType==0)
        {
            pop->deltaT = m3 * pop->deltaT ;
        }
        if (mcmcOptions.verbose>= 3)
        {fprintf (stderr, "\n proposed deltaT  %.10Lf, old %.10Lf \n",pop->deltaT, pop->olddeltaT );
        }
        if (!mcmcOptions.noData)
        {
            pop->delta = pop->deltaT * pop->x;
        }
        if (mcmcOptions.verbose>=3)
            fprintf (stderr, "\n proposed delta %.10Lf, old %.10Lf \n", pop->delta, pop->olddelta );
        if (!mcmcOptions.noData)
        {
            pop->scaledtimeOriginInput = pop->timeOriginInput / pop->theta;
            pop->timeOriginSTD = pop->scaledtimeOriginInput / pop->x;
            if (mcmcOptions.verbose>=2)
                fprintf (stderr, "\n time origin input   %.10Lf  , scaled by theta %.10Lf, std %.10Lf \n",  pop->timeOriginInput,  pop->scaledtimeOriginInput, pop->timeOriginSTD  );
            if (pop->popSize < pop->sampleSize)
            {
                isInvalidMove = true;
                //            allPopulationPopSizesSet=false;
                //break;
            }
            pop->numCompletedCoalescences=0;
            pop->CoalescentEventTimes.clear();
            pop->immigrantsPopOrderedByModelTime.clear();
            pop->numIncomingMigrations=0;
            //pop->oldFatherPop = pop->FatherPop;
        }
    //}
    if (!mcmcOptions.noData)
    {
        chain->initPopulationMigration();//after setting the timeSTD
        chain->initPopulationsCoalescentAndMigrationEventsFromRootedTree(chain->rMRCAPopulation, programOptions.healthyTipLabel);
        chain->filterSortPopulationsCoalescentEvents();
        bool allPopulationPopSizesSet = chain->checkMigrationsOrder();
    }
}
void NewGlobalScaledGrowthRateForPopulationMove::safeCurrentValue(ProgramOptions &programOptions,MCMCoptions &mcmcOptions)
{
    Chain *chain=getChain();
    //for (unsigned int i = 0; i < chain->numClones; ++i)
   // {
        //auto popI =  chain->populations[i];
        pop->olddeltaT= pop->deltaT;
    
        if (!mcmcOptions.noData)
        {
            pop->oldTheta =pop->theta;
            pop->olddelta= pop->delta;
            pop->oldScaledTimeOriginInput= pop->scaledtimeOriginInput ;
            pop->oldTimeOriginSTD = pop->timeOriginSTD;
            pop->oldFatherPop = pop->FatherPop;
            
            pop->oldCoalescentEventTimes = pop->CoalescentEventTimes;
            pop->CoalescentEventTimes.clear();
            pop->oldimmigrantsPopOrderedByModelTime= pop->immigrantsPopOrderedByModelTime;
            pop->immigrantsPopOrderedByModelTime.clear();
        }
    //}
}

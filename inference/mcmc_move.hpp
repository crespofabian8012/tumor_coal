/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the MCMC C++ library.
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
 * MCMC move class
 */

#ifndef mcmc_move_hpp
#define mcmc_move_hpp

#include "mcmc_chain.hpp"
#include "mcmc_parameter.hpp"


//class Chain;
//class MCMCoptions;
//class ProgramOptions;
//class Population;


class MCMCParameterWithKernel;
class MCMCVectorParameter;

class MCMCmove{
    Chain *chain;
    std::string nameMove;
protected:
    int numberAccept;
    int numberReject;
    int numberAttemps;
    vector<MCMCParameterWithKernel *> mcmcParamKernels;
    vector<IMCMCParameter<long double> *> affectedParameters;
    
    //MCMCVectorParameter * vectorParam;
    MCMCVectorParameterWithKernel *vectorParam;
    
    MCMCParameter<vector<long double>> *affectedVectorParameters;
    //MCMCParameter<vector<long double>> sampleSizesVector;
    
    long double newLogConditionalLikelihoodTree;
    long double newLogConditionalLikelihoodSequences;
public:
    bool isInvalidMove;
    std::string name();
    MCMCmove(Chain *chain, std::string nameMove, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam );
    virtual void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions,gsl_rng *randomGenerator)=0;
    virtual void rollbackMove( ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    virtual void saveCurrentValue( ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    virtual long  double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    virtual void increaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    virtual void decreaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    double numberAccepted();
    void  resetNumberAccepted();
    void  resetNumberRejected();
    double numberRejected();
    Chain * getChain();
public:
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions,gsl_rng *randomGenerator);
    
    // virtual ~MCMCmove()=0;
};

class NewTotalEffectPopSizeMove:public MCMCmove{
public:
    NewTotalEffectPopSizeMove(Chain *chain, std::string nameMove, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions,gsl_rng *randomGenerator);
};

class NewProportionsVectorMove:public MCMCmove{
public:
    NewProportionsVectorMove(Chain *chain, std::string nameMove, vector<MCMCParameterWithKernel *> &mcmcParKernel, vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions,gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void increaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
class NewGrowthRateMoveForPopulation:public MCMCmove{
    Population *pop;
public:
    NewGrowthRateMoveForPopulation(Chain *chain, std::string nameMove,Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions,gsl_rng *randomGenerator);
};

class NewEffectPopSizeMoveForPopulation:public MCMCmove{//this is not used since we can change the total effective population size and the proportions vector to achieve the same
    Population *pop;
public:
    NewEffectPopSizeMoveForPopulation(Chain *chain, std::string nameMove, Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel, vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
};

class NewTimeOriginOnTreeforPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewTimeOriginOnTreeforPopulationMove(Chain *chain, std::string nameMove, Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void increaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
class NewTimeOriginOnEdgeforPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewTimeOriginOnEdgeforPopulationMove(Chain *chain, std::string nameMove, Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel, vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions,gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
};
class NewGlobalScaledMutationRateMove:public MCMCmove{
    Population *pop;
public:
    NewGlobalScaledMutationRateMove(Chain *chain, std::string nameMove, Population *pop, vector<MCMCParameterWithKernel *> &mcmcParKernel,  vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void increaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
class NewGlobalScaledGrowthRateForPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewGlobalScaledGrowthRateForPopulationMove(Chain *chain, std::string nameMove, Population *pop,vector<MCMCParameterWithKernel *> &mcmcParKernel, vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void rollbackMove(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions, gsl_rng *randomGenerator);
    void increaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
#endif /* mcmc_move_hpp */

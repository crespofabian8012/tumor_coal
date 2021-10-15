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


//#include "mcmc_chain.hpp"
#include "mcmc_parameter.hpp"


class MCMCParameterWithKernel;
class MCMCVectorParameter;
class Chain;
class MCMCmove{
    int chainId;
    //Chain *chain;
    std::string nameMove;
protected:
    int numberAccept;
    int numberReject;
    int numberAttemps;
    std::vector<MCMCParameterWithKernel *> mcmcParamKernels;
    std::vector<IMCMCParameter<long double> *> affectedParameters;
    
    //MCMCVectorParameter * vectorParam;
    MCMCVectorParameterWithKernel *vectorParam;
    
    MCMCParameter<std::vector<long double>> *affectedVectorParameters;
    //MCMCParameter<vector<long double>> sampleSizesVector;
    
    long double newLogConditionalLikelihoodTree;
    long double newLogConditionalLikelihoodSequences;
    long double chainLogConditionalLikelihoodTree;
    long double chainLogConditionalLikelihoodSequences;
public:
    bool isInvalidMove;
    std::string name();
    MCMCmove(int chainId, std::string nameMove, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam,  long double currentlogConditionalLikelihoodTree,
                      long double currentlogConditionalLikelihoodSequences);
    virtual void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator, Chain* chain)=0;
    virtual void rollbackMove( ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)=0;
    virtual void saveCurrentValue( ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)=0;
    virtual long  double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, Chain* chain)=0;
    virtual void increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)=0;
    virtual void decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions)=0;
    long double getKernelValue(){return mcmcParamKernels[0]->getKernelParameter();}
    void setKernelValue(long double newKernelParameterValue){ mcmcParamKernels[0]->modifyKernelParameter(newKernelParameterValue);}
    void  setCurrentLogLikelihoods(long double currentlogConditionalLikelihoodTree,
                                             long double currentlogConditionalLikelihoodSequences);
    long double  getCurrentLogLikelihoodSequences() const;
    long double  getCurrentLogLikelihoodTree() const;
    double numberAccepted();
    void  resetNumberAccepted();
    void  resetNumberRejected();
    double numberRejected();
   // Chain * getChain();
    int  getChainId();
public:
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain);
    
    virtual ~MCMCmove() = default; // make dtor virtual
    MCMCmove(MCMCmove&&) = default;  // support moving
    MCMCmove& operator=(MCMCmove&&) = default;
    MCMCmove(const MCMCmove&) = default; // support copying
    MCMCmove& operator=(const MCMCmove&) = default;
    
};

class NewProportionsVectorMove:public MCMCmove{
public:
    NewProportionsVectorMove(int chainId, std::string nameMove, std::vector<MCMCParameterWithKernel *> &mcmcParKernel, std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences);
    void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator,Chain* chain);
    void rollbackMove(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost,Chain* chain);
    void increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
};
class NewGrowthRateMoveForPopulation:public MCMCmove{
    Population *pop;
public:
    NewGrowthRateMoveForPopulation(int chainId, std::string nameMove,Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences);
    void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator,Chain* chain);
    void rollbackMove(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl, boost::mt19937* rngBoost, Chain* chain);
};

class NewTimeOriginOnTreeforPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewTimeOriginOnTreeforPopulationMove(int chainId, std::string nameMove, Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences);
    void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const  gsl_rng *randomGenerator,Chain* chain);
    void rollbackMove(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain);
    void increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
};
class NewTimeOriginOnEdgeforPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewTimeOriginOnEdgeforPopulationMove(int chainId, std::string nameMove, Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel, std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences);
    void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator,Chain* chain);
    void rollbackMove(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,const  gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain);
};
class NewGlobalScaledMutationRateMove:public MCMCmove{
    Population *pop;
public:
    NewGlobalScaledMutationRateMove(int chainId, std::string nameMove, Population *pop, std::vector<MCMCParameterWithKernel *> &mcmcParKernel,  std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences);
    void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator,Chain* chain);
    void rollbackMove(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain);
    void increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
};
class NewGlobalScaledGrowthRateForPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewGlobalScaledGrowthRateForPopulationMove(int chainId, std::string nameMove, Population *pop,std::vector<MCMCParameterWithKernel *> &mcmcParKernel, std::vector<IMCMCParameter<long double> *> &affectedParameters, MCMCVectorParameterWithKernel &vectorParam, long double currentlogConditionalLikelihoodTree,
    long double currentlogConditionalLikelihoodSequences);
    void makeProposal(ProgramOptions &programOptions, MCMCOptions &mcmcOptions, const gsl_rng *randomGenerator,Chain* chain);
    void rollbackMove(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    long double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void saveCurrentValue(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,Chain* chain);
    void  move(ProgramOptions &programOptions, MCMCOptions &mcmcOptions,const gsl_rng *rngGsl,boost::mt19937* rngBoost, Chain* chain);
    void increaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
    void decreaseParameter(ProgramOptions &programOptions, MCMCOptions &mcmcOptions);
};
#endif /* mcmc_move_hpp */

//
//  mcmc_move.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 28/10/19.
//

#ifndef mcmc_move_hpp
#define mcmc_move_hpp

#include "mcmc_chain.hpp"

class MCMCmove{
    Chain *chain;
    string nameMove;
protected:
    int numberAccept;
    int numberReject;
    int numberAttemps;
    double newLogConditionalLikelihoodTree;
    double newLogConditionalLikelihoodSequences;
public:
    MCMCmove(Chain *chain, string nameMove);
    virtual void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    virtual void rollbackMove()=0;
    virtual void safeCurrentValue()=0;
    virtual double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    double numberAccepted();
    double numberRejected();
    Chain * getChain();
public:
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
   
   // virtual ~MCMCmove()=0;
};

class NewTotalEffectPopSizeMove:public MCMCmove{
public:
    NewTotalEffectPopSizeMove(Chain *chain, string nameMove);
     void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
     void rollbackMove();
     double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void safeCurrentValue();
     void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};

class NewProportionsVectorMove:public MCMCmove{
public:
    NewProportionsVectorMove(Chain *chain, string nameMove);
     void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
     void rollbackMove();
     double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
     void safeCurrentValue();
     void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
class NewGrowthRateMoveForPopulation:public MCMCmove{
    Population *pop;
public:
    NewGrowthRateMoveForPopulation(Chain *chain, string nameMove, Population *pop);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void rollbackMove();
    double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void safeCurrentValue();
     void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};

class NewEffectPopSizeMoveForPopulation:public MCMCmove{//this is not used since we can change the total effective population size and the proportions vector to achieve the same
    Population *pop;
public:
    NewEffectPopSizeMoveForPopulation(Chain *chain, string nameMove, Population *pop);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void rollbackMove();
    double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void safeCurrentValue();
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};

class NewTimeOriginOnTreeforPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewTimeOriginOnTreeforPopulationMove(Chain *chain, string nameMove, Population *pop);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void rollbackMove();
    double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void safeCurrentValue();
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
class NewTimeOriginOnEdgeforPopulationMove:public MCMCmove{
    Population *pop;
public:
    NewTimeOriginOnEdgeforPopulationMove(Chain *chain, string nameMove, Population *pop);
    void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void rollbackMove();
    double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void safeCurrentValue();
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
};
#endif /* mcmc_move_hpp */

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
public:
    MCMCmove(Chain *chain);
    virtual void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    virtual void rollbackMove()=0;
    virtual void safeCurrentValue()=0;
    virtual double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions)=0;
    void  move(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    Chain * getChain();
   // virtual ~MCMCmove()=0;
};

class NewTotalEffectPopSizeMove:MCMCmove{
public:
     void makeProposal(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
     void rollbackMove();
     double computeLogAcceptanceProb(ProgramOptions &programOptions, MCMCoptions &mcmcOptions);
    void safeCurrentValue();
};

//class NewProportionsVectorMove:MCMCmove{
//public:
//     void makeProposal();
//     void rollbackMove();
//     void computeAcceptanceProb();
//     void safeCurrentValue();
//};
//class NewGrowthRateMove:MCMCmove{
//public:
//    void makeProposal();
//    void rollbackMove();
//    void computeAcceptanceProb();
//    void safeCurrentValue();
//};
#endif /* mcmc_move_hpp */

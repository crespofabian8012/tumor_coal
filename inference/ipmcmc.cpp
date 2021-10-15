//
//  pmcmc.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 29/09/2021.
//

#include "ipmcmc.hpp"


IPMCMC::IPMCMC(size_t numCSMC,size_t numSMC, SMCOptions &smcOptions, PosetSMCParams &posetSMCParams, MCMCOptions &mcmcOptions,PMCMCOptions &pmcmcOptions, PosetSMC &posetSMC ):
smcOptions(smcOptions),pmcmcOptions(pmcmcOptions), mcmcOptions(mcmcOptions),  posetSMCParams(posetSMCParams), numCSMC(numCSMC), numSMC(numSMC)
{
    
    for(size_t i =0; i< numSMC; ++i){
        SMCs.emplace_back(new SMC<State, PosetSMCParams>(posetSMC, smcOptions) );
    }
    for(size_t i =0; i< numCSMC; ++i){
        
        condSMCs.emplace_back(new ConditionalSMC<State, PosetSMCParams>(posetSMC, smcOptions) );
    }
    
}
//void initialize(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions,               MCMCoptions &mcmcOptions, FilePaths &filePaths,
//                   std::vector<std::vector<int> > &ObservedData,
//                   char* ObservedCellNames[],
//                   MSA &msa,
//                   std::string& healthyTipLabel
//                  );

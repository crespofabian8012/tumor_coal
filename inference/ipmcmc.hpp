//
//  pmcmc.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 29/09/2021.
//

#ifndef pmcmc_hpp
#define pmcmc_hpp


#include "smc.hpp"
#include "csmc.hpp"
#include "poset_smc.hpp"
#include "pmcmc_options.hpp"
#include "mcmc_parameter.hpp"
#include "chain_manager.hpp"
#include "smc_options.hpp"

#include "poset_smc_params.hpp"
#include "random.h"
#include <vector>

class State;
class PosetSMCParams;

using CondSMC = ConditionalSMC<State,PosetSMCParams>;
using UncondSMC = SMC<State,PosetSMCParams>;
using namespace std;
class IPMCMC {
    
    size_t numMCMCiterations = 0;
    SMCOptions smcOptions;
    PMCMCOptions pmcmcOptions;
    MCMCOptions mcmcOptions;
    PosetSMCParams posetSMCParams;
   // std::vector<MCMCParameterWithKernel> mcmcParams;//MCMCParameterWithKernel includes the proposal distribution
    //ModelLikelihood modelLik;
    std::vector<std::shared_ptr<CondSMC>> condSMCs;
    std::vector<std::shared_ptr<UncondSMC>> SMCs;
    size_t numCSMC = 4;
    size_t numSMC = 4;
    
public:
    IPMCMC(size_t numCSMC,size_t numSMC, SMCOptions &smcOptions, PosetSMCParams &posetSMCParams, MCMCOptions &mcmcOptions,PMCMCOptions & pmcmcOptions, PosetSMC &posetSMC );
    
    void initialize(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions,               MCMCOptions &mcmcOptions, FilePaths &filePaths,
                    std::vector<std::vector<int> > &ObservedData,
                    char* ObservedCellNames[],
                    MSA &msa,
                    std::string& healthyTipLabel
                   );
    
    void runIPMCMC(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCOptions &mcmcOptions, FilePaths &filePaths);
    
   
    
    ~IPMCMC(){
        
    };
};

#endif /* pmcmc_hpp */

//
//  mcmc.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 15/03/2021.
//

#ifndef mcmc_hpp
#define mcmc_hpp

#include "data_types.hpp"
#include "chain_manager.hpp"

extern "C"
    {
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
#include <libpll/pll.h>
    
    }

class MCMC {
    
    ChainManager * chainManager;
    
public:
    MCMC(int numChains);
    
    void initialize(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions,               MCMCoptions &mcmcOptions, FilePaths &filePaths,
                    std::vector<int> &sampleSizes,
                    std::vector<std::vector<int> > &ObservedData,
                    char* ObservedCellNames[],
                    pll_msa_t *msa,
                    pll_rtree_t * initialRootedTree,
                    std::vector<StructuredCoalescentTree *> structuredCoalTrees,
                    std::string& healthyTipLabel,
                    const std::vector<pll_rtree_t *> &trueTrees,        const  std::vector<long double> &trueThetas,
                    const std::vector<std::vector<long double>> &trueDeltaTs,
                    const  std::vector<std::vector<long double>> &trueTs, std::vector<Partition *> &partitions);
    
    void runMCMC(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths);
    
    
    ~MCMC();
};

#endif /* mcmc_hpp */

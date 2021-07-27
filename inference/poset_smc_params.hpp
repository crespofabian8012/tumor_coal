//
//  poset_smc_params.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 17/03/2021.
//

#ifndef poset_smc_params_hpp
#define poset_smc_params_hpp

#include  "data_types.hpp"
#include  "mcmc_parameter.hpp"
#include  "pll_buffer_manager.hpp"
#include  "genotype_error_model.hpp"

#include  <vector>


extern "C"
{
#include "libpll/pll.h"
}

class PosetSMCParams{
       
    
public:
    int numberClones;
    int sampleSize;
    unsigned int num_sites;
    pll_msa_t *msa;
    const pll_partition_t *partition;
    PLLBufferManager *const pll_buffer_manager;
    std::vector<int> positions;
    ProgramOptions *programOptions;
    GenotypeErrorModel *gtErrorModel;
    
    std::shared_ptr<MCMCParameterWithKernel> theta;
    
    std::shared_ptr<MCMCParameterWithKernel> seqError;
    std::shared_ptr<MCMCParameterWithKernel> dropoutError;
    std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationDeltaTs;
    std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationToriginSTDs;
    std::shared_ptr<MCMCVectorParameterWithKernel>proportions;

    PosetSMCParams(int numberClones,
                   int sampleSize,
                   unsigned int num_sites,
                   pll_msa_t *msa,
                   const pll_partition_t *partition,
                   PLLBufferManager *const pll_buffer_manager,
                   std::vector<int> &positions,
                   ProgramOptions &programOptions,
                   GenotypeErrorModel *gtErrorModel);
    
    void set(std::shared_ptr<MCMCParameterWithKernel> thetaPar,
             std::shared_ptr<MCMCParameterWithKernel> seqErrorPar,
             std::shared_ptr<MCMCParameterWithKernel> dropoutErrorPar,
             std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationDeltaTPar,
             std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationToriginSTDPar,
             std::shared_ptr<MCMCVectorParameterWithKernel>proportionsPar
             );
    ProgramOptions& getProgramOptions();
    std::shared_ptr<MCMCParameterWithKernel> getTheta();
    std::shared_ptr<MCMCParameterWithKernel> getSeqError();
    std::shared_ptr<MCMCParameterWithKernel> getDropoutError() ;
    std::shared_ptr<MCMCParameterWithKernel> getPopulationDeltaT(int i);
    std::shared_ptr<MCMCParameterWithKernel> getPopulationToriginSTD(int i);
    std::shared_ptr<MCMCVectorParameterWithKernel> getProportionVector() ;
    int  getSampleSize() const;
    int  getNumClones() const;
    
};

#endif /* poset_smc_params_hpp */

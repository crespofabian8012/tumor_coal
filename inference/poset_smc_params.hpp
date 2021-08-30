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
#include  "partition.hpp"
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
    const Partition *partition;
    PLLBufferManager *const pll_buffer_manager;
    std::vector<int> positions;
    ProgramOptions *programOptions;
    GenotypeErrorModel *gtErrorModel;
    std::vector<int> sampleSizes;
    bool doPriorPost;
    bool doFixedEventimes;
    
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    std::shared_ptr<double> theta;
    std::shared_ptr<double> seqError;
    std::shared_ptr<double> dropoutError;
    std::vector<double> populationDeltaTs;
    std::vector<double> populationToriginSTDs;
    std::vector<double>proportions;
    
    PosetSMCParams(int numberClones,
                   int sampleSize,
                   std::vector<int> &sampleSizes,
                   unsigned int num_sites,
                   pll_msa_t *msa,
                   const Partition *partition,
                   PLLBufferManager *const pll_buffer_manager,
                   std::vector<int> &positions,
                   ProgramOptions &programOptions,
                   GenotypeErrorModel *gtErrorModel,
                   double theta,
                   std::vector<double> populationDeltaTs,
                   std::vector<double> populationToriginSTDs,
                   std::vector<double>proportions,
                   std::vector<std::vector<double>> coalTimesModelTimePerPopulation);
    
    //    void set(std::shared_ptr<MCMCParameterWithKernel> thetaPar,
    //             std::shared_ptr<MCMCParameterWithKernel> seqErrorPar,
    //             std::shared_ptr<MCMCParameterWithKernel> dropoutErrorPar,
    //             std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationDeltaTPar,
    //             std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationToriginSTDPar,
    //             std::shared_ptr<MCMCVectorParameterWithKernel>proportionsPar
    //             );
    ProgramOptions& getProgramOptions();
    //    std::shared_ptr<MCMCParameterWithKernel> getTheta();
    //    std::shared_ptr<MCMCParameterWithKernel> getSeqError();
    //    std::shared_ptr<MCMCParameterWithKernel> getDropoutError() ;
    //    std::shared_ptr<MCMCParameterWithKernel> getPopulationDeltaT(int i);
    //    std::shared_ptr<MCMCParameterWithKernel> getPopulationToriginSTD(int i);
    //    std::shared_ptr<MCMCVectorParameterWithKernel> getProportionVector() ;
    int  getSampleSize() const;
    int  getNumClones() const;
    double getPopulationEvent(int idx_population, int idx_event);
    void buildListEventTimesPerPopulation();
};

#endif /* poset_smc_params_hpp */

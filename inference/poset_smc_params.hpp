//
//  poset_smc_params.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 17/03/2021.
//

#ifndef poset_smc_params_hpp
#define poset_smc_params_hpp

#include  "msa.hpp"
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
    MSA  *msa;
    //pll_msa_t *msa;
    const Partition *partition;
    PLLBufferManager *const pll_buffer_manager;
    std::vector<int> positions;
    ProgramOptions *programOptions;
    GenotypeErrorModel *gtErrorModel;
    std::vector<int> sampleSizes;
    bool doPriorPost;
    bool doFixedEventimes;
    bool usePriorInSMC1 ;
    int verbose;
      
    std::vector<double> trueCoalTimes;
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    double theta;
    double seqError;
    double dropoutError;
    std::vector<double> populationDeltaTs;
    std::vector<double> populationToriginSTDs;
    std::vector<double>proportions;
    
    PosetSMCParams(int numberClones,
                   int sampleSize,
                   std::vector<int> &sampleSizes,
                   unsigned int num_sites,
                   MSA *msa,
                   //pll_msa_t *msa,
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
    
  // PosetSMCParams(const PosetSMCParams &original);
      
  // PosetSMCParams &operator=(const PosetSMCParams &original);
      
   //PosetSMCParams(PosetSMCParams&& rhs);
    

    ProgramOptions& getProgramOptions();
  
    int  getSampleSize() const;
    int  getNumClones() const;
    double getPopulationEvent(int idx_population, int idx_event);
    void buildListEventTimesPerPopulation();
    
    double getNumberPopulations() const{return numberClones;};
    double getTheta() const{return theta;};
    std::vector<double> getScaledGrowthRates() const{return populationDeltaTs;};
    double getScaledGrowthRatesAt(size_t i) const{return populationDeltaTs[i];};
    std::vector<double> getTimeOrigins() const{return populationToriginSTDs;};
    double getTimeOriginAt(size_t i) const{return populationToriginSTDs[i];};
    std::vector<double> getProportionsVector() const{return proportions;};
    double getProportionPopulation(size_t i) const{return proportions[i];};
    ~PosetSMCParams(){};
};

#endif /* poset_smc_params_hpp */

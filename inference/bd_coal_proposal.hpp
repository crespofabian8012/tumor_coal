//
//  bd_coal_proposal.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#ifndef bd_coal_proposal_hpp
#define bd_coal_proposal_hpp


#include  "msa.hpp"
#include "pg_proposal.hpp"
#include "state.hpp"
#include "bd_coal_model_params.hpp"
#include "bd_coal_state.hpp"
#include "poset_smc_params.hpp"
#include "pmmh_proposal.hpp"


class BDCoalModelPGProposal : public PGProposal<State, PosetSMCParams>
{
    size_t numPopulations;
    int  sampleSize;
    std::vector<int>sampleSizes;
    std::vector<int> positions;
    
    MSA *msa;
    const Partition *partition;
    GenotypeErrorModel *gtErrorModel;
    
    PLLBufferManager *const pll_buffer_manager;
    
    double theta;
    std::vector<double> populationDeltaTs;
    std::vector<double> populationToriginSTDs;
    std::vector<double>proportions;
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    BDCoalPriorParams &priorParams;
    ProgramOptions& programOptions;
public:
    shared_ptr<PosetSMCParams> sample_from_prior(gsl_rng *random);
    shared_ptr<PosetSMCParams> propose(gsl_rng *random, const PosetSMCParams &curr, shared_ptr<ParticleGenealogy<State> > genealogy);
    double log_prior(const PosetSMCParams &curr); // log p(curr)
    
    BDCoalModelPGProposal(BDCoalPriorParams &priorParams, size_t numPopulations,int sampleSize,
                                                 unsigned int num_sites,
                                                 MSA *msa,
                                                 const Partition *partition,
                                                 PLLBufferManager *const pll_buffer_manager,
                                                 ProgramOptions &programOptions);
    
    
};

class BDCoalModelRandomWalkProposal : public PMMHProposal<PosetSMCParams>
{
    size_t numPopulations;
    int  sampleSize;
    std::vector<int>sampleSizes;
    std::vector<int> positions;
    
    MSA *msa;
    const Partition *partition;
    GenotypeErrorModel *gtErrorModel;
    
    PLLBufferManager *const pll_buffer_manager;
    
    double theta;
    std::vector<double> populationDeltaTs;
    std::vector<double> populationToriginSTDs;
    std::vector<double>proportions;
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    BDCoalPriorParams &priorParams;
    ProgramOptions& programOptions;
public:
    PosetSMCParams* sample_from_prior(gsl_rng *random);
    PosetSMCParams* propose(gsl_rng *random,  PosetSMCParams *curr);
    double log_proposal(PosetSMCParams *curr, PosetSMCParams *prev);
    double log_prior(PosetSMCParams *curr); // log p(curr)
    void adapt(size_t num_accepts, size_t curr_iter);
    
    BDCoalModelRandomWalkProposal(BDCoalPriorParams &priorParams, size_t numPopulations,int sampleSize,
                                                 unsigned int num_sites,
                                                 MSA *msa,
                                                 const Partition *partition,
                                                 PLLBufferManager *const pll_buffer_manager,
                                                 ProgramOptions &programOptions);
    
    
};

#endif /* bd_coal_proposal_hpp */

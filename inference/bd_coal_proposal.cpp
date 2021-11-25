//
//  bd_coal_proposal.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#include "bd_coal_proposal.hpp"
#include "poset_smc.hpp"


BDCoalModelPGProposal::BDCoalModelPGProposal(BDCoalPriorParams &priorParams,                                             size_t numPopulations,int sampleSize,
                                             unsigned int num_sites,
                                             MSA *msa,
                                             const Partition *partition,
                                             PLLBufferManager *const pll_buffer_manager,
                                             ProgramOptions &programOptions):
numPopulations(numPopulations), sampleSize(sampleSize), msa(msa), partition(partition), 
pll_buffer_manager(pll_buffer_manager), priorParams(priorParams), programOptions(programOptions)
{
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        populationDeltaTs.push_back(0.0);
        populationToriginSTDs.push_back(0.0);
        proportions.push_back(0.0);
    }
    
}
shared_ptr<PosetSMCParams> BDCoalModelPGProposal::sample_from_prior(gsl_rng *random){
    
    sampleSizes = {sampleSize}; //TODO for p populations divide
    //sampleSize in p sums
    gtErrorModel = new GenotypeErrorModel("GT20", 0.0,  0.0, 16);
    theta = Random::RandomExponential(priorParams.lambdaExpPriorTheta, NULL, true, random, NULL);
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        populationDeltaTs[i] = Random::RandomExponential(priorParams.lambdaExpPriorScaledGrowthRate, NULL, true, random, NULL);
        
        populationToriginSTDs[i] = Random::RandomDensityModelTimeOrigin(populationDeltaTs[i], true, 0.0, random, NULL);
        
        proportions[i] = sampleSizes[i] / (1.0*sampleSize);
    }
    
    std::vector<int> positions(sampleSize);
    std::iota( std::begin( positions ), std::end( positions ), 1 );
    coalTimesModelTimePerPopulation = {};
    
    
    shared_ptr<PosetSMCParams> result = std::make_shared<PosetSMCParams>(numPopulations, sampleSize,  sampleSizes, msa->getLength(), msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                                                                         theta,
                                                                         populationDeltaTs,
                                                                         populationToriginSTDs,
                                                                         proportions,
                                                                         coalTimesModelTimePerPopulation);
    
    
    return result;
};
shared_ptr<PosetSMCParams> BDCoalModelPGProposal::propose(gsl_rng *random, const PosetSMCParams &curr, shared_ptr<ParticleGenealogy<State> > genealogy){
    
    // size_t T = genealogy->size();
    long double randomNumber= Random::randomUniformFromGsl2(random);
    long double m = exp(2 * log(priorParams.paramMultiplierProposalTheta) *(randomNumber -0.5));
    double theta = curr.getTheta() * m;
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        randomNumber= Random::randomUniformFromGsl2(random);
        m = exp(2 * log(priorParams.paramMultiplierProposalScaledGrowthRate) *(randomNumber -0.5));
        
        populationDeltaTs[i] = curr.getScaledGrowthRatesAt(i) *m;
        populationToriginSTDs[i] = Random::RandomDensityModelTimeOrigin(populationDeltaTs[i], true, 0.0, random, NULL);
        
        proportions[i] = sampleSizes[i] / (1.0*sampleSize);
    }
    
    
    shared_ptr<PosetSMCParams> result = std::make_shared<PosetSMCParams>(numPopulations, sampleSize,  sampleSizes, msa->getLength(), msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                                                                         theta,
                                                                         populationDeltaTs,
                                                                         populationToriginSTDs,
                                                                         proportions,
                                                                         coalTimesModelTimePerPopulation);
    
    result->doPriorPost = true;
    result->doFixedEventimes = false;
    result->usePriorInSMC1 = true;
    bool normalizedCLVs = false;
    
    //         num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    //         posetSMC = new PosetSMC(programOptions.numClones,  num_iter, doPlots);
    
    PosetSMC::PosetSMCKernel kernelType = PosetSMC::PosetSMCKernel::TSMC1;
    
    if (kernelType== PosetSMC::PosetSMCKernel::TSMC1){
        
        if (!(result->usePriorInSMC1)){
            
            if (normalizedCLVs){
                programOptions.normalizeLeavesClv =true;
                programOptions.normalizeClv =true;
            }
            else{
                programOptions.normalizeLeavesClv =false;
                programOptions.normalizeClv =false;
            }
        }
        else{
            programOptions.normalizeLeavesClv =false;
            programOptions.normalizeClv =false;
        }
    }
    
    
    return result;
}

double BDCoalModelPGProposal::log_prior(const PosetSMCParams &curr){
    
    double log_prior = 0.0;
    
    
    return log_prior;
} // log p(curr)

///
BDCoalModelRandomWalkProposal::BDCoalModelRandomWalkProposal(BDCoalPriorParams &priorParams, size_t numPopulations,int sampleSize,
                                                             unsigned int num_sites,
                                                             MSA *msa,
                                                             const Partition *partition,
                                                             PLLBufferManager *const pll_buffer_manager,
                                                             ProgramOptions &programOptions):
numPopulations(numPopulations), sampleSize(sampleSize), msa(msa), partition(partition),
pll_buffer_manager(pll_buffer_manager), priorParams(priorParams), programOptions(programOptions)
{
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        populationDeltaTs.push_back(0.0);
        populationToriginSTDs.push_back(0.0);
        proportions.push_back(0.0);
    }
    
}

PosetSMCParams* BDCoalModelRandomWalkProposal::sample_from_prior(gsl_rng *random){
    
    sampleSizes = {sampleSize}; //TODO for p populations divide
    //sampleSize in p sums
    gtErrorModel = new GenotypeErrorModel("GT20", 0.0,  0.0, 16);
    theta = Random::RandomExponential(priorParams.lambdaExpPriorTheta, NULL, true, random, NULL);
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        populationDeltaTs[i] = Random::RandomExponential(priorParams.lambdaExpPriorScaledGrowthRate, NULL, true, random, NULL);
        
        populationToriginSTDs[i] = Random::RandomDensityModelTimeOrigin(populationDeltaTs[i], true, 0.0, random, NULL);
        
        proportions[i] = sampleSizes[i] / (1.0*sampleSize);
    }
    
    std::vector<int> positions(sampleSize);
    std::iota( std::begin( positions ), std::end( positions ), 1 );
    coalTimesModelTimePerPopulation = {};
    
    
    PosetSMCParams* result = new PosetSMCParams(numPopulations, sampleSize,  sampleSizes, msa->getLength(), msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                                                theta,
                                                populationDeltaTs,
                                                populationToriginSTDs,
                                                proportions,
                                                coalTimesModelTimePerPopulation);
    
    result->doPriorPost = true;
    result->doFixedEventimes = false;
    result->usePriorInSMC1 = true;
    bool normalizedCLVs = false;
    
    //         num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    //         posetSMC = new PosetSMC(programOptions.numClones,  num_iter, doPlots);
    
    PosetSMC::PosetSMCKernel kernelType = PosetSMC::PosetSMCKernel::TSMC1;
    
    if (kernelType== PosetSMC::PosetSMCKernel::TSMC1){
        
        if (!(result->usePriorInSMC1)){
            
            if (normalizedCLVs){
                programOptions.normalizeLeavesClv =true;
                programOptions.normalizeClv =true;
            }
            else{
                programOptions.normalizeLeavesClv =false;
                programOptions.normalizeClv =false;
            }
        }
        else{
            programOptions.normalizeLeavesClv =false;
            programOptions.normalizeClv =false;
        }
    }
    
    
    return result;
}
PosetSMCParams* BDCoalModelRandomWalkProposal::propose(gsl_rng *random,  PosetSMCParams *curr){
    
    long double randomNumber= Random::randomUniformFromGsl2(random);
    long double m = exp(2 * log(priorParams.paramMultiplierProposalTheta) *(randomNumber -0.5));
    double theta = curr->getTheta() * m;
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        randomNumber= Random::randomUniformFromGsl2(random);
        m = exp(2 * log(priorParams.paramMultiplierProposalScaledGrowthRate) *(randomNumber -0.5));
        
        populationDeltaTs[i] = curr->getScaledGrowthRatesAt(i) *m;
        populationToriginSTDs[i] = Random::RandomDensityModelTimeOrigin(populationDeltaTs[i], true, 0.0, random, NULL);
        
        proportions[i] = sampleSizes[i] / (1.0*sampleSize);
    }
    
    PosetSMCParams* result =  new PosetSMCParams(numPopulations, sampleSize,  sampleSizes, msa->getLength(), msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                                                 theta,
                                                 populationDeltaTs,
                                                 populationToriginSTDs,
                                                 proportions,
                                                 coalTimesModelTimePerPopulation);
    
    result->doPriorPost = true;
    result->doFixedEventimes = false;
    result->usePriorInSMC1 = true;
    bool normalizedCLVs = false;
    
    //         num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    //         posetSMC = new PosetSMC(programOptions.numClones,  num_iter, doPlots);
    
    PosetSMC::PosetSMCKernel kernelType = PosetSMC::PosetSMCKernel::TSMC1;
    
    if (kernelType== PosetSMC::PosetSMCKernel::TSMC1){
        
        if (!(result->usePriorInSMC1)){
            
            if (normalizedCLVs){
                programOptions.normalizeLeavesClv =true;
                programOptions.normalizeClv =true;
            }
            else{
                programOptions.normalizeLeavesClv =false;
                programOptions.normalizeClv =false;
            }
        }
        else{
            programOptions.normalizeLeavesClv =false;
            programOptions.normalizeClv =false;
        }
    }
    
    
    return result;
}
double BDCoalModelRandomWalkProposal::log_proposal(PosetSMCParams *curr, PosetSMCParams *prev){
    
    double log_proposal = 0.0;
    log_proposal+= log(curr->getTheta() / prev->getTheta());
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        log_proposal+= log(curr->getScaledGrowthRatesAt(i) / prev->getScaledGrowthRatesAt(i));
    }
    return log_proposal;
}
double BDCoalModelRandomWalkProposal::log_prior(PosetSMCParams *curr){
    
    double log_prior = 0.0;
    
    long double currTheta = curr->getTheta();
    long double zero = 0.0;
    long double lambdaPriorTheta = priorParams.lambdaExpPriorTheta;
    log_prior += Distributions::LogExponentialDensity(currTheta, zero, lambdaPriorTheta);
    
    for(size_t i = 0 ; i < numPopulations; ++i)
    {
        long double scaledGrowthRate = curr->getScaledGrowthRatesAt(i);
        long double lambdaPriorDelta = priorParams.lambdaExpPriorScaledGrowthRate;
        log_prior += Distributions::LogExponentialDensity(scaledGrowthRate, zero, lambdaPriorDelta);
    }
    
    return log_prior;
} // log p(curr)
void BDCoalModelRandomWalkProposal::adapt(size_t num_accepts, size_t curr_iter){
    
    
    
}

/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/

/*
 * MCMC inference of scaled growth rates
 */




#include <iostream>
#include <string>

#include <map>
#include <unordered_map>
#include <string>
#include <omp.h>

#include "random.h"
#include "data_utils.hpp"
#include "population.hpp"
#include "mcmc_chain.hpp"
#include "chain_manager.hpp"
#include "definitions.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "output_functions.hpp"
#include "mcmc.hpp"

#include <gsl/gsl_randist.h>
#include <boost/program_options.hpp>


#include "smc.hpp"
#include "state.hpp"
#include "poset_smc.hpp"




extern "C"
    {
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
#include <libpll/pll.h>
    }
#include <stdarg.h>
#include <search.h>
#include <time.h>



int main(int argc, char* argv[] )
{
    FILE *input_file;
    const char* input_path;
    std::string config_file;
    
    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("config_file,c", po::value<std::string>(&config_file)->required(), "path to configuration file.")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    input_path = config_file.c_str();
    
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    bool doMCMC = false;
    
    if (config_file.empty())
    {
        std::cout << "\nERROR: No parameters specified (use command  parameter file)"<< std::endl;
        Output::PrintUsage();
    }
    else
    {
        input_path = config_file.c_str();
        // 1. call function to parse the MCMC configuration file
        if ((input_file = freopen(input_path, "r", stdin)) != NULL)
        {
            ReadMCMCParametersFromFile(programOptions, filePaths, mcmcOptions);
        }
        else
        {
            if (!mcmcOptions.noData)
            {
                std::cout << "\nERROR: No parameters specified (use command line or parameter file)"<< std::endl;
                exit(-1);
            }
        }
    }
    printProgramHeader();
    setDefaultOptions(programOptions, mcmcOptions);
    
    std::vector<StructuredCoalescentTree *> structuredCoalTrees;
    std::vector<pll_rtree_t *> trueTrees;
    std::vector<long double> trueThetas;
    std::vector<std::vector<long double>> trueDeltaTs;
    std::vector<std::vector<long double>> trueTs;
    
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip ;//= argv[3];
    char *treefileName; //= argv[4];
    
    std::vector<int> sampleSizes;
    
    //2. initialize data structures
    //    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    //
    //    //3. do inference
    
    pll_msa_t *msa;
    
    //    if (mcmcOptions.useSequencesLikelihood ==1)
    //   {
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    Utils::ReadParametersFromFastaFile(fileNameFasta,  programOptions.numCells, programOptions.TotalNumSequences, programOptions.numSites);
    std::vector<std::vector<int> > ObservedData;
    
    char *ObservedCellNames[programOptions.numCells];
    Utils::ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    programOptions.numCells = programOptions.TotalNumSequences;
    
    programOptions.TotalTumorSequences=programOptions.TotalNumSequences-1;
    
    if (programOptions.numClones==1)
        sampleSizes.push_back(programOptions.numCells-1); //minus the healthycell
    
    fileNamePhylip =filePaths.inputGenotypeFilePhylip;
    msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
    if (!msa)
        std::cout << "Error reading phylip file \n"<< std::endl;
    // }
    unsigned long stats_mask = PLLMOD_MSA_STATS_DUP_TAXA;
    pllmod_msa_stats_t * stats = pllmod_msa_compute_stats(msa,
                                                          10,
                                                          pll_map_gt10, // map is not used here
                                                          NULL,
                                                          stats_mask);
    assert(stats);
    if (stats->dup_taxa_pairs_count > 0){
        
        std::cout << "Error, Duplicate sequence names found: \n"<< stats->dup_taxa_pairs_count <<  std::endl;
    }
    stats_mask = PLLMOD_MSA_STATS_DUP_SEQS;
    vector<pair<unsigned long,unsigned long> > dup_seqs;
    std::set<size_t> gap_seqs;
    pllmod_msa_destroy_stats(stats);
    pllmod_msa_stats_t * dup_stats = pllmod_msa_compute_stats(msa,
                                                              10,
                                                              pll_map_gt10,
                                                              NULL,
                                                              stats_mask);
    
    
    assert(dup_stats);
    
    for (unsigned long c = 0; c < dup_stats->dup_seqs_pairs_count; ++c)
    {
        std::cout << dup_stats->dup_seqs_pairs[c]<< std::endl;
        //dup_seqs.emplace_back(dup_stats->dup_seqs_pairs[c*2]);
        dup_seqs.emplace_back(dup_stats->dup_seqs_pairs[c*2],
                              dup_stats->dup_seqs_pairs[c*2+1]);
        
        std::cout << "Error, Duplicate sequences found: \n" <<  std::endl;
        
    }
    pllmod_msa_destroy_stats(dup_stats);
    
    
    
    treefileName = filePaths.inputTreeFile;
    if ((treefileName != NULL) && (treefileName[0] == '\0'))
    {
        std::cout << "\nERROR: No tree file  specified (use command line or parameter file)" << std::endl;
        exit(-1);
    }
    std::string healthyTipLabel = "healthycell";
    programOptions.healthyTipLabel ="healthycell";
    
    //int sampleEvery = mcmcOptions.thinning;
    
    std::vector<gsl_rng * > randomGenerators(mcmcOptions.numChains);
    
    // int maxNumberThreads = omp_get_max_threads();
    Random::allocateListRandomNumbersGenerators(randomGenerators);
    std::vector<boost::mt19937 *> randomGeneratorsBoost = Random::allocateListRandomNumbersGeneratorsBoost(mcmcOptions.numChains);
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    
    programOptions.seqErrorRate = 0.01;
    programOptions.dropoutRate =0.2;
    mcmcOptions.noData = false;
    mcmcOptions.useGSLRandomGenerator = true;
    mcmcOptions.splitThetaDeltaTmoves=true;
    programOptions.doUseFixedTree=NO;
    //int iterationToComputeThinnig= floor(mcmcOptions.percentIterationsToComputeThinnig * mcmcOptions.Niterations);
    mcmcOptions.doMCMCMoveTimeOriginInputOldestPop=false;
    mcmcOptions.iterationToComputeThinnig= floor(0.25 *mcmcOptions.Niterations);
    mcmcOptions.numberWarmUpIterations = 10000;//floor(0.20 *mcmcOptions.Niterations);
    mcmcOptions.iterationsToMonitorChain= floor(0.05*mcmcOptions.iterationToComputeThinnig);
    mcmcOptions.iterationsToMonitorChain =   floor(mcmcOptions.numberWarmUpIterations / 20);
    mcmcOptions.verbose=0;
    mcmcOptions.burnInIterations= 5000;
    mcmcOptions.Niterations= 500000;
    mcmcOptions.maxNumberIndependentPosteriorValues=200;
    mcmcOptions.maxNumberIndependentPosteriorValues=500;
    mcmcOptions.numberIterationsAfterConvergence = 10000;
    
    
    programOptions.doUsefixedMutationRate = false;
    programOptions.K=0.8;
    
    
    
    if (doMCMC){
        std::vector<Chain*> chainsFortheSameTree(mcmcOptions.numberChainsPerTree);
        int numberTrees = floor((double)mcmcOptions.numChains / mcmcOptions.numberChainsPerTree);
        
        mcmcOptions.fixedValuesForSimulation = true;
        programOptions.doUsefixedMutationRate=true;
        /*simulate trees to initialize the chains*/
        if (programOptions.doUseFixedTree == NO){
            std::cout << "\nSimulating " << numberTrees << " trees .... \n"<<std::endl;
            simulateTrees( numberTrees,structuredCoalTrees,  trueTrees,         trueThetas,
                          trueDeltaTs,
                          trueTs,
                          sampleSizes, programOptions,mcmcOptions, randomGenerators,randomGeneratorsBoost, ObservedData,ObservedCellNames, msa,  healthyTipLabel,   programOptions.seqErrorRate,
                          programOptions.dropoutRate  );
            std::cout << "\nSimulation trees ended. \n"<< std::endl;
            
        }
        mcmcOptions.printChainStateEvery=mcmcOptions.iterationsToMonitorChain;
        std::vector<Partition *> partitionList(mcmcOptions.numChains);
        
        MCMC * mcmc =new MCMC(mcmcOptions.numChains);
        
        mcmc->initialize(randomGenerators, programOptions, mcmcOptions, filePaths,
                         sampleSizes,
                         ObservedData,
                         ObservedCellNames,
                         msa,
                         initialRootedTree,
                         structuredCoalTrees,
                         healthyTipLabel,
                         trueTrees,    trueThetas,
                         trueDeltaTs,
                         trueTs, partitionList);
        
        float start = clock();
        clock_t begin = omp_get_wtime();
        //bool converged =false;
        
        mcmc->runMCMC( randomGenerators, programOptions, mcmcOptions, filePaths);
        
        clock_t end = omp_get_wtime();
        if (!mcmcOptions.noData && mcmcOptions.useSequencesLikelihood ==1)
            pll_msa_destroy(msa);
        
        
        double elapsed_time = double(end - begin);
        std::cout << "\n\n*** Program finished ***" <<std::endl;
        /* execution time */
        double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
        std::cout << "\n\n_________________________________________________________________" <<std::endl;
        std::cout << "\n sSerial Time processing: " << secs << " seconds and parallel Time processing " << elapsed_time  <<std::endl;
        
    }
    else/* SMC */
    {
        
        programOptions.meanADOsite = 0.1;
        programOptions.varADOsite=0.01;
        programOptions.meanADOcell = 0.1;
        programOptions.varADOcell=0.01;
        programOptions.meanGenotypingError= 0.01;
        programOptions.varGenotypingError=0.001;
        programOptions.fixedADOrate=0.01;
        
        
        std::vector<int> positions(programOptions.TotalTumorSequences);
        std::iota( std::begin( positions ), std::end( positions ), 1 );
        //std::random_shuffle( positions.begin(), positions.end());
        
        
        GenotypeErrorModel *gtErrorModel= new GenotypeErrorModel("GT20", programOptions.meanGenotypingError,  1.0 - sqrt (1.0 - programOptions.fixedADOrate), 16);
        PLLBufferManager *pll_buffer_manager = new PLLBufferManager;
        //        const pll_partition_t* pll_partition= pll_utils::createGTReferencePartition(msa);
        
        Partition * partition = new Partition(msa,
                                              16,// model->states,//numberStates
                                              1,//RATE_CATS, // unsigned  int  numberRateCats
                                              0, //int statesPadded
                                              true,//PLL_ATTRIB_ARCH_SSE
                                              false, false, false, false, false);
        
        
        PosetSMCParams psParams(programOptions.numClones, programOptions.TotalNumSequences,  sampleSizes,programOptions.numSites, msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel);
        
        size_t num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
        PosetSMC posetSMC(programOptions.numClones,  num_iter);
        SMCOptions smcOptions;
        
        smcOptions.num_threads = 1;
        smcOptions.use_SPF = false;
        smcOptions.ess_threshold = 1;
        smcOptions.num_particles = 7000;
        smcOptions.resample_last_round = true;
        //smcOptions.resampling_scheme= 0;
        
        
        smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::MULTINOMIAL;
        //smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::STRATIFIED;
        //smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::SYSTEMATIC;
        // smcOptions.main_seed = 7482984;
        //smcOptions.resampling_seed = 849284;
        smcOptions.track_population = true;
        smcOptions.init();
        smcOptions.debug = true;
        
        SMC<State, PosetSMCParams>  smc(posetSMC, smcOptions);
        
        std::cout<< "Running Sequential Monte Carlo(SMC) with "<<smcOptions.num_particles <<" particles....\n" << std::endl;
        smc.run_smc(psParams);
        
        ParticlePopulation<State> *currenPop = smc.get_curr_population();
        vector<shared_ptr<State>> *particles = currenPop->get_particles();
        double   log_marginal = smc.get_log_marginal_likelihood();
        ParticlePopulation<State> *pop0 = smc.get_population(0);
        
        
        int num_particles = pop0[0].get_num_particles();
        
        std::vector<long double> deltas0;
        std::vector<long double> Ts0;
        std::vector<long double> Thetas;
        std::vector<long double> SeqError;
        std::vector<long double> ADOError;
        
        std::vector<long double> currentDeltas;
        std::vector<long double> currentTs;
        std::vector<long double> currentThetas;
        std::vector<long double> currentSeqError;
        std::vector<long double> currentADOError;
        
        vector<double> *normalized_weights = currenPop->get_normalized_weights();
        double max = -DOUBLE_INF;
        shared_ptr<State> best_particle;
        for (size_t i=0; i<num_particles; i++){
            State &s = pop0[0].get_particle(i);
            shared_ptr<State> currents = particles->at(i);
            // Population *pop= s.getPopulationByIndex(0);
            deltas0.push_back(s.getPopulationByIndex(0)->delta);
            Ts0.push_back(s.getPopulationByIndex(0)->timeOriginSTD);
            Thetas.push_back(s.getTheta());
            ADOError.push_back(s.getErrorModel().getADOErrorRate());
            SeqError.push_back(s.getErrorModel().getSeqErrorRate());
            
            currentDeltas.push_back(currents->getPopulationByIndex(0)->delta);
            currentTs.push_back(currents->getPopulationByIndex(0)->timeOriginSTD);
            currentThetas.push_back(currents->getTheta());
            currentSeqError.push_back(currents->getErrorModel().getADOErrorRate());
            currentADOError.push_back(currents->getErrorModel().getSeqErrorRate());
            if ((*normalized_weights)[i] > max) {
                max = (*normalized_weights)[i]  ;
                best_particle = currents;
            }
        }
        std::sort(deltas0.begin(), deltas0.end());
        std::cout<< "Prior distribution" << std::endl;
        std::cout<< "Delta, mean: " << Utils::mean(deltas0) <<" var: " <<Utils::variance(deltas0) << std::endl;
        std::cout<< "T, mean: " << Utils::mean(Ts0) <<" var: " <<Utils::variance(Ts0) <<std::endl;
        std::cout<< "Theta, mean: " << Utils::mean(Thetas) <<" var: " <<Utils::variance(Thetas) <<std::endl;
        std::cout<< "SeqError, mean: " << Utils::mean(SeqError) <<" var: " <<Utils::variance(SeqError) <<std::endl;
        std::cout<< "ADOError, mean: " << Utils::mean(ADOError) <<" var: " <<Utils::variance(ADOError) <<std::endl;
        
        std::cout<< "Posterior distribution" << std::endl;
        std::cout<< "Delta, mean: " << Utils::mean(currentDeltas) <<" var: " <<Utils::variance(currentDeltas) << std::endl;
        std::cout<< "T, mean: " << Utils::mean(currentTs) <<" var: " <<Utils::variance(currentTs) <<std::endl;
        std::cout<< "Theta, mean: " << Utils::mean(currentThetas) <<" var: " <<Utils::variance(currentThetas) <<std::endl;
        std::cout<< "SeqError, mean: " << Utils::mean(currentSeqError) <<" var: " <<Utils::variance(currentSeqError) <<std::endl;
        std::cout<< "ADOError, mean: " << Utils::mean(currentADOError) <<" var: " <<Utils::variance(currentADOError) <<std::endl;
        
        
        
        double log_marginal_lik = smc.get_log_marginal_likelihood();
        cout << "Estimate log marginal " << log_marginal  << endl;
        cout << "Estimate log P(y)= " << log_marginal_lik  << endl;
        
        assert(best_particle->getRoots().size() == 1);
        best_particle->printTree(best_particle->getRoots()[0], std::cerr);
        
    }
    
    /* clean memory*/
    
    Random::freeListRandomNumbersGenerators(randomGenerators);
    trueTrees.clear();
    structuredCoalTrees.clear();
    trueThetas.clear();
    trueDeltaTs.clear();
    trueTs.clear();
    sampleSizes.clear();
    
    pll_msa_destroy(msa);
    std::cout << "\nIf you need help type '-?' in the command line of the program\n"<<std::endl;
    
    return 0;
    
}


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



//#include <stan⁩/lib⁩/stan_math⁩/stan/math/prim/fun/Eigen.hpp>

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
#include "msa.hpp"
#include "tree.hpp"
#include "pmmh.hpp"
#include "RF_distance_calculator.hpp"

#include <gsl/gsl_randist.h>
#include <boost/program_options.hpp>

#include "smc.hpp"
#include "state.hpp"
#include "poset_smc.hpp"
#include "bd_coal_proposal.hpp"
#include "ipmcmc.hpp"
#include "pg.hpp"
#include "pmmh.hpp"

//#include "stan_model_1_population_JC_genotypes.hpp"




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

enum InferenceMethod {
    PARTICLEGIBBS = 0,
    IPMCMC = 1,
    PMMH = 2,
    PHMC = 3
};
namespace po = boost::program_options;
bool process_command_line(int argc, char** argv,
                          std::string& config_file)
{
    po::options_description desc("Program options");
    po::variables_map vm;
    try
    {
        desc.add_options()
        ("help", "produce help message.")
        ("config_file,c", po::value<std::string>(&config_file), "path to configuration file.")
        ;
        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        
        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }
        
        
        po::notify(vm);
    }
    catch (const boost::program_options::required_option & e) {
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 1;
        } else {
            throw e;
        }
    }
    catch(...)
    {
        std::cerr << "Unknown error!" << "\n";
        return false;
    }
    
    std::stringstream ss;
    
    
    return true;
}


int main(int argc, char* argv[] )
{
    FILE *input_file;
    const char* input_path;
    std::string config_file;
    
    bool result = process_command_line(argc, argv,config_file);
    if (!result)
        return 1;
    
    input_path = config_file.c_str();
    
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCOptions mcmcOptions;
    
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
    
    pll_msa_t *msa;
    
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    Utils::ReadParametersFromFastaFile(fileNameFasta,  programOptions.numCells, programOptions.TotalNumSequences, programOptions.numSites);
    std::vector<std::vector<int> > ObservedData;
    
    //fasta input
    //  char *ObservedCellNames[programOptions.numCells];
    //    Utils::ReadFastaFile(fileNameFasta, ObservedData,  ObservedCellNames, programOptions);
    //    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    //    programOptions.numCells = programOptions.TotalNumSequences;
    //
    //    programOptions.TotalTumorSequences=programOptions.TotalNumSequences-1;
    
    fileNamePhylip =filePaths.inputGenotypeFilePhylip;
    msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
    if (msa==NULL)
        std::cout << "Error reading phylip file \n"<< std::endl;
    // }
    MSA msaWrapper(msa);
    if (!MSA::checkMSA(msaWrapper)){
        
        std::cout << "Error cheking the sequences alignment \n"<< std::endl;
    }
    
    programOptions.numClones = 1;
    programOptions.TotalNumSequences = msa->count;
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    programOptions.numCells = programOptions.TotalNumSequences;
    
    programOptions.TotalTumorSequences=programOptions.TotalNumSequences-1;
    
    if (programOptions.numClones==1)
        sampleSizes.push_back(programOptions.numCells-1); //minus the healthycell
    
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
    
    //pll_rnode_t * healtyTip = initialRootedTree->root->right;
    char * newickWithoutHealthy = pll_rtree_export_newick(initialRootedTree->root->left, NULL);
    std::cout << "\n True tree newick: \n " << newickWithoutHealthy<< std::endl;
    RootedTree rootedTree(newickWithoutHealthy, false);
    
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
    
    programOptions.meanADOsite = 0.1;
    programOptions.varADOsite=0.01;
    programOptions.meanADOcell = 0.1;
    programOptions.varADOcell=0.01;
    programOptions.meanGenotypingError= 0.05;
    programOptions.varGenotypingError=0.001;
    programOptions.fixedADOrate=0.1;
    
    std::vector<int> positions(programOptions.TotalTumorSequences);
    std::iota( std::begin( positions ), std::end( positions ), 1 );
    //std::random_shuffle( positions.begin(), positions.end());
    
    // GenotypeErrorModel *gtErrorModel= new GenotypeErrorModel("GT20", programOptions.meanGenotypingError,  1.0 - sqrt (1.0 - programOptions.fixedADOrate), 16);
    
    GenotypeErrorModel *gtErrorModel= new GenotypeErrorModel("GT20", 0.0, 0.0, 16);
    
    PLLBufferManager *pll_buffer_manager = new PLLBufferManager;
    //        const pll_partition_t* pll_partition= pll_utils::createGTReferencePartition(msa);
    
    std::vector<double> coalTimes = pll_utils::getOrderedCoalTimesFromRootedTree(initialRootedTree, healthyTipLabel);
    
    pll_utils::printChronologicalOrderRootedTree(initialRootedTree,healthyTipLabel, std::cerr);
    
    Partition * partition = new Partition(msa,
                                          16,// model->states,//numberStates
                                          1,//RATE_CATS, // unsigned  int  numberRateCats
                                          0, //int statesPadded
                                          true,//PLL_ATTRIB_ARCH_SSE
                                          false, false, false, false, false);
    double theta = 1;
    
    double TscaledByTheta = coalTimes[coalTimes.size()-1];
    
    std::vector<double> timeOriginSTDs = { TscaledByTheta/theta };
    
    std::vector<double> deltas  = {100};//{21.519};
    
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
    coalTimesModelTimePerPopulation.push_back(coalTimes);
    
    
    PosetSMCParams psParams(programOptions.numClones, programOptions.TotalNumSequences,  sampleSizes,programOptions.numSites, &msaWrapper, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                            theta,
                            deltas,
                            timeOriginSTDs,
                            {1.0},
                            coalTimesModelTimePerPopulation);
    psParams.doPriorPost= false;
    psParams.doFixedEventimes= false;
    psParams.usePriorInSMC1 = true;
    bool normalizedCLVs = false;
    
    
    size_t num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    PosetSMC posetSMC(programOptions.numClones,  num_iter, false);
    
    posetSMC.kernelType = PosetSMC::PosetSMCKernel::TSMC1;

    
    if (posetSMC.kernelType== PosetSMC::PosetSMCKernel::TSMC1){
        
        if (!(psParams.usePriorInSMC1)){
            
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
    
    //  smc options
    SMCOptions smcOptions;
    smcOptions.num_threads = 5;
    smcOptions.use_SPF = false;
    smcOptions.ess_threshold = 1;
    smcOptions.num_particles = 100;
    smcOptions.resample_last_round = false;
    smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::MULTINOMIAL;
    //smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::STRATIFIED;
    //smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::SYSTEMATIC;
    //  smcOptions.main_seed = 346435;
    // smcOptions.resampling_seed = 2345666;
    smcOptions.track_population = false;
    smcOptions.init();
    smcOptions.debug = true;
    
    
    PMCMCOptions pmcmc_options(22441453521, 10000);
    pmcmc_options.burn_in = 1000;
    
    BDCoalPriorParams priorParams(10, 1, programOptions.meanGenotypingError,
                                  programOptions.varGenotypingError,
                                  programOptions.meanAmplificationError,
                                  programOptions.varAmplificationError,
                                  0.7, 0.7);
    
    //3. do inference
    InferenceMethod method = PMMH;
    
    if (method == IPMCMC){
        //Interactive PMCMC
        //TODO
        
        
    
    }
    else if(method == PARTICLEGIBBS){
        //Particle Gibbs
        ConditionalSMC<State, PosetSMCParams> csmc(posetSMC, smcOptions);
        
        BDCoalModelPGProposal param_proposal(priorParams, programOptions.numClones, programOptions.TotalNumSequences,
                                             msaWrapper.getLength(),
                                             &msaWrapper,
                                             partition,
                                             pll_buffer_manager,
                                             programOptions);
        

        ParticleGibbs<State, PosetSMCParams> pg(pmcmc_options, csmc, param_proposal);
        
        pg.run();
        vector<shared_ptr<PosetSMCParams> > &samples = pg.get_parameters();
        //   compute the posterior mean for the first scaled hrowth rate
        double mean = 0.0;
        size_t count = 0;
        for (size_t i = pmcmc_options.burn_in; i < samples.size(); i+=10) {
            mean += samples[i]->populationDeltaTs[0];
            count++;
        }
        mean /= count;
        cout << mean << ", " << psParams.populationDeltaTs[0] << endl;
    }
    else if(method == PMMH){
        //Particle Marginal Metropolis-Hastings
        smcOptions.num_particles = 8;
     
        
        ConditionalSMC<State, PosetSMCParams> csmc(posetSMC, smcOptions);
         BDCoalModelRandomWalkProposal rw_param_proposal(priorParams, programOptions.numClones, programOptions.TotalNumSequences,
                                                    msaWrapper.getLength(),
                                                    &msaWrapper,
                                                    partition,
                                                    pll_buffer_manager,
                                                    programOptions);
        
          csmc.initialize(psParams);
          double logZ = csmc.get_log_marginal_likelihood();
          std::cout << logZ << std::endl; // logZ at true params

          PMCMCOptions pmcmc_options(346575392, 200);
          pmcmc_options.burn_in = 100;
         
          ParticleMMH<State, PosetSMCParams> pmmh(pmcmc_options, csmc, rw_param_proposal);
          pmmh.run();
          vector<PosetSMCParams *> *samples = pmmh.get_parameters();
          // compute the posterior mean for beta
          double mean = 0.0;
          size_t count = 0;
          for (size_t i = pmcmc_options.burn_in; i < samples->size(); i+=10) {
              mean += (*samples)[i]->populationDeltaTs[0];
              count++;
          }
          mean /= count;
        cout << mean << ", " << psParams.populationDeltaTs[0] << endl;
        
        
    }
    else{
        //Particle Hamiltonian Monte Carlo
        //TODO
        
        
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
    pll_rtree_destroy(initialRootedTree,NULL);
    std::cout << "\nIf you need help type '-?' in the command line of the program\n"<<std::endl;
    

    return 0;
    
}


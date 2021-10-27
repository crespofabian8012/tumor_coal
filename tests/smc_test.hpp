//
//  smc_test.hpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 31/08/2021.
//

#ifndef smc_test_hpp
#define smc_test_hpp

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
//#include "chain_manager.hpp"
#include "definitions.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "output_functions.hpp"
//#include "mcmc.hpp"
#include "msa.hpp"
#include "tree.hpp"
#include "RF_distance_calculator.hpp"

#include <gsl/gsl_randist.h>
#include <Eigen/Eigen>
#include <boost/program_options.hpp>
#include <chrono>
#include <sys/time.h>

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

#include "gtest/gtest.h"

using SMC_State_PosetSMCParams = SMC<State, PosetSMCParams>;
class PosetSMCOnePopulationTest : public ::testing::TestWithParam<std::tuple<std::string, std::string, size_t,SMCOptions::ResamplingScheme, PosetSMC::PosetSMCKernel, double, double, bool  >> {
protected:
    
public:
    // PosetSMCTest(std::string inputGenotypeFilePhylip, std::string inputTreePath):
    PosetSMCOnePopulationTest()
    {
        inputGenotypeFilePhylipPath = std::get<0>(GetParam());
        inputTreePath = std::get<1>(GetParam());
        bool doPlots = std::get<7>(GetParam());
        
        smcOptions.num_threads = 5;
        smcOptions.use_SPF = false;
        smcOptions.ess_threshold = 1;
        smcOptions.num_particles = std::get<2>(GetParam());
        smcOptions.resample_last_round = false;
        smcOptions.resampling_scheme =  std::get<3>(GetParam());
        smcOptions.track_population = false;
        using namespace std::chrono;
        milliseconds ms = duration_cast< milliseconds >(
                                                        system_clock::now().time_since_epoch());
        long seed = ms.count();;
        random = Random::generateRandomObject( seed);
        struct timeval tv;
        gettimeofday(&tv,0);
        seed = tv.tv_sec + tv.tv_usec;
       // smcOptions.main_seed = seed;
        
        const char* input_path;
        const std::string config_file = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/parametersMCMC.txt";
        
        input_path = config_file.c_str();
        
        setDefaultOptions(programOptions, mcmcOptions);
        InitFilesPathsOptions(filePaths, programOptions);
        
        
        char *fileNameFasta ;//= argv[2];
        char *fileNamePhylip;
        char *treefileName; //= argv[4];
        
        pll_msa_t *msa;
        fileNameFasta = filePaths.inputGenotypeFileFasta;
        
        strcpy(filePaths.inputGenotypeFilePhylip, inputGenotypeFilePhylipPath.c_str());
        
        fileNamePhylip =filePaths.inputGenotypeFilePhylip;
        msa = pll_phylip_load(fileNamePhylip, PLL_FALSE);
        if (msa==NULL)
            std::cout << "Error reading phylip file \n"<< std::endl;
        // }
        MSA *msaWrapper = new MSA(msa);
        //MSA msaWrapper(msa);
        if (!MSA::checkMSA(*msaWrapper)){
            
            std::cout << "Error cheking the sequences alignment \n"<< std::endl;
        }
        
        programOptions.numClones = 1;
        programOptions.TotalNumSequences = msa->count;
        programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
        programOptions.numCells = programOptions.TotalNumSequences;
        programOptions.TotalTumorSequences=programOptions.TotalNumSequences-1;
        
        if (programOptions.numClones==1)
            sampleSizes.push_back(programOptions.numCells-1); //minus the healthycell
        
        strcpy(filePaths.inputTreeFile, inputTreePath.c_str());
        treefileName = filePaths.inputTreeFile;
        
        if ((treefileName != NULL) && (treefileName[0] == '\0'))
        {
            std::cout << "\nERROR: No tree file  specified (use command line or parameter file)" << std::endl;
            exit(-1);
        }
        std::string healthyTipLabel = "healthycell";
        programOptions.healthyTipLabel ="healthycell";
        
        std::vector<gsl_rng * > randomGenerators(mcmcOptions.numChains);
        
        
        initialRootedTree = pll_rtree_parse_newick(treefileName);
        char * newickWithoutHealthy = pll_rtree_export_newick(initialRootedTree->root->left, NULL);
        std::cout << "\n True tree newick: \n " << newickWithoutHealthy<< std::endl;
        rootedTree = new RootedTree(newickWithoutHealthy, false);
        
        programOptions.seqErrorRate = 0.01;
        programOptions.dropoutRate =0.2;
        mcmcOptions.noData = false;
        mcmcOptions.useGSLRandomGenerator = true;
        mcmcOptions.splitThetaDeltaTmoves=true;
        programOptions.doUseFixedTree=NO;
        
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
        
        GenotypeErrorModel *gtErrorModel= new GenotypeErrorModel("GT20", 0.0, 0.0, 16);
        
        PLLBufferManager *pll_buffer_manager = new PLLBufferManager;
        //        const pll_partition_t* pll_partition= pll_utils::createGTReferencePartition(msa);
        
        std::vector<double> coalTimes = pll_utils::getOrderedCoalTimesFromRootedTree(initialRootedTree, healthyTipLabel);
        
        std::cout << "\n True tree chronological order: \n " << std::endl;
        pll_utils::printChronologicalOrderRootedTree(initialRootedTree,healthyTipLabel, std::cerr);
        
        Partition * partition = new Partition(msa,
                                              16,// model->states,//numberStates
                                              1,//RATE_CATS, // unsigned  int  numberRateCats
                                              0, //int statesPadded
                                              true,//PLL_ATTRIB_ARCH_SSE
                                              false, false, false, false, false);
        true_theta =  std::get<5>(GetParam());
        double TscaledByTheta = coalTimes[coalTimes.size()-1];
        std::vector<double> deltas  = {std::get<6>(GetParam())};
        std::vector<double> timeOriginSTDs = { TscaledByTheta/true_theta };
        std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
         
        const double theta = true_theta;
        transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
        
        coalTimesModelTimePerPopulation.push_back(coalTimes);
        
        
        psParams = new  PosetSMCParams(programOptions.numClones, programOptions.TotalNumSequences,  sampleSizes,programOptions.numSites, msaWrapper, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                                true_theta,
                                deltas,
                                timeOriginSTDs,
                                {1.0},
                                coalTimesModelTimePerPopulation);
        psParams->doPriorPost = true;
        psParams->doFixedEventimes = false;
        psParams->usePriorInSMC1 = false;
        bool normalizedCLVs = true;
        
        num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
        posetSMC = new PosetSMC(programOptions.numClones,  num_iter, doPlots);
        
        posetSMC->kernelType = std::get<4>(GetParam());
        
        if (posetSMC->kernelType== PosetSMC::PosetSMCKernel::TSMC1){
            
            if (!(psParams->usePriorInSMC1)){
                
                if (normalizedCLVs){
                    programOptions.normalizeLeavesClv =true;
                    programOptions.normalizeClv =true;
                     
                    
                }
                else{
                    programOptions.normalizeLeavesClv =false;
                    programOptions.normalizeClv =false;
                    
                }
            }
        }
            
            
        smcOptions.init();
        smcOptions.debug = true;
        
        smc = new SMC<State, PosetSMCParams>(*posetSMC, smcOptions);
    }
    
    ~PosetSMCOnePopulationTest() override {
      //  pll_msa_destroy(msa);
      //  pll_rtree_destroy(initialRootedTree,NULL);
    }
  
    //    void SetUp() override {
    //
    //
    //    }
    
   
    
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip;
    char *treefileName; //= argv[4];
    std::string inputGenotypeFilePhylipPath;
    std::string inputTreePath ;
    std::vector<int> sampleSizes;
    gsl_rng * random;
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCOptions mcmcOptions;
    pll_msa_t *msa;
    pll_rtree_t *initialRootedTree;
    SMC_State_PosetSMCParams *smc;
    int num_particles;
    size_t  num_iter;
    SMCOptions smcOptions;
    PosetSMCParams *psParams;
    PosetSMC *posetSMC;
    //MSA *msa;
    RootedTree *rootedTree;
    double true_theta;
};


#endif /* smc_test_hpp */

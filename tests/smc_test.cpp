
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



void test_smc_1pop_4tips( )
{
    std::cout << "\nSMC tests with 4 tips\n" << std::endl;
    const char* input_path;
    std::string config_file = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/parametersMCMC.txt";
    
    input_path = config_file.c_str();
    
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    
    
    setDefaultOptions(programOptions, mcmcOptions);
    
    
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip;
    char *treefileName; //= argv[4];
    
    std::vector<int> sampleSizes;
    
    //2. initialize data structures
    //    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    
    pll_msa_t *msa;
    
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    
    
    std::string inputGenotypeFilePhylip ="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001.phylip";
    
    strcpy(filePaths.inputGenotypeFilePhylip, inputGenotypeFilePhylip.c_str());
    
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
    
    
    std::string inputTreePath = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Model_time_0002_0010.tre";
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
    
    
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    char * newickWithoutHealthy = pll_rtree_export_newick(initialRootedTree->root->left, NULL);
    std::cout << "\n True tree newick: \n " << newickWithoutHealthy<< std::endl;
    RootedTree rootedTree(newickWithoutHealthy, false);
    
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
    double theta = 2.5;//0.019;
    double TscaledByTheta = coalTimes[coalTimes.size()-1];
    std::vector<double> deltas  = {100};//{21.519};
    std::vector<double> timeOriginSTDs = { TscaledByTheta/theta };
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
    
    coalTimesModelTimePerPopulation.push_back(coalTimes);
    
    
    PosetSMCParams psParams(programOptions.numClones, programOptions.TotalNumSequences,  sampleSizes,programOptions.numSites, msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                            theta,
                            deltas,
                            timeOriginSTDs,
                            {1.0},
                            coalTimesModelTimePerPopulation);
    psParams.doPriorPost= true;
    psParams.doFixedEventimes= false;
    
    size_t num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    PosetSMC posetSMC(programOptions.numClones,  num_iter);
    SMCOptions smcOptions;
    
    smcOptions.num_threads = 5;
    smcOptions.use_SPF = false;
    smcOptions.ess_threshold = 1;
    smcOptions.num_particles = 1000;
    smcOptions.resample_last_round = false;
    
    smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::MULTINOMIAL;
    
    smcOptions.track_population = false;
    smcOptions.init();
    smcOptions.debug = true;
    
    SMC<State, PosetSMCParams>  smc(posetSMC, smcOptions);
    
    std::cout<< "\n\nRunning Sequential Monte Carlo(SMC)" << " with "<<smcOptions.num_particles <<" particles....\n" << std::endl;
    smc.run_smc(psParams);
    
    ParticlePopulation<State> *currenPop = smc.get_curr_population();
    vector<shared_ptr<State>> *particles = currenPop->get_particles();
    double   log_marginal = smc.get_log_marginal_likelihood();
    //ParticlePopulation<State> *pop0 = smc.get_population(0);
    
    int num_particles = currenPop->get_num_particles();
    
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
    std::vector<long double> weights;
    std::vector<long double> rootLogLiks;
    
    std::vector<double> *normalized_weights = currenPop->get_normalized_weights();
    double max = -DOUBLE_INF;
    shared_ptr<State> best_particle;
    
    Eigen::ArrayXXd allCoalTimesScaledByTheta(num_particles,num_iter-1);
    
    std::vector<RootedTree> lastPopulationTrees;
    for (size_t i=0; i<num_particles; i++){
        
        shared_ptr<State> currents = particles->at(i);
        currentDeltas.push_back(currents->getPopulationByIndex(0)->delta);
        currentTs.push_back(currents->getPopulationByIndex(0)->timeOriginSTD);
        currentThetas.push_back(currents->getTheta());
        currentSeqError.push_back(currents->getErrorModel().getADOErrorRate());
        currentADOError.push_back(currents->getErrorModel().getSeqErrorRate());
        
        //currents->printTree(currents->getRoots()[0], std::cerr);
        
        // std::cout << "Particle weight "<< i << " " <<(*normalized_weights)[i] <<std::endl;
        assert(currents->getRoots().size() == 1);
        weights.push_back((*normalized_weights)[i]);
        rootLogLiks.push_back(currents->getRootAt(0)->ln_likelihood);
        //std::cout << "Root loglik "<< i << " " <<currents->getRootAt(0)->ln_likelihood <<std::endl;
        allCoalTimesScaledByTheta.row(i) = Eigen::Map<Eigen::RowVectorXd>(currents->getCoalEventTimes().data(),num_iter-1 );
        //currents->printTreeChronologicalOrder(currents->getRoots()[0],std::cerr);
        
        if ((*normalized_weights)[i] > max) {
            
            max = (*normalized_weights)[i]  ;
            best_particle = currents;
        }
        std::string newick = currents->getNewick(currents->getRootAt(0));
        lastPopulationTrees.emplace_back(RootedTree(newick, false));
    }
    lastPopulationTrees.emplace_back(rootedTree);
    
    std::unique_ptr<RFDistanceCalculator> rfCalculator;
    rfCalculator.reset(new RFDistanceCalculator(lastPopulationTrees, false));
    
    
    std::cout<< "Posterior distribution" << std::endl;
    std::cout<< "Delta, mean: " << Utils::mean(currentDeltas) <<" var: " <<Utils::variance(currentDeltas) << std::endl;
    std::cout<< "T, mean: " << Utils::mean(currentTs) <<" var: " <<Utils::variance(currentTs) <<std::endl;
    std::cout<< "Theta, mean: " << Utils::mean(currentThetas) <<" var: " <<Utils::variance(currentThetas) <<std::endl;
    std::cout<< "SeqError, mean: " << Utils::mean(currentSeqError) <<" var: " <<Utils::variance(currentSeqError) <<std::endl;
    std::cout<< "ADOError, mean: " << Utils::mean(currentADOError) <<" var: " <<Utils::variance(currentADOError) <<std::endl;
    std::cout<< "Normalized weight " << Utils::mean(weights) <<" var: " <<Utils::variance(weights) <<std::endl;
    std::cout<< "Root log liks " << Utils::mean(rootLogLiks) <<" var: " <<Utils::variance(rootLogLiks) <<std::endl;
    
    double avgRF = rfCalculator->avgRF();
    
    size_t numUniqueTopologies = rfCalculator->numUniqueTrees();
    
    std::cout << "Avg RF " << avgRF  << std::endl;
    std::cout << "Number unique topologies " << numUniqueTopologies  << std::endl;
    
    EXPECT_EQ(avgRF,0.0);
    assert(abs(avgRF-0) <0.0000000001);
    assert(numUniqueTopologies ==1 );
    
    Eigen::ArrayXXd allCoalTimes =allCoalTimesScaledByTheta *(1./theta) ;
    
    Eigen::Array<double, 1, Eigen::Dynamic> means = allCoalTimes.colwise().mean();
    Eigen::Array<double, 1, Eigen::Dynamic> mins = allCoalTimes.colwise().minCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> maxs = allCoalTimes.colwise().maxCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> std_dev = ((allCoalTimes.rowwise() - means.array()).square().colwise().sum()/(num_iter-2)).sqrt();
    for (size_t i=0; i<(num_iter-1); i++){
        
        std::cout << "coal time " << i << " min, max: "<< mins(i)<< ","<< maxs(i) << " mean: "<< means(i) <<" std: " << std_dev(i) << " true "<< coalTimes[i]<< "\n"<<std::endl;
        
    }
    
    
    double log_marginal_lik = smc.get_log_marginal_likelihood();
    std::cout << "Estimate log marginal " << log_marginal  << endl;
    std::cout << "Estimate log P(y)= " << log_marginal_lik  << endl;
    
    assert(best_particle->getRoots().size() == 1);
    std::cout << "Tree topology with the greater loglik: "  << endl;
    best_particle->printTree(best_particle->getRoots()[0], std::cerr);
    
    /* clean memory*/
    
    pll_msa_destroy(msa);
    pll_rtree_destroy(initialRootedTree,NULL);
    
}

void test_smc_1pop_15tips( )
{
    
    std::cout << "\nSMC tests with 15 tips\n" << std::endl;
    const char* input_path;
    std::string config_file = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/parametersMCMC.txt";
    
    input_path = config_file.c_str();
    
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    
    
    setDefaultOptions(programOptions, mcmcOptions);
    
    
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip;
    char *treefileName; //= argv[4];
    
    std::vector<int> sampleSizes;
    
    //2. initialize data structures
    //    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    
    pll_msa_t *msa;
    
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    
    
    std::string inputGenotypeFilePhylip ="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001_n15.phylip";
    
    strcpy(filePaths.inputGenotypeFilePhylip, inputGenotypeFilePhylip.c_str());
    
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
    
    
    std::string inputTreePath = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n15.tre";
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
    
    
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    char * newickWithoutHealthy = pll_rtree_export_newick(initialRootedTree->root->left, NULL);
    std::cout << "\n True tree newick: \n " << newickWithoutHealthy<< std::endl;
    RootedTree rootedTree(newickWithoutHealthy, false);
    
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
    double theta = 2.5;//0.019;
    double TscaledByTheta = coalTimes[coalTimes.size()-1];
    std::vector<double> deltas  = {100};//{21.519};
    std::vector<double> timeOriginSTDs = { TscaledByTheta/theta };
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
    
    coalTimesModelTimePerPopulation.push_back(coalTimes);
    
    
    PosetSMCParams psParams(programOptions.numClones, programOptions.TotalNumSequences,  sampleSizes,programOptions.numSites, msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                            theta,
                            deltas,
                            timeOriginSTDs,
                            {1.0},
                            coalTimesModelTimePerPopulation);
    psParams.doPriorPost= true;
    psParams.doFixedEventimes= false;
    
    size_t num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    PosetSMC posetSMC(programOptions.numClones,  num_iter);
    SMCOptions smcOptions;
    
    smcOptions.num_threads = 5;
    smcOptions.use_SPF = false;
    smcOptions.ess_threshold = 1;
    smcOptions.num_particles = 1000;
    smcOptions.resample_last_round = false;
    
    smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::MULTINOMIAL;
    
    smcOptions.track_population = false;
    smcOptions.init();
    smcOptions.debug = true;
    
    SMC<State, PosetSMCParams>  smc(posetSMC, smcOptions);
    
    std::cout<< "\n\nRunning Sequential Monte Carlo(SMC)" << " with "<<smcOptions.num_particles <<" particles....\n" << std::endl;
    smc.run_smc(psParams);
    
    ParticlePopulation<State> *currenPop = smc.get_curr_population();
    vector<shared_ptr<State>> *particles = currenPop->get_particles();
    double   log_marginal = smc.get_log_marginal_likelihood();
    //ParticlePopulation<State> *pop0 = smc.get_population(0);
    
    int num_particles = currenPop->get_num_particles();
    
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
    std::vector<long double> weights;
    std::vector<long double> rootLogLiks;
    
    std::vector<double> *normalized_weights = currenPop->get_normalized_weights();
    double max = -DOUBLE_INF;
    shared_ptr<State> best_particle;
    
    Eigen::ArrayXXd allCoalTimesScaledByTheta(num_particles,num_iter-1);
    
    std::vector<RootedTree> lastPopulationTrees;
    for (size_t i=0; i<num_particles; i++){
        
        shared_ptr<State> currents = particles->at(i);
        currentDeltas.push_back(currents->getPopulationByIndex(0)->delta);
        currentTs.push_back(currents->getPopulationByIndex(0)->timeOriginSTD);
        currentThetas.push_back(currents->getTheta());
        currentSeqError.push_back(currents->getErrorModel().getADOErrorRate());
        currentADOError.push_back(currents->getErrorModel().getSeqErrorRate());
        
        //currents->printTree(currents->getRoots()[0], std::cerr);
        // std::cout << "Particle weight "<< i << " " <<(*normalized_weights)[i] <<std::endl;
        assert(currents->getRoots().size() == 1);
        weights.push_back((*normalized_weights)[i]);
        rootLogLiks.push_back(currents->getRootAt(0)->ln_likelihood);
        // std::cout << "Root loglik "<< i << " " <<currents->getRootAt(0)->ln_likelihood <<std::endl;
        allCoalTimesScaledByTheta.row(i) = Eigen::Map<Eigen::RowVectorXd>(currents->getCoalEventTimes().data(),num_iter-1 );
        //currents->printTreeChronologicalOrder(currents->getRoots()[0],std::cerr);
        
        if ((*normalized_weights)[i] > max) {
            
            max = (*normalized_weights)[i]  ;
            best_particle = currents;
        }
        std::string newick = currents->getNewick(currents->getRootAt(0));
        lastPopulationTrees.emplace_back(RootedTree(newick, false));
    }
    lastPopulationTrees.emplace_back(rootedTree);
    
    std::unique_ptr<RFDistanceCalculator> rfCalculator;
    rfCalculator.reset(new RFDistanceCalculator(lastPopulationTrees, false));
    
    
    std::cout<< "Posterior distribution" << std::endl;
    std::cout<< "Delta, mean: " << Utils::mean(currentDeltas) <<" var: " <<Utils::variance(currentDeltas) << std::endl;
    std::cout<< "T, mean: " << Utils::mean(currentTs) <<" var: " <<Utils::variance(currentTs) <<std::endl;
    std::cout<< "Theta, mean: " << Utils::mean(currentThetas) <<" var: " <<Utils::variance(currentThetas) <<std::endl;
    std::cout<< "SeqError, mean: " << Utils::mean(currentSeqError) <<" var: " <<Utils::variance(currentSeqError) <<std::endl;
    std::cout<< "ADOError, mean: " << Utils::mean(currentADOError) <<" var: " <<Utils::variance(currentADOError) <<std::endl;
    std::cout<< "Normalized weight " << Utils::mean(weights) <<" var: " <<Utils::variance(weights) <<std::endl;
    std::cout<< "Root log liks " << Utils::mean(rootLogLiks) <<" var: " <<Utils::variance(rootLogLiks) <<std::endl;
    
    double avgRF = rfCalculator->avgRF();
    
    size_t numUniqueTopologies = rfCalculator->numUniqueTrees();
    
    std::cout << "Avg RF " << avgRF  << std::endl;
    std::cout << "Number unique topologies " << numUniqueTopologies  << std::endl;
    
    
    Eigen::ArrayXXd allCoalTimes =allCoalTimesScaledByTheta *(1./theta) ;
    
    Eigen::Array<double, 1, Eigen::Dynamic> means = allCoalTimes.colwise().mean();
    Eigen::Array<double, 1, Eigen::Dynamic> mins = allCoalTimes.colwise().minCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> maxs = allCoalTimes.colwise().maxCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> std_dev = ((allCoalTimes.rowwise() - means.array()).square().colwise().sum()/(num_iter-2)).sqrt();
    for (size_t i=0; i<(num_iter-1); i++){
        
        std::cout << "coal time " << i << " min, max: "<< mins(i)<< ","<< maxs(i) << " mean: "<< means(i) <<" std: " << std_dev(i) << " true "<< coalTimes[i]<< "\n"<<std::endl;
        
    }
    
    
    double log_marginal_lik = smc.get_log_marginal_likelihood();
    std::cout << "Estimate log marginal " << log_marginal  << endl;
    std::cout << "Estimate log P(y)= " << log_marginal_lik  << endl;
    
    
    assert(numUniqueTopologies ==2 );
    
    assert(best_particle->getRoots().size() == 1);
    std::cout << "Tree topology with the greater loglik: "  << endl;
    best_particle->printTree(best_particle->getRoots()[0], std::cerr);
    
    /* clean memory*/
    
    pll_msa_destroy(msa);
    pll_rtree_destroy(initialRootedTree,NULL);
    
}


void test_smc_1pop_25tips( )
{
    
    std::cout << "\nSMC tests with 25 tips\n" << std::endl;
    const char* input_path;
    std::string config_file = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/parametersMCMC.txt";
    
    input_path = config_file.c_str();
    
    ProgramOptions programOptions;
    FilePaths filePaths;
    MCMCoptions mcmcOptions;
    
    
    setDefaultOptions(programOptions, mcmcOptions);
    
    
    char *fileNameFasta ;//= argv[2];
    char *fileNamePhylip;
    char *treefileName; //= argv[4];
    
    std::vector<int> sampleSizes;
    
    //2. initialize data structures
    //    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    
    pll_msa_t *msa;
    
    fileNameFasta = filePaths.inputGenotypeFileFasta;
    
    
    std::string inputGenotypeFilePhylip ="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0001_0001.phylip";
    
    strcpy(filePaths.inputGenotypeFilePhylip, inputGenotypeFilePhylip.c_str());
    
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
    
    
    std::string inputTreePath = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Model_time_0002_0001.tre";
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
    
    
    
    pll_rtree_t * initialRootedTree = pll_rtree_parse_newick(treefileName);
    char * newickWithoutHealthy = pll_rtree_export_newick(initialRootedTree->root->left, NULL);
    std::cout << "\n True tree newick: \n " << newickWithoutHealthy<< std::endl;
    RootedTree rootedTree(newickWithoutHealthy, false);
    
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
    double theta = 2.5;//0.019;
    double TscaledByTheta = coalTimes[coalTimes.size()-1];
    std::vector<double> deltas  = {100};//{21.519};
    std::vector<double> timeOriginSTDs = { TscaledByTheta/theta };
    std::vector<std::vector<double>> coalTimesModelTimePerPopulation;
    
    transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
    
    coalTimesModelTimePerPopulation.push_back(coalTimes);
    
    
    PosetSMCParams psParams(programOptions.numClones, programOptions.TotalNumSequences,  sampleSizes,programOptions.numSites, msa, partition, pll_buffer_manager, positions, programOptions, gtErrorModel,
                            theta,
                            deltas,
                            timeOriginSTDs,
                            {1.0},
                            coalTimesModelTimePerPopulation);
    psParams.doPriorPost= true;
    psParams.doFixedEventimes= false;
    
    size_t num_iter = programOptions.TotalTumorSequences +programOptions.numClones-1;
    PosetSMC posetSMC(programOptions.numClones,  num_iter);
    SMCOptions smcOptions;
    
    smcOptions.num_threads = 5;
    smcOptions.use_SPF = false;
    smcOptions.ess_threshold = 1;
    smcOptions.num_particles = 1000;
    smcOptions.resample_last_round = false;
    
    smcOptions.resampling_scheme =  SMCOptions::ResamplingScheme::MULTINOMIAL;
    
    smcOptions.track_population = false;
    smcOptions.init();
    smcOptions.debug = true;
    
    SMC<State, PosetSMCParams>  smc(posetSMC, smcOptions);
    
    std::cout<< "\n\nRunning Sequential Monte Carlo(SMC)" << " with "<<smcOptions.num_particles <<" particles....\n" << std::endl;
    smc.run_smc(psParams);
    
    ParticlePopulation<State> *currenPop = smc.get_curr_population();
    vector<shared_ptr<State>> *particles = currenPop->get_particles();
    double   log_marginal = smc.get_log_marginal_likelihood();
    //ParticlePopulation<State> *pop0 = smc.get_population(0);
    
    int num_particles = currenPop->get_num_particles();
    
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
    std::vector<long double> weights;
    std::vector<long double> rootLogLiks;
    
    std::vector<double> *normalized_weights = currenPop->get_normalized_weights();
    double max = -DOUBLE_INF;
    shared_ptr<State> best_particle;
    
    Eigen::ArrayXXd allCoalTimesScaledByTheta(num_particles,num_iter-1);
    
    std::vector<RootedTree> lastPopulationTrees;
    for (size_t i=0; i<num_particles; i++){
        
        shared_ptr<State> currents = particles->at(i);
        currentDeltas.push_back(currents->getPopulationByIndex(0)->delta);
        currentTs.push_back(currents->getPopulationByIndex(0)->timeOriginSTD);
        currentThetas.push_back(currents->getTheta());
        currentSeqError.push_back(currents->getErrorModel().getADOErrorRate());
        currentADOError.push_back(currents->getErrorModel().getSeqErrorRate());
        
        //currents->printTree(currents->getRoots()[0], std::cerr);
        // std::cout << "Particle weight "<< i << " " <<(*normalized_weights)[i] <<std::endl;
        assert(currents->getRoots().size() == 1);
        weights.push_back((*normalized_weights)[i]);
        rootLogLiks.push_back(currents->getRootAt(0)->ln_likelihood);
        // std::cout << "Root loglik "<< i << " " <<currents->getRootAt(0)->ln_likelihood <<std::endl;
        allCoalTimesScaledByTheta.row(i) = Eigen::Map<Eigen::RowVectorXd>(currents->getCoalEventTimes().data(),num_iter-1 );
        //currents->printTreeChronologicalOrder(currents->getRoots()[0],std::cerr);
        
        if ((*normalized_weights)[i] > max) {
            
            max = (*normalized_weights)[i]  ;
            best_particle = currents;
        }
        std::string newick = currents->getNewick(currents->getRootAt(0));
        lastPopulationTrees.emplace_back(RootedTree(newick, false));
    }
    lastPopulationTrees.emplace_back(rootedTree);
    
    std::unique_ptr<RFDistanceCalculator> rfCalculator;
    rfCalculator.reset(new RFDistanceCalculator(lastPopulationTrees, false));
    
    
    std::cout<< "Posterior distribution" << std::endl;
    std::cout<< "Delta, mean: " << Utils::mean(currentDeltas) <<" var: " <<Utils::variance(currentDeltas) << std::endl;
    std::cout<< "T, mean: " << Utils::mean(currentTs) <<" var: " <<Utils::variance(currentTs) <<std::endl;
    std::cout<< "Theta, mean: " << Utils::mean(currentThetas) <<" var: " <<Utils::variance(currentThetas) <<std::endl;
    std::cout<< "SeqError, mean: " << Utils::mean(currentSeqError) <<" var: " <<Utils::variance(currentSeqError) <<std::endl;
    std::cout<< "ADOError, mean: " << Utils::mean(currentADOError) <<" var: " <<Utils::variance(currentADOError) <<std::endl;
    std::cout<< "Normalized weight " << Utils::mean(weights) <<" var: " <<Utils::variance(weights) <<std::endl;
    std::cout<< "Root log liks " << Utils::mean(rootLogLiks) <<" var: " <<Utils::variance(rootLogLiks) <<std::endl;
    
    double avgRF = rfCalculator->avgRF();
    
    size_t numUniqueTopologies = rfCalculator->numUniqueTrees();
    
    std::cout << "Avg RF " << avgRF  << std::endl;
    std::cout << "Number unique topologies " << numUniqueTopologies  << std::endl;
    
    
    Eigen::ArrayXXd allCoalTimes =allCoalTimesScaledByTheta *(1./theta) ;
    
    Eigen::Array<double, 1, Eigen::Dynamic> means = allCoalTimes.colwise().mean();
    Eigen::Array<double, 1, Eigen::Dynamic> mins = allCoalTimes.colwise().minCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> maxs = allCoalTimes.colwise().maxCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> std_dev = ((allCoalTimes.rowwise() - means.array()).square().colwise().sum()/(num_iter-2)).sqrt();
    for (size_t i=0; i<(num_iter-1); i++){
        
        std::cout << "coal time " << i << " min, max: "<< mins(i)<< ","<< maxs(i) << " mean: "<< means(i) <<" std: " << std_dev(i) << " true "<< coalTimes[i]<< "\n"<<std::endl;
        
    }
    
    
    double log_marginal_lik = smc.get_log_marginal_likelihood();
    std::cout << "Estimate log marginal " << log_marginal  << endl;
    std::cout << "Estimate log P(y)= " << log_marginal_lik  << endl;
    
    
    assert(numUniqueTopologies ==2 );
    
    assert(best_particle->getRoots().size() == 1);
    std::cout << "Tree topology with the greater loglik: "  << endl;
    best_particle->printTree(best_particle->getRoots()[0], std::cerr);
    
    /* clean memory*/
    
    pll_msa_destroy(msa);
    pll_rtree_destroy(initialRootedTree,NULL);
    
}

void test_smc( ){
    std::cout << "\n\nRunning  SMC tests...!\n" << std::endl;
    test_smc_1pop_4tips();
    test_smc_1pop_15tips();
    test_smc_1pop_25tips( );
    std::cout << "\n\nEnd of  SMC tests...!\n" << std::endl;
}

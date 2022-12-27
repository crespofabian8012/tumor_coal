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



#include "stan_model_1_population_JC_genotypes.hpp"
#include <stan/services/sample/hmc_nuts_diag_e_adapt.hpp>
#include <stan/services/sample/hmc_nuts_dense_e.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/logger.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/io/dump.hpp>

#include <iostream>
#include <fstream>

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

using SMC_State_PosetSMCParams = SMC<State, PosetSMCParams>;
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
    
    unsigned int random_seed = 64764894795;//= dynamic_cast<u_int_argument*>(parser.arg("random")->arg("seed"))->value();
    
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
    Files files;
    InitFilesPathsOptions(filePaths, programOptions);
    InitFiles(files);
    
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
    programOptions.healthyTipLabel = healthyTipLabel;
    
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
    
    posetSMC.kernelType =  PosetSMC::PosetSMCKernel::TSMC1;
    posetSMC.incrementProposalDist = PosetSMC::PosetSMCIncrementProposal::LBD_COAL;
    
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
    //smcOptions.main_seed = 346435;
    //smcOptions.resampling_seed = 2345666;
    smcOptions.track_population = false;
    smcOptions.init();
    smcOptions.debug = true;
    
    
    PMCMCOptions pmcmc_options(random_seed, 10000);
    pmcmc_options.burn_in = 1000;
    
    BDCoalPriorParams priorParams(10, 1, programOptions.meanGenotypingError,
                                  programOptions.varGenotypingError,
                                  programOptions.meanAmplificationError,
                                  programOptions.varAmplificationError,
                                  0.7, 0.7);
    
    //3. do inference
    InferenceMethod method = PHMC;
    
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
        std::vector<shared_ptr<PosetSMCParams> > &samples = pg.get_parameters();
        //   compute the posterior mean for the first scaled growth rate
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
        
        PMCMCOptions pmcmc_options(random_seed, 200);
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
        //1-sample initial parameters Gammas, Ts, Thetas,
        //2-run SMC to approximate the distributions over
        //topologies and branch lengths
        //3-run HMC
        //4-if HMC not converged, then run SMC with the median parameter values of current HMC(got to  step 2)
        
        int numIterCondictionalHMC = 10;
        int numTotalIter = 5;
        std::stringstream model_log;
        std::stringstream myfile;
        
        //auto filePath = File::getSpecialLocation (userHomeDirectory).getChildFile ("output250.R");
        //myfile.open(filePath.getFullPathName.toRawUTF8());
        
        string outputfile     = "⁨‎output.R";
        string diagnosticfile = "⁨‎diagnostic.R";
        
        PrepareTempFileInputStan(filePaths, programOptions, files, 0);
        SMC_State_PosetSMCParams smc(posetSMC, smcOptions);
        
        //for (int i = 0; i < num_chains; ++i)
        
        for( size_t iter = 0; iter < numTotalIter; ++iter){
            
            smc.run_smc(psParams);
            
            //get the best tree topology
            ParticlePopulation<State> *currenPop = smc.get_curr_population();
            vector<shared_ptr<State>> *particles = currenPop->get_particles();
            double   log_marginal = smc.get_log_marginal_likelihood();
            std::cout << "Marginal likelihood "<< " " <<log_marginal <<std::endl;
            
            std::vector<RootedTree> bestTrees;
            std::vector<long double> weights;
            std::vector<long double> rootLogLiks;
            
            std::vector<double> *normalized_weights = currenPop->get_normalized_weights();
            double max = -DOUBLE_INF;
            shared_ptr<State> best_particle;
            
            Eigen::ArrayXXd allCoalTimesScaledByTheta(smcOptions.num_particles,num_iter-1);
            for (size_t i= 0; i< smcOptions.num_particles; i++){
                shared_ptr<State> currents = particles->at(i);
                // std::cout << "Particle weight "<< i << " " <<(*normalized_weights)[i] <<std::endl;
                assert(currents->getRoots().size() == 1);
                weights.push_back((*normalized_weights)[i]);
                rootLogLiks.push_back(currents->getRootAt(0)->ln_likelihood);
                //std::cout << "Root loglik "<< i << " " <<currents->getRootAt(0)->ln_likelihood <<std::endl;
                std::vector<double> sortedCoalTimes=currents->getCoalEventTimesScaledBytheta();
                std::sort(sortedCoalTimes.begin(),sortedCoalTimes.end() );
                allCoalTimesScaledByTheta.row(i) = Eigen::Map<Eigen::RowVectorXd>(sortedCoalTimes.data(),num_iter-1 );
                //currents->printTreeChronologicalOrder(currents->getRoots()[0],std::cerr);
                if ((*normalized_weights)[i] > max) {
                    max = (*normalized_weights)[i]  ;
                    best_particle = currents;
                }
            }
            
            std::cout << "Tree topology with the greater loglik: " <<"\n" << endl;
            best_particle->printTree(best_particle->getRoots()[0], std::cerr);
            std::string newick = best_particle->getNewick(best_particle->getRootAt(0));
            bestTrees.emplace_back(RootedTree(newick, false));
            // if the RF distance between of best tree of current iteration and the best tree of the pevious iteration
            //is greater than 0 then run HMC
            
            PrepareTempFileInputStan(filePaths, programOptions, files, iter);
            
            //save info to a file files.fpStanDump->path
            fstream stan_input_stream (files.fpStanDump->path, std::fstream::out);
            
            //open the file files.fpStanDump->path
            string line;
            ifstream data_stream (files.fpStanDump->path, std::fstream::in);
            if (data_stream.is_open())
            {
                //           while ( getline (data_stream,line) )
                //           {
                //             cout << line << '\n';
                //           }
                
            }
            
            if (data_stream.rdstate() & std::ifstream::failbit) {
                std::stringstream msg;
                msg << "Can't open specified file, \"" << files.fpStanDump->path << "\"" << std::endl;
                data_stream.close();
                throw std::invalid_argument(msg.str());
            }
            stan::io::dump data_var_context(data_stream);
            
            data_stream.close();
            
            stan::callbacks::writer init_writer;
            stan::callbacks::interrupt interrupt;
            std::fstream output_stream(outputfile.c_str(),
                                       std::fstream::out);
            stan::callbacks::stream_writer sample_writer(output_stream, "# ");
            
            stan::callbacks::stream_writer info(std::cout);
            stan::callbacks::stream_writer err(std::cout);
            
            std::fstream diagnostic_stream(diagnosticfile.c_str(),
                                           std::fstream::out);
            stan::callbacks::stream_writer diagnostic_writer(diagnostic_stream, "# ");
            stan::io::empty_var_context context;
            //stan::io::dump metric_context(get_var_context(metric_file->value()));
            
            stan::callbacks::stream_logger logger(std::cout, std::cout, std::cout,std::cerr, std::cerr);
            stan_model model(data_var_context, random_seed, &std::cout);
            
            unsigned int chain = 1;
            double init_radius = 0;
            int num_warmup = numIterCondictionalHMC / 2;
            int num_samples = numIterCondictionalHMC / 2;
            int num_thin = 1;
            bool save_warmup = true;
            int refresh = 0;
            double stepsize = 0.1;
            double stepsize_jitter = 0;
            int max_depth = 8;
            double delta = .1;
            double gamma = .1;
            double kappa = .1;
            double t0 = .1;
            unsigned int init_buffer = 50;
            unsigned int term_buffer = 50;
            unsigned int window = 100;
            int return_code;
            
            return_code = stan::services::sample::hmc_nuts_diag_e_adapt(
                                                                        model, context, random_seed, chain, init_radius, num_warmup, num_samples,
                                                                        num_thin, save_warmup, refresh, stepsize, stepsize_jitter, max_depth,
                                                                        delta, gamma, kappa, t0, init_buffer, term_buffer, window, interrupt,
                                                                        logger, init_writer, sample_writer, diagnostic_writer);
            
            //draw   psParams.theta
            //         psParams.deltas,
            //          psParams.timeOriginSTDs from conditional posterior
            
            int num_output_lines = (num_warmup + num_samples) / num_thin;
            std::vector<std::string> init_values;
//            init_values = init_writer.string_values();
//            std::vector<std::vector<std::string> > parameter_names;
//            parameter_names = sample_writer.;
//            std::vector<std::vector<double> > parameter_values;
//            parameter_values = sample_writer.vector_double_values();
            //        std::vector<std::vector<std::string> > diagnostic_names;
            //        diagnostic_names = diagnostic_writer.vector_string_values();
            //        std::vector<std::vector<double> > diagnostic_values;
            //        diagnostic_values = diagnostic_writer.vector_double_values();
            
            
        }
        
        
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

//stan::io::dump data_var_context(get_var_context(dynamic_cast<string_argument*>(parser.arg("data")->arg("file"))->value()));
//
//  std::fstream output_stream(dynamic_cast<string_argument*>(parser.arg("output")->arg("file"))->value().c_str(),
//                             std::fstream::out);
//  stan::callbacks::stream_writer sample_writer(output_stream, "# ");
//
//  std::fstream diagnostic_stream(dynamic_cast<string_argument*>(parser.arg("output")->arg("diagnostic_file"))->value().c_str(),
//                                 std::fstream::out);
//  stan::callbacks::stream_writer diagnostic_writer(diagnostic_stream, "# ");
//
//  unsigned int random_seed = dynamic_cast<u_int_argument*>(parser.arg("random")->arg("seed"))->value();




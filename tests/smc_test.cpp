#include  "smc_test.hpp"
#include <iostream>
#include <string>

#include <map>
#include <unordered_map>
#include <string>
#include <omp.h>


#include "gtest/gtest.h"


TEST_P(PosetSMCOnePopulationTest,  smcTest4tipsDelta100Theta1)
{
    
    std::string method =  posetSMC->getPosetSMCKernelName();
    std::string healthyTipLabel = "healthycell";
    std::vector<double> coalTimes = pll_utils::getOrderedCoalTimesFromRootedTree(initialRootedTree, healthyTipLabel);
    
    double theta = true_theta;
    psParams->trueCoalTimes = coalTimes;
    
    std::cout<< "\n\nRunning Sequential Monte Carlo(SMC)" << " with "<<smcOptions.num_particles <<" particles and "<< method << "....\n" << std::endl;
    smc->run_smc(*psParams);
    
    ParticlePopulation<State> *currenPop = smc->get_curr_population();
    vector<shared_ptr<State>> *particles = currenPop->get_particles();
    double   log_marginal = smc->get_log_marginal_likelihood();
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
    
    Eigen::ArrayXXd allCoalTimesScaledByTheta(smcOptions.num_particles,num_iter-1);
    
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
        std::vector<double> sortedCoalTimes=currents->getCoalEventTimesScaledBytheta();
        std::sort(sortedCoalTimes.begin(),sortedCoalTimes.end() );
        allCoalTimesScaledByTheta.row(i) = Eigen::Map<Eigen::RowVectorXd>(sortedCoalTimes.data(),num_iter-1 );
        //currents->printTreeChronologicalOrder(currents->getRoots()[0],std::cerr);
        
        if ((*normalized_weights)[i] > max) {
            
            max = (*normalized_weights)[i]  ;
            best_particle = currents;
        }
        std::string newick = currents->getNewick(currents->getRootAt(0));
        lastPopulationTrees.emplace_back(RootedTree(newick, false));
    }
    lastPopulationTrees.emplace_back(*rootedTree);
    
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
    
    std::cout << "Avg RF(in the last SMC population + true tree) " << avgRF  << std::endl;
    std::cout << "Number unique topologies(counting the true topology)" << numUniqueTopologies  << std::endl;
    std::cout << "Number unique topologies(in the last SMC population)" << numUniqueTopologies-1  << std::endl;
    
    Eigen::Array<double, 1, Eigen::Dynamic> means = allCoalTimesScaledByTheta.colwise().mean();
    Eigen::Array<double, 1, Eigen::Dynamic> mins = allCoalTimesScaledByTheta.colwise().minCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> maxs = allCoalTimesScaledByTheta.colwise().maxCoeff();
    Eigen::Array<double, 1, Eigen::Dynamic> std_dev = ((allCoalTimesScaledByTheta.rowwise() - means.array()).square().colwise().sum()/(num_iter-2)).sqrt();
    
    
    transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
    for (size_t i=0; i<(num_iter-1); i++){
        
        std::cout << "coal time scaled by theta " << i << " min, max: "<< mins(i)<< ","<< maxs(i) << " mean: "<< means(i) <<" std: " << std_dev(i) << " true "<< coalTimes[i]*true_theta<< "\n"<<std::endl;
        
    }
    
    
    double log_marginal_lik = smc->get_log_marginal_likelihood();
    std::cout << "Estimate log marginal " << log_marginal  << endl;
    std::cout << "Estimate log P(y)= " << log_marginal_lik  << endl;
    
    assert(best_particle->getRoots().size() == 1);
    std::cout << "Tree topology with the greater loglik: "  << endl;
    best_particle->printTree(best_particle->getRoots()[0], std::cerr);
    
    
    //EXPECT_EQ(avgRF,0.0);
    EXPECT_TRUE(avgRF<1);
    //assert(abs(avgRF-0) <0.0000000001);
    assert(numUniqueTopologies >=1 );
    
}
INSTANTIATE_TEST_SUITE_P(
                         PosetSMCOnePopulationTests,
                         PosetSMCOnePopulationTest,
                         ::testing::Values(
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0009_0001_n4.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0009_n4.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0006_0001_n15_theta1.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0006_n15_theta1.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001_n15.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n15.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0001_0001.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Model_time_0002_0001.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false )
                                           ));
//the order of the parameters in the tuple is inputPhylipPath,
//inputPhylipTree, num_partcicles, resampling_scheme, smc_method,
//true_theta, true_delta

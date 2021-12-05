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
    string usePriorDensityasProposal = ( psParams->usePriorInSMC1)?"with prior proposal": "with post proposal";
    string PriorPostIncrementDensity = posetSMC->getPosetSMCPriorName();
    
    
    std::cout<< "\n\nRunning Sequential Monte Carlo(SMC)" << " with "<<smcOptions.num_particles <<" particles and "<< method << " "<< usePriorDensityasProposal << " " <<"...." << std::endl;
    if (posetSMC->kernelType== PosetSMC::PosetSMCKernel::PRIORPOST)
        std::cout<<"The proposal density for the increments is " << PriorPostIncrementDensity << std::endl;
    if (programOptions.normalizeClv){
        std::cout<<"The clv vectors will be normalized!(inner product clv vector and stable proportion vector equals 1)" <<  std::endl;

    }
    
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
    if (rootedTree)
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
    
    double log_marginal_lik = smc->get_log_marginal_likelihood();
    std::cout << "Estimate log marginal " << log_marginal  << endl;
    std::cout << "Estimate log P(y)= " << log_marginal_lik  << endl;
    
    assert(best_particle->getRoots().size() == 1);
    std::cout << "Tree topology with the greater loglik: " <<"\n" << endl;
    best_particle->printTree(best_particle->getRoots()[0], std::cerr);
    
    if (coalTimes.size()>0){
        transform(coalTimes.begin(), coalTimes.end(), coalTimes.begin(), [theta](double &c){ return c/theta; });
        for (size_t i=0; i<(num_iter-1); i++){
         if (i==0)
            std::cout <<"\n"<< endl;
        std::cout << "coal time scaled by theta " << i << " min, max: "<< mins(i)<< ","<< maxs(i) << " mean: "<< means(i) <<" std: " << std_dev(i) << " true "<< coalTimes[i]*true_theta<< "\n"<<std::endl;
        
        }
    }
    

    //EXPECT_EQ(avgRF,0.0);
    EXPECT_TRUE(avgRF<1);
    assert(numUniqueTopologies >=1 );
    
}
INSTANTIATE_TEST_SUITE_P(
                         PosetSMCOnePopulationTests,
                         PosetSMCOnePopulationTest,
                         ::testing::Values(
                                          std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0009_0001_n4.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0009_n4.tre", 1000, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0,0.057393051,  false, true ),
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_haplotypes_dir/true_hap_Delta=10.618_T=0.309_theta=0.034_nMU=3.420_n=72_0002_0012_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=10.618_T=0.309_theta=0.034_nMU=3.420_n=72_0002_0012_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.034,10.618,0.309, false,  false ),
                                               
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_haplotypes_dir/true_hap_Delta=10.804_T=0.255_theta=0.057_nMU=3.321_n=86_0002_0056_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=10.804_T=0.255_theta=0.057_nMU=3.321_n=86_0002_0056_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.057,10.804,0.255, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_haplotypes_dir/true_hap_Delta=104.770_T=0.054_theta=0.694_nMU=0.529_n=18_0002_0085_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=104.770_T=0.054_theta=0.694_nMU=0.529_n=18_0002_0085_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.694,104.770,0.054, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_haplotypes_dir/true_hap_Delta=107.926_T=0.046_theta=0.961_nMU=0.704_n=31_0002_0045_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=107.926_T=0.046_theta=0.961_nMU=0.704_n=31_0002_0045_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.961,107.926, 0.046,false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=110.864_T=0.042_theta=1.672_nMU=0.536_n=19_0002_0087_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=110.864_T=0.042_theta=1.672_nMU=0.536_n=19_0002_0087_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.672, 110.864,0.042, false,  false ),
                                           
                                          std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=12.393_T=0.334_theta=0.166_nMU=3.089_n=84_0002_0011_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=12.393_T=0.334_theta=0.166_nMU=3.089_n=84_0002_0011_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.166, 12.393,0.334, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=12.841_T=0.273_theta=1.107_nMU=3.225_n=73_0002_0044_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=12.841_T=0.273_theta=1.107_nMU=3.225_n=73_0002_0044_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.107, 12.841,0.273, false,  false ),
                                           
                                           
                                            std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=124.018_T=0.035_theta=1.442_nMU=0.784_n=52_0002_0092_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=124.018_T=0.035_theta=1.442_nMU=0.784_n=52_0002_0092_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.442, 124.018,0.035, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=130.106_T=0.037_theta=0.069_nMU=1.078_n=77_0002_0066_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=130.106_T=0.037_theta=0.069_nMU=1.078_n=77_0002_0066_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.069, 130.106, 0.037,false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=131.724_T=0.086_theta=0.335_nMU=0.810_n=39_0002_0008_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=131.724_T=0.086_theta=0.335_nMU=0.810_n=39_0002_0008_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.335,131.724,0.086, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=135.606_T=0.060_theta=1.514_nMU=0.563_n=22_0002_0033_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=135.606_T=0.060_theta=1.514_nMU=0.563_n=22_0002_0033_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.514,135.606, 0.060, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=138.055_T=0.047_theta=2.911_nMU=0.609_n=32_0002_0088_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=138.055_T=0.047_theta=2.911_nMU=0.609_n=32_0002_0088_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.911,138.055,0.047, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=136.037_T=0.051_theta=4.541_nMU=1.104_n=84_0002_0006_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=136.037_T=0.051_theta=4.541_nMU=1.104_n=84_0002_0006_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 4.541,136.037,0.051, false,  false ),
                                           
                    
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=139.706_T=0.032_theta=1.174_nMU=1.187_n=92_0002_0021_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=139.706_T=0.032_theta=1.174_nMU=1.187_n=92_0002_0021_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.174,139.706, 0.032,false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=143.540_T=0.031_theta=0.085_nMU=0.940_n=77_0002_0067_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=143.540_T=0.031_theta=0.085_nMU=0.940_n=77_0002_0067_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.085,143.540, 0.031,false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=143.569_T=0.035_theta=0.238_nMU=0.991_n=79_0002_0091_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=143.569_T=0.035_theta=0.238_nMU=0.991_n=79_0002_0091_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.238,143.569,0.035, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=146.162_T=0.032_theta=1.734_nMU=0.835_n=66_0002_0055_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=146.162_T=0.032_theta=1.734_nMU=0.835_n=66_0002_0055_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.734,146.162,0.032, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=146.397_T=0.034_theta=1.508_nMU=0.505_n=32_0002_0004_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=146.397_T=0.034_theta=1.508_nMU=0.505_n=32_0002_0004_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.508,146.397,0.034, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=156.620_T=0.024_theta=0.147_nMU=0.914_n=91_0002_0079_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=156.620_T=0.024_theta=0.147_nMU=0.914_n=91_0002_0079_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.147,156.620, 0.024, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=16.119_T=0.294_theta=2.006_nMU=2.799_n=51_0002_0075_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=16.119_T=0.294_theta=2.006_nMU=2.799_n=51_0002_0075_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.006,16.119,0.294, false,  false ),
                                           
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=162.411_T=0.029_theta=1.050_nMU=0.742_n=53_0002_0038_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=162.411_T=0.029_theta=1.050_nMU=0.742_n=53_0002_0038_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.050,162.411,0.029, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=171.666_T=0.044_theta=0.173_nMU=1.049_n=97_0002_0093_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=171.666_T=0.044_theta=0.173_nMU=1.049_n=97_0002_0093_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.173,171.666, 0.044, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=172.294_T=0.026_theta=1.192_nMU=0.873_n=83_0002_0003_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=172.294_T=0.026_theta=1.192_nMU=0.873_n=83_0002_0003_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.192,172.294,0.026, false,  false ),
                                           
                                          std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=174.658_T=0.027_theta=2.991_nMU=0.368_n=19_0002_0047_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=174.658_T=0.027_theta=2.991_nMU=0.368_n=19_0002_0047_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.991,174.658,0.027, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=178.238_T=0.052_theta=0.005_nMU=0.948_n=70_0002_0090_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=178.238_T=0.052_theta=0.005_nMU=0.948_n=70_0002_0090_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::PRIORPOST, 0.005,178.238,0.052, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=183.603_T=0.034_theta=0.737_nMU=0.951_n=94_0002_0073_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=183.603_T=0.034_theta=0.737_nMU=0.951_n=94_0002_0073_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.737,183.603,0.034, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=186.718_T=0.026_theta=0.702_nMU=0.654_n=48_0002_0022_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=186.718_T=0.026_theta=0.702_nMU=0.654_n=48_0002_0022_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.702,186.718, 0.026, false,  false ),
                                           
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=188.507_T=0.026_theta=1.724_nMU=0.686_n=53_0002_0071_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=188.507_T=0.026_theta=1.724_nMU=0.686_n=53_0002_0071_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.724,188.507, 0.026, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=19.179_T=0.118_theta=0.269_nMU=2.579_n=92_0002_0095_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=19.179_T=0.118_theta=0.269_nMU=2.579_n=92_0002_0095_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.269,19.179, 0.118, false,  false ),
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=19.204_T=0.111_theta=0.577_nMU=2.128_n=65_0002_0086_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=19.204_T=0.111_theta=0.577_nMU=2.128_n=65_0002_0086_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 0.577,19.204,0.111, false,  false ),
                                           
                                           
                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/true_hap_Delta=392.693_T=0.017_theta=4.901_nMU=0.589_n=96_0002_0094_0001.txt", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/Data_benchmark2/Cellphy_ML_trees/Delta=392.693_T=0.017_theta=4.901_nMU=0.589_n=96_0002_0094_0001.phy.raxml.bestTree", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 4.901,392.693, 0.017, false,  false )
//                                             std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0009_0001_n4.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0009_n4.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0006_0001_n15_theta1.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0006_n15_theta1.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001_n15.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n15.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0001_0001.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Model_time_0002_0001.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001._n4_delta10.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n4_delta10.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,10.0, false ),
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0009_0001_n4.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0009_n4.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0006_0001_n15_theta1.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0006_n15_theta1.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001_n15.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n15.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0001_0001.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Model_time_0002_0001.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false )
//                                          std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001._n4_delta10.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n4_delta10.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,10.0, false ),
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0009_0001_n4.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0009_n4.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0006_0001_n15_theta1.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0006_n15_theta1.tre", 256, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 1.0,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0010_0001_n15.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Scaled_physical_time_0002_0010_n15.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false ),
//
//                                           std::make_tuple("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/true_hap_0002_0001_0001.phylip", "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/data/input_data/trees_Model_time_0002_0001.tre", 128, SMCOptions::ResamplingScheme::STRATIFIED, PosetSMC::PosetSMCKernel::TSMC1, 2.5,100.0, false )
                                           ));
//the order of the parameters in the tuple is inputPhylipPath,
//inputPhylipTree, num_partcicles, resampling_scheme, smc_method,
//true_theta, true_delta, true_T doPlots, isUltrametricTree

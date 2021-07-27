//
//  mcmc.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 15/03/2021.
//

#include "mcmc.hpp"

MCMC::MCMC(int  numChains){
    
    if (numChains >0)
        chainManager= new ChainManager(numChains);
    
}
void MCMC::initialize(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths,
                      std::vector<int> &sampleSizes,
                      std::vector<std::vector<int> > &ObservedData,
                      char* ObservedCellNames[],
                      pll_msa_t *msa,
                      pll_rtree_t * initialRootedTree,
                      std::vector<StructuredCoalescentTree *> structuredCoalTrees,
                      std::string& healthyTipLabel,
                      const std::vector<pll_rtree_t *> &trueTrees,
                      const  std::vector<long double> &trueThetas,
                      const std::vector<std::vector<long double>> &trueDeltaTs,
                      const  std::vector<std::vector<long double>> &trueTs){
    
    chainManager->initializeChains(randomGenerators, programOptions, mcmcOptions, filePaths,
                                   sampleSizes, ObservedData, ObservedCellNames, msa, initialRootedTree,
                                   structuredCoalTrees, healthyTipLabel, trueTrees, trueThetas, trueDeltaTs, trueTs);
    
}
void MCMC::runMCMC(std::vector< gsl_rng * > &randomGenerators, ProgramOptions &programOptions, MCMCoptions &mcmcOptions, FilePaths &filePaths){
    
    bool converged =false;
    std::cout << "\n Tuning MCMC parameters... \n"<<std::endl;
    chainManager->performWarmUp(programOptions, mcmcOptions, randomGenerators, filePaths);
    std::cout << "\n MCMC parameters tuned! \n"<<std::endl;
    chainManager->resetChainValues(programOptions, mcmcOptions);
    
    for (size_t currentIteration = 0; currentIteration < mcmcOptions.Niterations; ++currentIteration)
    {
        chainManager->stepAllChainsNoSave( currentIteration,programOptions, mcmcOptions, randomGenerators  );
        
        if (currentIteration >=  mcmcOptions.burnInIterations ){
            chainManager->saveChainState( currentIteration, programOptions, mcmcOptions );

        }
        chainManager->writeChainsOutput(currentIteration, programOptions, mcmcOptions, filePaths);
        if ( (currentIteration >= (mcmcOptions.burnInIterations +mcmcOptions.maxNumberIndependentPosteriorValues)) &&  currentIteration % 1000 ==0 ){
            converged =  chainManager->checkConvergence(mcmcOptions.numberChainsPerTree, programOptions, mcmcOptions, 0.01) ;
        }
        if (chainManager->numberFinishedChains >= 0.95*mcmcOptions.numChains){
            std::cout << "\nMCMC chains converged .... \n"<<std::endl;
            break;
        }
    }
    
}
//for (size_t currentIteration = 0; currentIteration < mcmcOptions.Niterations; ++currentIteration)
//         {
//             chainManager->stepAllChainsNoSave( currentIteration,programOptions, mcmcOptions, randomGenerators  );
//             
//            if (currentIteration >=  mcmcOptions.burnInIterations ){
//                   chainManager->saveChainState( currentIteration, programOptions, mcmcOptions );
//     //              //chainManager->saveChainState( currentIteration, programOptions, mcmcOptions );
//     //              chainManager->writeChainsOutput(currentIteration, programOptions, mcmcOptions, filePaths);
//             }
//              chainManager->writeChainsOutput(currentIteration, programOptions, mcmcOptions, filePaths);
//             if ( (currentIteration >= (mcmcOptions.burnInIterations +mcmcOptions.maxNumberIndependentPosteriorValues)) &&  currentIteration % 1000 ==0 ){
//                   chainManager->checkConvergence(mcmcOptions.numberChainsPerTree, programOptions, mcmcOptions, 0.01) ;
//             }
//             if (chainManager->numberFinishedChains >= 0.98*mcmcOptions.numChains)
//                 break;
//     //        if (converged){
//     //            std::cout << "\nMCMC chains converged .... \n"<<std::endl;
//     //            //break;
//     //           // mcmcOptions.Niterations = mcmcOptions.Niterations + 10000;
//     //        }
//         }
//

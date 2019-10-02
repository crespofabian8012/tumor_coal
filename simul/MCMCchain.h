//
//  MCMCchain.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 29/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef MCMCchain_h
#define MCMCchain_h
static void InitializeChains(Chain **chains,  const ProgramOptions *programOptions, const MCMCoptions *mcmcOptions, int * sampleSizes, long int *seed, char* ObservedCellNames[], pll_msa_t * msa);
static void InitChainPopulations(Population **populations, int numClones, int noisy,  int TotalNumSequences );
static void FillChainPopulationsFromPriors(Population **populations,Chain *chain, const ProgramOptions *programOptions, const MCMCoptions *mcmcOptions,int *sampleSizes, long int *seed);
static void runChain(Chain *chain,  const  MCMCoptions *opt,  long int *seed, const FilePaths *filePaths, Files *files, const ProgramOptions *programOptions,
              double  *varTimeGMRCA,char* ObservedCellNames[], pll_msa_t * msa, int *sampleSizes
              );
void setChainPopulationSampleSizes(Chain *chain, int *sampleSizes,  ProgramOptions *programOptions);
#endif /* MCMCchain_h */

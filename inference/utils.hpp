//
//  utils.hpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 5/10/19.
//

#ifndef utils_hpp
#define utils_hpp

#include "data_utils.hpp"
#include "Population.hpp"
#include "Chain.hpp"
#include "libpll/pll_optimize.h"
#include "libpll/pll_tree.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pll_msa.h"
#include "pllmod_common.h"
#include "definitions.h"

#include <stdarg.h>
#include <search.h>
#include <time.h>

void ReadParametersFromFastaFile(char *fileName, ProgramOptions *programOptions);

void ReadFastaFile(char *fileName, int** ObservedData,  char **ObservedCellNames, ProgramOptions *programOptions);

double DensityTime(double delta, double u);
double LogDensityCoalescentTimesForPopulation(pll_unode_t  *tree, pll_unode_t *nodes,  population **populations, int numClones);

double LogConditionalLikelihoodTree(pll_unode_t  *tree, pll_unode_t *nodes, population **populations, int numClones);

double LogProbNoCoalescentEventBetweenTimes(population *popI,double from, double to, int numberActiveInd);

double  LogConditionalLikelihoodSequences(pll_msa_t * msa, char* NewickString, ProgramOptions *programOptions, double seqError,double dropoutError);

void dealloc_data(pll_unode_t * node, void (*cb_destroy)(void *));

int set_tipclv1(pll_partition_t * partition,
                unsigned int tip_index,
                const pll_state_t * map,
                const char * sequence,
                double _seq_error_rate,
                double _dropout_rate);
void  destroyTree(pll_utree_t * tree, void (*cb_destroy)(void *));

double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);

void set_partition_tips( pll_partition_t * partition, pll_msa_t * msa, ProgramOptions *programOptions, double seqError, double dropoutError);



double RandomLogUniform( double from, double to, long int *seed);
void  RandomDirichlet (double s, int vectorSize, double **outputVector, long int *seed);

void InitPopulationSampleSizes(population **populations, int TotalSampleSize, int numClones, double *proportionsVector, long int *seed);
#endif /* utils_hpp */

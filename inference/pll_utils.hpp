//
//  pll_utils.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 15/03/2021.
//

#ifndef pll_utils_hpp
#define pll_utils_hpp

#include  "data_types.hpp"

extern "C"
    {
    #include "libpll/pll.h"
    #include "libpll/pll_optimize.h"
    #include "libpll/pllmod_common.h"
    #include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
    #include <libpll/pll_tree.h>
    #include <libpll/pllmod_util.h>
    #include "libpll/pll_optimize.h"
    #include <libpll/pllmod_common.h>
    #include <libpll/pllmod_algorithm.h>
}

namespace pll_utils{
double  * expandUniqRates(int states, const  double  * uniq_rates, const int * rate_sym);

long double  LogConditionalLikelihoodSequences(pll_msa_t * msa,
                                               char* NewickString,
                                               ProgramOptions &programOptions,
                                               long double seqError,
                                               long double dropoutError);
long double  LogConditionalLikelihoodSequences(pll_msa_t * msa,
                                               pll_utree_t *unrootedTree,
                                               ProgramOptions &programOptions,
                                               long double seqError,
                                               long double dropoutError);
long double  LogConditionalLikelihoodSequences(pll_msa_t * msa,
                                          pll_rtree_t *rootedTree,
                                          ProgramOptions &programOptions,
                                          long double seqError,
                                          long double dropoutError);
long double  LogConditionalLikelihoodSequencesRootedTree(pll_msa_t * msa,
                                                                    pll_rtree_t *rootedTree,
                                                                    ProgramOptions &programOptions,
                                                                    long double seqError,
                                                                    long double dropoutError);
void setPartitionTipsCostum( pll_partition_t * partition,
                               pll_msa_t * msa,
                               ProgramOptions &programOptions,
                               long double  seqError,
                               long double  dropoutError);


int setTipClvCustomErrorModel(pll_partition_t * partition,
                                  unsigned int tip_index,
                                  const pll_state_t * map,
                                  const char * sequence,
                                  long double  _seq_error_rate,
                                  long double  _dropout_rate
                                  );
void  destroyUTree(pll_utree_t * tree, void (*cb_destroy)(void *));
void  destroyRTree(pll_rtree_t * tree, void (*cb_destroy)(void *));
void deallocUDataCostum(pll_unode_t * node, void (*cb_destroy)(void *));
void deallocRDataCostum(pll_rnode_t * node, void (*cb_destroy)(void *));
pll_rnode_t * cloneRNode(const pll_rnode_t * node);
pll_rtree_t * cloneRTree(const pll_rtree_t * tree);
char * xstrdup(const char * s);
const pll_partition_t *createReferencePartition(pll_msa_t * msa) ;
const pll_partition_t *createGTReferencePartition(pll_msa_t * msa);
double computeLogLikelihood(double *clv, unsigned int *scale_buffer,
                                        const pll_partition_t *p);
int cb_rfull_traversal(pll_rnode_t * node);
double meanHeterozygosity(long seq_count, long site_count, pll_msa_t* msa);
double meanDistance(long seq_count, long site_count, pll_msa_t* msa);
void computeCoalTimesInsideNode(pll_rnode_t *node,std::string& healthyTipLabel, unsigned int &consecutiveIndex );
void printPostOrderRootedTree(pll_rtree_t *tree,std::string& healthyTipLabel, std::ostream &stream);
void printChronologicalOrderRootedTree(pll_rtree_t *tree,std::string& healthyTipLabel, std::ostream &stream);
bool compareByTime (pll_rnode_t* x,pll_rnode_t* y);
bool isLeaf(pll_rnode_t *node);
std::vector<double> getOrderedCoalTimesFromRootedTree(pll_rtree_t *rootedTree, std::string& healthyTipLabel);

}
#endif /* pll_utils_hpp */

//
//  treeLikelihood.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/06/2021.
//

#ifndef treeLikelihood_hpp
#define treeLikelihood_hpp

#include "pll_utils.hpp"
#include "partition.hpp"
#include "genotype_error_model.hpp"
#include <stdio.h>


class TreeLikelihood{
private:
    pll_rtree_t * rooted_tree;
    Partition* reference_partition;
    pll_msa_t * msa;
    GenotypeErrorModel *gtErrorModel;
    pll_operation_t * operations ;
    unsigned int ops_count;
    double *branch_lengths;
    unsigned int * matrix_indices;
    int number_genotype_states;
    int number_phased_genotype_states;
    int number_states;
public:
    TreeLikelihood(Partition *reference_partition, pll_rtree_t * rooted_tree, pll_msa_t * msa, GenotypeErrorModel *gtErrorModel);
    
    TreeLikelihood(int numberTips,
       int clvBuffers,
       int numberStates,
       int numberSites,
       int numberRateMatrices,
       int probMatrices,
       int numberRateCats,
       int numberScaleBuffers,
       int statesPadded,
       bool sse,
       bool avx,
       bool  avx2,
       bool avx512,
       bool  asc,
       bool tipPatternCompression, pll_rtree_t * rooted_tree, pll_msa_t * msa, GenotypeErrorModel *gtErrorModel);
   
    double computeRootLogLikelihood();
    int setPartitionTips(bool doUseGenotypes);
    int initOperations();
    void fillTipClv(unsigned int tip_id, std::vector<double>& clv) const;
    void recomputeTipClvs();
    void createHashTable();
    void destroyHashTable();
    void changeGenotypeErrorModel(GenotypeErrorModel *gtErrorModel);
    static void fillTipClv(unsigned int tip_id, std::vector<double>& clv, pll_msa_t *msaP, GenotypeErrorModel &gtError, unsigned int number_statesP);
    
    ~TreeLikelihood()
       {
           //pll_msa_destroy(msa);
           delete reference_partition;
          // pll_rtree_destroy (rooted_tree, NULL);
           
           free( branch_lengths);
           free( matrix_indices);
           free(operations);
           //destroyHashTable();
       }
};
#endif /* treeLikelihood_hpp */

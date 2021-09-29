//
//  partition.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 27/06/2021.
//

#include "partition.hpp"
#include "core_likelihood.hpp"
#include <memory.h>

Partition::Partition(unsigned int numberTips,
                     unsigned int clvBuffers,
                     unsigned int numberStates,
                     unsigned int numberSites,
                     unsigned int numberRateMatrices,
                     unsigned int probMatrices,
                     unsigned int numberRateCats,
                     unsigned int numberScaleBuffers,
                     unsigned int statesPadded,
                     bool sse=false,
                     bool avx=false,
                     bool  avx2=false,
                     bool avx512=false,
                     bool  asc=false,
                     bool tipPatternCompression=false):
numberTips(numberTips),
clvBuffers(clvBuffers),
numberStates(numberStates),
numberSites(numberSites),
numberRateMatrices(numberRateMatrices),
probMatrices(probMatrices),
numberRateCats(numberRateCats),
numberScaleBuffers(numberScaleBuffers),
statesPadded(statesPadded),
sse(sse),
avx(avx),
avx2(avx2),
avx512(avx512),
asc(asc),
tipPatternCompression(tipPatternCompression)
{
    
    if (sse){
        attributes |= PLL_ATTRIB_ARCH_SSE ;
    }
    else if (avx)
        attributes |= PLL_ATTRIB_ARCH_AVX;
    else if (avx2)
        attributes |= PLL_ATTRIB_ARCH_AVX2;
    else if (avx512)
        attributes |= PLL_ATTRIB_ARCH_AVX512;
    else
        attributes |= PLL_ATTRIB_ARCH_CPU;
    
    if (asc)
        attributes |= PLL_ATTRIB_AB_LEWIS;
    
    unsigned int attr;
    if (sse)
        attr= PLL_ALIGNMENT_SSE;
    else if (avx || avx2)
        attr=PLL_ALIGNMENT_AVX;
    else
        attr=PLL_ALIGNMENT_CPU;
    
    attributes &= ~PLL_ATTRIB_SITE_REPEATS;
    attributes &= ~PLL_ATTRIB_PATTERN_TIP;
    
    partition =  pll_partition_create(numberTips, clvBuffers, numberStates, numberSites, numberRateMatrices, probMatrices, numberRateCats, numberScaleBuffers, attributes);
    
    if(!partition)
    {
        printf("Fail creating partition \n");
        //pll_partition_destroy(partition);
        exit(pll_errno);
        
    }
    model = pllmod_util_model_info_genotype("GT16JC");
    sumtable = (double *)pll_aligned_alloc(partition->sites * partition->rate_cats * partition->states_padded * sizeof(double), partition->alignment);
    
    if (!sumtable)
    {
        printf("Fail creating sumtable \n");
        //pll_partition_destroy(partition);
        exit(pll_errno);
    }
}


Partition::Partition(pll_msa_t *msa,
                     int numberStates,
                     int numberRateCats,
                     int statesPadded,
                     bool sse,
                     bool avx,
                     bool  avx2,
                     bool avx512,
                     bool  asc,
                     bool tipPatternCompression):
numberStates(numberStates),
numberRateCats(numberRateCats),
statesPadded(statesPadded),
sse(sse),
avx(avx),
avx2(avx2),
avx512(avx512),
asc(asc),
tipPatternCompression(tipPatternCompression)
{
    
    if (sse){
        attributes |= PLL_ATTRIB_ARCH_SSE ;
    }
    else if (avx)
        attributes |= PLL_ATTRIB_ARCH_AVX;
    else if (avx2)
        attributes |= PLL_ATTRIB_ARCH_AVX2;
    else if (avx512)
        attributes |= PLL_ATTRIB_ARCH_AVX512;
    else
        attributes |= PLL_ATTRIB_ARCH_CPU;
    
    if (asc)
        attributes |= PLL_ATTRIB_AB_LEWIS;
    
    unsigned int attr;
    if (sse)
        attr= PLL_ALIGNMENT_SSE;
    else if (avx || avx2)
        attr=PLL_ALIGNMENT_AVX;
    else
        attr=PLL_ALIGNMENT_CPU;
    
    attributes &= ~PLL_ATTRIB_SITE_REPEATS;
    attributes &= ~PLL_ATTRIB_PATTERN_TIP;
    
    numberTips=msa->count;
    numberSites= msa->length;
    
    
    //const unsigned int rate_category_count = 1;
    //unsigned int number_states=16;
    const unsigned int subst_model_count = 1;
    
    //const double subst_params[6] = {1, 1, 1, 1, 1, 1};
    //const unsigned int rate_category_count = 4;
    // double rate_categories[4] = {0, 0, 0, 0};
    double rate_categories[1] = {0};
    
    //model =  Model(DataType::autodetect, "GT16JC");

    model = pllmod_util_model_info_genotype("GT16JC");
 
    //assert(model->states ==numberStates);
    
    // double  rate_categories[number_states];
    // std::fill_n (rate_categories, number_states*(number_states-1) /2, 1);
    
    
    //double rate_categories[1] = {1};
    pll_compute_gamma_cats(1, 1, rate_categories, PLL_GAMMA_RATES_MEAN);
    //pll_compute_gamma_cats(1, 4, rate_categories, PLL_GAMMA_RATES_MEAN);
    
    //numberRateCats=4;
    numberRateCats=1;
    partition = pll_partition_create(
                                     msa->count,
                                     0, // Don't allocate any inner CLV's.
                                     numberStates, msa->length, subst_model_count,
                                     0, // Don't allocate any pmatrices.
                                     numberRateCats,
                                     0, // Don't allocate any scale buffers.
                                     PLL_ATTRIB_ARCH_SSE);
    
    if(!partition)
    {
        printf("Fail creating partition \n");
        //pll_partition_destroy(partition);
        exit(pll_errno);
        
    }
    
    
    double subst_params[ numberStates*(numberStates-1) /2] ;
    std::fill_n (subst_params,numberStates*(numberStates-1) /2, 1);
    
    
    
    double  category_weights[numberStates];
    std::fill_n (category_weights, numberStates*(numberStates-1) /2, 1);
    
    double genotype_frequencies[numberStates];
    std::fill_n (genotype_frequencies, numberStates, 1.0/(numberStates));
    
    
    pll_set_frequencies(partition, 0, genotype_frequencies);
    pll_set_category_rates(partition, rate_categories);
    pll_set_subst_params(partition, 0, subst_params);
    
    //   for (unsigned int i = 0; i < msa->count; i++) {
    //
    ////     if (pll_set_tip_states(partition, i, pll_map_gt10, msa->sequence[i]) != PLL_SUCCESS)
    ////       std::vector<double> tmp_clv;
    ////
    ////        computeCLV(i, msa, gtError, tmp_clv);
    ////
    ////            reference_partition->initTipCLV(i, tmp_clv.data());
    //       if (pll_set_tip_states(partition, i, pll_map_gt10, msa->sequence[i]) != PLL_SUCCESS)
    //        {
    //          fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
    //          exit(pll_errno);
    //        }
    //    }
    //
    pll_update_eigen(partition, 0);
    
    
    sumtable = (double *)pll_aligned_alloc(partition->sites * partition->rate_cats * partition->states_padded * sizeof(double), partition->alignment);
    
    if (!sumtable)
    {
        printf("Fail creating sumtable \n");
        //pll_partition_destroy(partition);
        exit(pll_errno);
    }
}


int  Partition::setTipStates(unsigned int tipClvIndex,std::string  sequence){
    
    if (pll_set_tip_states(partition, tipClvIndex, pll_map_gt10, sequence.c_str()) != PLL_SUCCESS)
    {
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
    }
    
    return PLL_SUCCESS;
}
int  Partition::setTipStates(unsigned int tipClvIndex, char * sequence){
    
    
    if (pll_set_tip_states(partition, tipClvIndex, pll_map_gt10, sequence) != PLL_SUCCESS)
    {
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
    }
    
    return PLL_SUCCESS;
}
int Partition::setTipCLV(unsigned int tipClvIndex, double * clv){
    
    if (pll_set_tip_clv(partition, tipClvIndex,  clv, PLL_TRUE) != PLL_SUCCESS)
    {
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
    }
    
    return PLL_SUCCESS;
    
}
int Partition::initTipCLV(unsigned int tipClvIndex, double * clv)const{
    
    if (pll_set_tip_clv(partition, tipClvIndex,  clv, PLL_TRUE) != PLL_SUCCESS)
    {
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
    }
    
    return PLL_SUCCESS;
    
}
void Partition::buildCLV(int tip_id, pll_msa_t *msa, GenotypeErrorModel *gtErrorModel, std::vector<double> &clv, bool normalize)const{
    
   // auto clv_size = msa->length * numberStates;
    //auto clv_size = msa->length * numberStates*numberRateCats;
    //here we assume numberRateCats=1
    auto clvp = clv.begin();
    
    auto seq = msa->sequence[tip_id];
    //auto charmap = _model.charmap();
    auto charmap = pll_map_gt10;
    
    for (size_t j = 0; j < msa->length; ++j)
    {
        auto charstate = (pll_state_t) seq[j];
        pll_state_t state = charmap ? charmap[(int) charstate] : charstate;
        
        
        double sum_lh = gtErrorModel->computeStateErrorProbPT20(state, clvp);
        
        if (j == 0 && 0)
        {
            printf("state: %llu ", state);
            printf("char: %c ", seq[j]);
            for (size_t k = 0; k < numberStates; ++k)
                printf("%lf ", clvp[k]);
            printf("\n");
            for (size_t k = 0; k < numberStates; ++k)
                printf("%lf ", clv[k]);
            printf("\n");
        }
        
        clvp += numberStates;
        //  clvp += numberStates*numberRateCats;
    }
    
    assert(clvp == clv.end());
}

void Partition::setPatternWeights( const unsigned int * patternWeights){
    
    pll_set_pattern_weights(partition, patternWeights);
    
}

void Partition::setSubstParams( unsigned int paramsIndex,
                               const double * params){
    
    pll_set_subst_params(partition,  paramsIndex,params);
    
}

void Partition::setFrequencies( unsigned int paramsIndex,
                               const double *frequencies){
    
    pll_set_frequencies(partition,  paramsIndex ,  frequencies);
    
}

void Partition::setCategoryRates(const double *rates){
    
    pll_set_category_rates(partition,  rates);
    
}

void Partition::setCategoryWeights(const double *rateWeights){
    
    pll_set_category_weights(partition,  rateWeights);
    
}

void Partition::updateEigen(unsigned int paramsIndex){
    
    pll_update_eigen(partition,  paramsIndex);
    
}

void Partition::updateProbMatrices(
                                   const  unsigned int * paramsIndexes,
                                   const  unsigned int * matrixIndexes,
                                   const double * branchLengths,
                                   unsigned int count){
    
    
    pll_update_prob_matrices(partition,
                             paramsIndexes,
                             matrixIndexes,
                             branchLengths,
                             count);
    
}

void Partition::updateInvariantSites(){
    pll_update_invariant_sites(partition);
}

void Partition::updateInvariantSitesProportion(int paramsIndex, double propInvar){
    pll_update_invariant_sites_proportion(partition,paramsIndex, propInvar);
}

void Partition::updatePartials(  pll_operation_t * operations,
                               unsigned int count){
    pll_update_partials(partition,operations, count);
}

//    def updatePartials(operations: Operations) = pll_update_partials(self, operations.pll_operation, operations.count)
//
double   Partition::computeRootLogLikelihood(int clvIndex, int  scalerIndex, unsigned int * freqsIndex, double * persiteLnL)
{
    
    double result= pll_compute_root_loglikelihood(partition,  clvIndex, scalerIndex, freqsIndex, persiteLnL);
    return result;
}
//
double  Partition::computeEdgeLogLikelihood(int parentCLVIndex, int parentScalerIndex, int childCLVIndex, int childScalerIndex, int matrixIndex, unsigned int* freqsIndex,double * persiteLnL){
    
    double result= pll_compute_edge_loglikelihood(partition, parentCLVIndex, parentScalerIndex, childCLVIndex, childScalerIndex, matrixIndex, freqsIndex, persiteLnL);
    return result;
    
}

int  Partition::updateSumtable(unsigned int parentCLVIndex,unsigned int childCLVIndex,int parentScalerIndex, int childScalerIndex, const unsigned int *paramsIndexes)
{
    
    if (pll_update_sumtable(partition, parentCLVIndex, childCLVIndex,parentScalerIndex, childScalerIndex, paramsIndexes,sumtable) != PLL_SUCCESS)
    {
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
    }
    
    return PLL_SUCCESS;
}


int  Partition::computeLikelihoodDerivatives(int parentScalerIndex, int childScalerIndex, double branchLength, const unsigned int * paramsIndices, double * df, double * ddf){
    
    
    if (pll_compute_likelihood_derivatives(partition, parentScalerIndex, childScalerIndex, branchLength, paramsIndices, sumtable, df, ddf) != PLL_SUCCESS)
    {
        fprintf(stderr, "PLL error %d: %s\n", pll_errno, pll_errmsg);
        exit(pll_errno);
    }
    
    return PLL_SUCCESS;
    
}

pll_partition_t * Partition::getPartition()const {
    
    return partition;
    
}
double* Partition::getCLV(int tipIndex)const {
    
    return(partition->clv[tipIndex]);
    
}
void Partition::showEigenDecomp(unsigned int float_precision) const{
    
    unsigned int i,j,k;
    double * ev;
    
    
    /* eigenvals */
    printf("Eigenvalues\n");
    printf("[");
    for (i = 0; i < partition->rate_matrices; ++i)
    {
        printf("(");
        for (j = 0; j < partition->states-1; ++j)
        {
            printf("%.*f,", float_precision, partition->eigenvals[i][j]);
        }
        printf("%.*f,", float_precision, partition->eigenvals[i][j]);
        printf(")");
        if (i < partition->rate_matrices - 1) printf(",");
    }
    printf("]\n");
    
    /* eigen vectors */
    printf("Eigenvectors\n");
    for (k = 0; k < partition->rate_matrices; ++k)
    {
        ev = partition->eigenvecs[k] + k*partition->states*partition->states;
        for (i = 0; i < partition->states; ++i)
        {
            for (j = 0; j < partition->states; ++j)
                printf("%+2.*f   ", float_precision, ev[i*partition->states+j]);
            printf("\n");
        }
        printf("\n");
    }
    
    printf("Inverse eigenvectors\n");
    for (k = 0; k < partition->rate_matrices; ++k)
    {
        ev = partition->inv_eigenvecs[k] + k*partition->states*partition->states;
        for (i = 0; i < partition->states; ++i)
        {
            for (j = 0; j < partition->states; ++j)
                printf("%+2.*f   ", float_precision, ev[i*partition->states+j]);
            printf("\n");
        }
        printf("\n");
    }
    
    
}

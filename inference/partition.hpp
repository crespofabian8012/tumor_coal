//
//  partition.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 27/06/2021.
//

#ifndef partition_hpp
#define partition_hpp

#include "genotype_error_model.hpp"
#include "model.hpp"

#include <stdio.h>
#include <string>
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

auto deletePartition = [](pll_partition_t* p)
  {
    pll_partition_destroy(p);

  };

class Partition{
private:
    //pll_partition_t * partition_rawPtr;
    std::shared_ptr<pll_partition_t> sharedptr_partition;
    std::unique_ptr<pll_partition_t, decltype(&pll_partition_destroy)>  partition;
    pllmod_subst_model_t* model;
    //Partition(){partition = nullptr;};
    
public:
    unsigned int attributes;
    unsigned int numberTips;
    unsigned int clvBuffers;
    unsigned int numberStates;
    unsigned int numberSites;
    unsigned int numberRateMatrices;
    unsigned int probMatrices;
    unsigned int numberRateCats;
    unsigned int numberScaleBuffers;
    unsigned int statesPadded;
    bool sse;
    bool avx;
    bool  avx2;
    bool avx512;
    bool asc;
    bool tipPatternCompression;
    double * sumtable;
    
    Partition(unsigned int numberTips,
              unsigned int clvBuffers,
              unsigned int numberStates,
              unsigned int numberSites,
              unsigned int numberRateMatrices,
              unsigned int probMatrices,
              unsigned int numberRateCats,
              unsigned int numberScaleBuffers,
              unsigned int statesPadded,
              bool sse,
              bool avx,
              bool  avx2,
              bool avx512,
              bool  asc,
              bool tipPatternCompression);
    /* contructor asociated with a msa*/
    Partition(pll_msa_t *msa,
                         int numberStates,
                         int numberRateCats,
                         int statesPadded,
                         bool sse=false,
                         bool avx=false,
                         bool  avx2=false,
                         bool avx512=false,
                         bool  asc=false,
              bool tipPatternCompression=false);
  
    
    Partition(const Partition &other)  = delete;
    void operator=(Partition const&)          = delete;
    // Partition(const Partition& original) ;
    //  Partition& operator=(const Partition &original) ;
    
    int  setTipStates(unsigned int tipClvIndex, std::string  sequence);
    int  setTipStates(unsigned int tipClvIndex, char * sequence);
    int setTipCLV(unsigned int tipClvIndex, double * clv);
    void setPatternWeights(const unsigned int * patternWeights);
    
    void setSubstParams( unsigned int paramsIndex,
                        const double * params);
    
    void setFrequencies(unsigned int paramsIndex, const double *frequencies);
    
    void setCategoryRates(const double *rates);
    
    void setCategoryWeights(const double *rateWeights);
    
    void updateEigen(unsigned int paramsIndex);
    
    void updateProbMatrices(const unsigned int * paramsIndexes,
                            const unsigned int * matrixIndexes,
                            const double * branchLengths,
                            unsigned int count);

    
    void updateInvariantSites();
    
    void updateInvariantSitesProportion(int paramsIndex, double propInvar);
    
    void updatePartials(pll_operation_t * operations,
                        unsigned int count);
    
    double   computeRootLogLikelihood(int clvIndex, int  scalerIndex, unsigned int * freqsIndex, double * persiteLnL);
    
    double  computeEdgeLogLikelihood(int parentCLVIndex, int parentScalerIndex, int childCLVIndex, int childScalerIndex, int matrixIndex, unsigned int* freqsIndex,double * persiteLnL);
    
    int  updateSumtable(unsigned int parentCLVIndex,unsigned int childCLVIndex,int parentScalerIndex, int childScalerIndex, const unsigned int *paramsIndexes);
    
    
    int  computeLikelihoodDerivatives(int parentScalerIndex, int childScalerIndex, double branchLength, const unsigned int * paramsIndices, double * df, double * ddf);
    
    pll_partition_t * getPartition()const ;
    
    int ascBiasCorrection() const{return sharedptr_partition->asc_bias_alloc; };
    
    unsigned int getStatesPadded() const{return sharedptr_partition->states_padded; };
    
    double* getCLV(int tipIndex) const;
    
    double* rates() const{return sharedptr_partition->rates; };
    
    double* propInvar() const{return sharedptr_partition->prop_invar; };
    
    double** eigenVals() const{return sharedptr_partition->eigenvals; };
    
    double** eigenVecs() const{return sharedptr_partition->eigenvecs; };
    
    double** invEigenVecs() const{return sharedptr_partition->inv_eigenvecs; };
    
    double** frequencies() const{return sharedptr_partition->frequencies; };
    
    double* rateWeights() const{return sharedptr_partition->rate_weights; };
    
    unsigned int* patternWeights() const{return sharedptr_partition->pattern_weights; };
    
    size_t alignment() const{return sharedptr_partition->alignment; };
    
    int* invariant() const{return sharedptr_partition->invariant; };
    
    void buildCLV(int tip_id, pll_msa_t *msa, GenotypeErrorModel *gtErrorModel, std::vector<double> &clv, bool normalize)const ;
    
    int initTipCLV(unsigned int tipClvIndex, double * clv)const;
    
    void showEigenDecomp(unsigned int float_precision) const;
    
    std::string getModelName() const {std::string str(model->name);
                                      return str;
                                      };
    
    ~Partition()
    {
        // if(sumtable)
        //         pll_aligned_free(sumtable);
        //   pll_partition_destroy(partition.get());
        
    }
    
};

#endif /* partition_hpp */





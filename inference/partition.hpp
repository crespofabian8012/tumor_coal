//
//  partition.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 27/06/2021.
//

#ifndef partition_hpp
#define partition_hpp

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


class Partition{
private:
    pll_partition_t * partition;
    Partition(){};
//    Partition(int numberTips,
//               int clvBuffers,
//               int numberStates,
//               int numberSites,
//               int numberRateMatrices,
//               int probMatrices,
//               int numberRateCats,
//               int numberScaleBuffers,
//               int statesPadded,
//               bool sse,
//               bool avx,
//               bool  avx2,
//               bool avx512,
//               bool  asc,
//               bool tipPatternCompression);
// 

public:
      
  
      unsigned int attributes;
      int numberTips;
      int clvBuffers;
      int numberStates;
      int numberSites;
      int numberRateMatrices;
      int probMatrices;
      int numberRateCats;
      int numberScaleBuffers;
      int statesPadded;
      bool sse;
      bool avx;
      bool  avx2;
      bool avx512;
      bool asc;
      bool tipPatternCompression;
      double * sumtable;
    
     Partition(int numberTips,
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
                   bool tipPatternCompression);
//    static Partition& getInstance(int numberTips,
//                                  int clvBuffers,
//                                  int numberStates,
//                                  int numberSites,
//                                  int numberRateMatrices,
//                                  int probMatrices,
//                                  int numberRateCats,
//                                  int numberScaleBuffers,
//                                  int statesPadded,
//                                  bool sse,
//                                  bool avx,
//                                  bool  avx2,
//                                  bool avx512,
//                                  bool  asc,
//                                  bool tipPatternCompression)
//           {
//               static Partition   instance( numberTips,
//                                      clvBuffers,
//                                      numberStates,
//                                      numberSites,
//                                      numberRateMatrices,
//                                      probMatrices,
//                                      numberRateCats,
//                                      numberScaleBuffers,
//                                      statesPadded,
//                                      sse,
//                                      avx,
//                                      avx2,
//                                      avx512,
//                                      asc,
//                                      tipPatternCompression); // Guaranteed to be destroyed.
//                                     // Instantiated on first use.
//               return instance;
//           }
    
    
 
   Partition(Partition const&)               = delete;
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
//void updateProbMatrices(
//     unsigned int * paramsIndexes,
//     unsigned int * matrixIndexes,
//     double * branchLengths,
//    unsigned int count);
    
void updateInvariantSites();

void updateInvariantSitesProportion(int paramsIndex, double propInvar);

void updatePartials(pll_operation_t * operations,
                    unsigned int count);

double   computeRootLogLikelihood(int clvIndex, int  scalerIndex, unsigned int * freqsIndex, double * persiteLnL);

double  computeEdgeLogLikelihood(int parentCLVIndex, int parentScalerIndex, int childCLVIndex, int childScalerIndex, int matrixIndex, unsigned int* freqsIndex,double * persiteLnL);
    
int  updateSumtable(unsigned int parentCLVIndex,unsigned int childCLVIndex,int parentScalerIndex, int childScalerIndex, const unsigned int *paramsIndexes);
    

int  computeLikelihoodDerivatives(int parentScalerIndex, int childScalerIndex, double branchLength, const unsigned int * paramsIndices, double * df, double * ddf);
 
pll_partition_t * getPartition();

~Partition()
{
   // if(sumtable)
   //         pll_aligned_free(sumtable);
        
    
    pll_partition_destroy(partition);
        
}
    
};

#endif /* partition_hpp */





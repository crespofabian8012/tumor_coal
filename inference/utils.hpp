/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/

/*
 * utils functions
 */
#ifndef utils_h
#define utils_h
//#include <vector>
//#include <algorithm>
//#include <numeric>
//#include <iostream>


#include "data_types.hpp"
#include "population.hpp"
#include "tree_node.hpp"


//class ProgramOptions;
//class TreeNode;
//class Population;
//using namespace std;

struct Utils
{
static void ReadParametersFromFastaFile(char *fileName, int &numCells, int &TotalNumSequences, int &numSites);
static void ReadFastaFile(char *fileName, vector<vector<int> > &ObservedData,  char **ObservedCellNames, ProgramOptions &programOptions);

static int CheckMatrixSymmetry(double matrix[4][4]);
static int WhichNucChar (char nucleotide);
static vector<int> SequenceToIntegers(char ** sequence, int length);
static vector<int> GenotypesToIntegers(char ** sequence, int length);
static int ChooseUniformState ( double *prob, long int *seed);
static char WhichIUPAC (int allele1, int allele2);
static char *WhichGenotypeFromIUPAC (int  iupac);
static int WhichGenotypeChar (char nucleotide);
static char WhichNuc (int nucleotide);
static char WhichConsensusBinary (int allele1, int allele2);
static char WhichMut (int state);
static void  normalizeVector(double *vector, int length);
static int compareIntDescending(const void *a, const void *b);
static int CompareGenotypes (int a1, int a2, int b1, int b2);

static double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym);

static void unscale(double * prob, unsigned int times);

static double DistanceBetweenSFS(int* SFS1, int*SFS, int numSNVs,  int numSites);


static int computeTajimaD(TreeNode **treeTips,  int numSites, int numCells);

static double ComputeESS(double *weights, int numberWeights);
//void computeUnfoldedISMSFS(int numSites,SiteStr* allSites,int numSNVs, int* SNVsites, int* SFS, int *numberDifferences);

static void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions &programOptions ) ;

static void computeGenotypesFreq(double freqs[10], pll_msa_t * msa);

static void InitNumberNodes(double *TotalBirthRate, double *TotalDeathRate, int *TotalN,  Population **populations, ProgramOptions &programOptions) ;


static double parameterMultiplierMCMCmove (double lengthInterval);
    
static   char * cb_serialize(const pll_rnode_t * node);

static  char * rTreeToNewick(const pll_rnode_t * root,
                             char * (*cb_serialize)(const pll_rnode_t *));
static char* appendCharToCharArray(char* array, char a);

   
static vector<long double>    potentialScaleReductionArray(int inner_size,int n, int m,  vector<vector<long double>> &means, vector<vector<long double>> &variances );
static    long double     potentialScaleReductionLongDouble(int inner_size,int n, int m,  vector<long double>& means, vector<long double> &variances );
    
static bool checkAllElementsGreaterThanZero(vector<long double> &values);
    
template<typename T, typename Iter_T>
static T computeAccum(Iter_T first, Iter_T last){
        
       T accum=T();
       for (auto it= first; it != last; ++it)
           accum += *it;
        
        return accum;
       // accum = std::accumulate(first, last, T());
    }
template<typename T>
static T mean2(const std::vector<T> &vec){
    
        size_t sz = vec.size();
    
        if (sz == 0)
            return T();
        if (sz == 1){
            T elem = vec[0];
            return elem;
         }
    
   // typename vector<T>::iterator first = &vec.begin();
    //typename vector<T>::iterator last = &vec.end();
    T mean = Utils::computeAccum(vec.begin(), vec.end());
    
    mean = mean / sz;

    return mean;
    }
    
template<typename InputIt, typename T, class BinaryOperation> constexpr
static   T accumulateOp(InputIt first, InputIt last,T init, BinaryOperation op)
    {
        
        for (auto it=first; it != last; ++it) {
            init = op(std::move(init), *it);
        }
        return init;
    }
template<typename T>
static T variance2(const std::vector<T> &vec)
    {
        size_t sz = vec.size();
        if (sz == 1)
            return 0.0;
        
        //  mean
        //typename vector<T>::iterator first = &vec.begin();
       // typename vector<T>::iterator last = &vec.end();
        T mean = Utils::computeAccum(vec.begin(), vec.end());
        
        mean = mean / sz;
        
        // variance
        auto variance_func = [&mean, &sz](T accumulator, const T& val)
        {
            return accumulator + ((val - mean)*(val - mean) / (sz - 1));
        };
        
       // return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
   
        return Utils::accumulateOp(vec.begin(), vec.end(), T(), variance_func);
        
    }

static long double mean(const std::vector<long double> &vec){
        
        size_t sz = vec.size();
        
        if (sz == 0)
            return 0.0;
        if (sz == 1){
            long double elem = vec[0];
            return elem;
        }
    
     long double accum=0.0;
     for (auto it= vec.begin(); it != vec.end(); ++it)
         accum += *it;

       long double mean = accum / sz;
        
        return mean;
    }
static long double variance(const std::vector<long double> &vec){
        size_t sz = vec.size();
        if (sz == 1)
            return 0.0;
    
       long double mean = Utils::mean(vec);
        
        // variance
        auto variance_func = [&mean, &sz](long double accumulator, const long double& val)
        {
            return accumulator + ((val - mean)*(val - mean) / (sz - 1));
        };
        
        long double result = accumulateOp(vec.begin(), vec.end(), 0.0, variance_func);
        return result;
        
    }
template <typename T>
static T* find(std::vector<T> &vec, const T& value)
    {
        auto it = std::find(vec.begin(), vec.end(), value);
        return (it != vec.end()) ? &(*it) : nullptr;
    }
template<typename T>
static T fastVariance(const std::vector<T> &vec)
    {
        int size = vec.size();
        
        double variance = 0;
        double t = vec[0];
        for (int i = 1; i < size; i++)
        {
            t += vec[i];
            double diff = ((i + 1) * vec[i]) - t;
            variance += (diff * diff) / ((i + 1.0) *i);
        }
        
        return variance / (size - 1);
    }
    
template <class OutputIter, class UnaryFunction>
static void apply(OutputIter first, OutputIter last, UnaryFunction f)
    {
        std::transform(first, last, first, f);
    }
template<typename T>
static vector<T>    potentialScaleReduction(size_t inner_size,size_t n, size_t m,  vector<vector<T>> &means, vector<vector<T>> &variances )
    {
        vector<T> elements;
        vector<T> varianceElements;
        vector<T> results;
        T meanOfMeans;
        T B;
        T W;
        T potentialScaleReduction ;
        for(size_t j=0; j <inner_size;++j)
        {
            
            std::transform( means.begin(), means.end(), elements.begin(), [&](vector<T> &vec){
                return vec.at(j);
            });
            std::transform( variances.begin(), variances.end(), varianceElements.begin(), [&](vector<T> &vec){
                return vec.at(j);
            });
            
            meanOfMeans = Utils::mean( elements);
            B = n * Utils::variance(elements) ;
            W = (1.0 /m)*Utils::mean(varianceElements) ;
            potentialScaleReduction= (1.0-1.0/n)* W+ (1.0/n)*B;
            results.push_back(potentialScaleReduction);
        }
        return results;
    }
template<typename T>
    static vector<T*>  vectorFromDoublePointer(T** doublePointer, size_t count)
    {
        vector<T*> result;
        T* current;
        for(size_t j=0; j <count;++j){
           
            current =doublePointer[j];
      
        
            result.push_back(doublePointer[j]);
        }
        return result;
    }
template <class DstType, class SrcType>
    bool IsType(const SrcType* src)
    {
        return dynamic_cast<const DstType*>(src) != nullptr;
    }
//template <class T>
//    void wrapArrayInVector( T *sourceArray, size_t arraySize, std::vector<T, std::allocator<T> > &targetVector ) {
//        typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
//        (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)((void *) &targetVector);
//        vectorPtr->_M_start = sourceArray;
//        vectorPtr->_M_finish = vectorPtr->_M_end_of_storage = vectorPtr->_M_start + arraySize;
//    }
//
//template <class T>
//    void releaseVectorWrapper( std::vector<T, std::allocator<T> > &targetVector ) {
//        typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
//        (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)((void *) &targetVector);
//        vectorPtr->_M_start = vectorPtr->_M_finish = vectorPtr->_M_end_of_storage = NULL;
//    }
    static vector<pll_rnode_t*>  filterHealthyTip(pll_rnode_t** doublePointer, size_t count, string &healthyTipLabel );
    static  void init_to_empty_str(char str[MAX_NAME]);
};

#endif /* utils_h */

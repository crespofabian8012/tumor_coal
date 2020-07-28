//
//  mcmc_parameter.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/7/20.
//

#ifndef mcmc_parameter_hpp
#define mcmc_parameter_hpp

#include <vector>
#include <stdio.h>

#include "data_utils.hpp"

template<typename T>
class MCMCParameter
{
public:
    T* currentValue;
    T* proposedValue;
    std::vector<T *> chainValues;
    MCMCParameter(T* currentValue, T* proposedValue);
    void safeCurrentValue();
    void rollabackCurrrentValue();
    ~MCMCParameter();
};
template<typename T>
class MCMCKernel
{
public:
    MCMCKernel();
    virtual T* proposedValue(T* currentValue, gsl_rng *randomGenerator)=0;
    virtual long double logKernel(T* currentValue, T* proposedValue)=0;
};
//template<typename T>
//class OptimalScaling:public MCMCKernel
//{
//public:
//    T* proposedValue(T* currentValue, gsl_rng *randomGenerator);
//    long double logKernel(T* currentValue, T* proposedValue);
//};

#endif /* mcmc_parameter_hpp */

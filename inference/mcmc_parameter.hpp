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
#include "random.h"

template <typename T>
class IMCMCKernel{
public:
    virtual T proposedValue(T& currentValue, gsl_rng *randomGenerator)=0;
    virtual  long double logKernel(T& currentValue, T& proposedValue)=0;
    virtual void setParameter(T& newParameter)=0;
    virtual T getParameter()=0;
    
};
template <typename T>
class MCMCKernel:public IMCMCKernel<T>
{
protected:
    T parameter;
public:
    MCMCKernel(T& argParameter){parameter=argParameter;};
    virtual T proposedValue(T& currentValue, gsl_rng *randomGenerator)=0;
    virtual  long double logKernel(T& currentValue, T& proposedValue)=0;
    T getParameter(){return parameter;}
    void setParameter(T& newParameter) {parameter = newParameter;}
    virtual ~MCMCKernel() = default;
};
template <typename T> class ProportionalScaling : public MCMCKernel<T> {
public:
    T parameter;
    ProportionalScaling(T& parameter):MCMCKernel<T>(parameter){};
    T proposedValue(T& currentValue, gsl_rng *randomGenerator) override{
        long double randomNumber= Random::randomUniformFromGsl2(randomGenerator);
        long double m = exp(2 * log(parameter) *(randomNumber -0.5));
        T result = m * currentValue ;
        return(result);
    };
    
   long double  logKernel(T& currentValue, T& proposedValue) override{
         long double result= log(proposedValue / currentValue);
         return(result);
    };
    
};
typedef ProportionalScaling<long double>    DoubleProportionalScaling;
typedef MCMCKernel<long double> DoubleMCMCKernel;
class DirichletKernel : public MCMCKernel<vector<long double>> {
public:
    vector<long double> parameter;
    DirichletKernel(vector<long double>& parameter):MCMCKernel<vector<long double>>(parameter){};
    vector<long double> proposedValue(vector<long double>& currentValue, gsl_rng *randomGenerator) override{
        vector< long double> newVector= currentValue;
        Random::randomDirichletFromVector (parameter, newVector);
        return(newVector);
    };
    
    long double  logKernel(vector<long double>& currentValue, vector<long double>& proposedValue) override{
        long double numerator= Distributions::DirichletDensity(currentValue, parameter, currentValue.size());
        long double  denominator = Distributions::DirichletDensity(proposedValue, parameter, currentValue.size());
        long double result = numerator / denominator;
        return(result);
    };
    
};
template<typename T>
class IMCMCParameter{
public:
    virtual void  saveCurrentValue()=0;
    virtual T  getCurrentValue()=0;
    virtual void  setValue(T& newValue )=0;
    virtual void  rollbackValue()=0;
    virtual void resetValues()=0;
    virtual void  saveCurrentValueToList()=0;
    
};
template<typename T>
class MCMCParameter:public IMCMCParameter<T>
{
public:
    T currentValue;
    T oldValue;
    std::vector<T> chainValues;
    T currentCumSum;
    T currentMean;
    T currenVariance;
    MCMCParameter(T& initialValue){
        currentValue = initialValue;
        oldValue = initialValue;
        chainValues.clear();
        currentCumSum = initialValue;
        currentMean = initialValue;
        currenVariance =T();//0.0 or 0 depending of the type T
    };
    void saveCurrentValue() override{
        oldValue=currentValue;
    };
    void saveCurrentValueToList() override{
        chainValues.push_back(currentValue);
        
        
        //currentCumSum += currentValue;
        
       // int path = chainValues.size() +1;
        //T term1 = (currentMean - currentValue) * (currentMean - currentValue);
       // T term2 = currenVariance + term1  / path;
        //currenVariance = path > 1 ? (path - 1.0) / path * (currenVariance + (currentMean - currentValue) * (currentMean - currentValue) / path) : 0.0;
        //currenVariance =  (path - 1.0) / path * (term2) ;
       // currentMean = ((path - 1.0) * currentMean + currentValue) / path;

    }
    void setValue( T& newValue) override{
        currentValue=newValue;
    };
    T getCurrentValue( ) override{
        return(currentValue);
    };
    void rollbackValue() override{
        currentValue=oldValue;
    };
    void resetValues() override{
        chainValues.clear();
    }
    T getCurrentMean(){
        return currentMean;
    }
    T getCurrentVariance(){
        return currenVariance;
    }
    T getCumSum(){
        return currentCumSum;
    }
    ~MCMCParameter();
};
template<typename T>
class IMCMCParameterWithKernel{
public:
    virtual void initKernel( T &kernelParameter)=0;
    virtual void initParameter(T &initialValue)=0;
    virtual void changeKernel(MCMCKernel<T>& newKernel)=0;
    virtual T getKernelParameter()=0;

    virtual void modifyKernelParameter(T &newParameter)=0;
    virtual T stepKernel(gsl_rng *randomGenerator)=0;
    virtual void resetParameterValues()=0;
    virtual void saveCurrentValue()=0;
    virtual void rollbackValue()=0;
    virtual void setParameterValue(T &newValue)=0;
};

typedef MCMCParameter<long double>  DoubleMCMCParameter;
typedef MCMCKernel<long double>  DoubleKernel;
class MCMCParameterWithKernel: public  IMCMCParameterWithKernel<long double>{
    //DoubleMCMCParameter *MCMCparameter;
    //DoubleMCMCKernel *kernel;
    IMCMCParameter<long double> *MCMCparameter;
    IMCMCKernel<long double > *kernel;
public:
    MCMCParameterWithKernel();
    MCMCParameterWithKernel(long  double initialValue, long double kernelParameter);
    void initKernel(long double &kernelParameter);
    void initParameter(long double &initialValue);
    void changeKernel(MCMCKernel<long double>& newKernel);
    long double getKernelParameter();
    IMCMCKernel<long double > *getKernel();
    IMCMCParameter<long double> *getParameter();
    void modifyKernelParameter(long double &newParameter);
    long double stepKernel(gsl_rng *randomGenerator);
    void resetParameterValues();
    void saveCurrentValue();
    void rollbackValue();
    void setParameterValue(long double & newValue);
};
typedef MCMCParameter<vector< long double>>  MCMCVector;
typedef MCMCKernel<vector< long double>> VectorKernel;
class MCMCVectorParameterWithKernel:public  IMCMCParameterWithKernel<vector<long double>>{
    IMCMCParameter<vector<long double>> *MCMCparameter;
    IMCMCKernel<vector<long double> > *kernel;
public:
    MCMCVectorParameterWithKernel();
    MCMCVectorParameterWithKernel(vector<long double> &nitialValue, vector<long double> &kernelParameter);
    
    void initKernel(vector<long double> &kernelParameter);
    void initParameter(vector<long double> &nitialValue);
    void changeKernel(MCMCKernel<vector<long double>>& newKernel);
    vector<long double> getKernelParameter();
    IMCMCKernel<vector< long double>> *getKernel();
    IMCMCParameter<vector<long double>> *getParameter();
    void modifyKernelParameter(vector<long double> &newParameter);
    vector<long double> stepKernel(gsl_rng *randomGenerator);
    void resetParameterValues();
    void saveCurrentValue();
    void rollbackValue();
    void setParameterValue(vector<long double> &newValue);
};
#endif /* mcmc_parameter_hpp */

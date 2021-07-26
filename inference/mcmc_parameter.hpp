//
//  mcmc_parameter.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/7/20.
//

#ifndef mcmc_parameter_hpp
#define mcmc_parameter_hpp


#include "data_utils.hpp"
#include "random.h"

#include <vector>
#include <stdio.h>


using std::placeholders::_1;

template <typename T>
class IMCMCKernel{
public:
    virtual T proposedValue(T& currentValue,const gsl_rng *randomGenerator)=0;
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
    virtual T proposedValue(T& currentValue,const gsl_rng *randomGenerator)=0;
    virtual  long double logKernel(T& currentValue, T& proposedValue)=0;
    T getParameter(){return parameter;}
    void setParameter(T& newParameter) {parameter = newParameter;}
    virtual ~MCMCKernel() = default;
};
template <typename T> class  SlidingWindow : public MCMCKernel<T> {
public:
    SlidingWindow(T& parameter):MCMCKernel<T>(parameter){};
    T proposedValue(T& currentValue, const gsl_rng *randomGenerator) override{
        long double randomNumber= Random::randomUniformFromGsl2(randomGenerator);
        T result = currentValue + (randomNumber-0.5) * 0.5 * this->getParameter() ;
        //this->getParameter() =windowSize is the distance  from the minimum to the  maximum value
        if (result <0 )
            result = -result;
        return result;
    };
   long double  logKernel(T& currentValue, T& proposedValue) override{
       long double result = 1.0;
         return(result);
    };
};
template <typename T> class  ProbabNormalProposal : public MCMCKernel<T> {
public:
   // T parameter;
    ProbabNormalProposal(T& parameter):MCMCKernel<T>(parameter){};
    T proposedValue(T& currentValue, const gsl_rng *randomGenerator) override{
        long double randomNumber= Random::randomNormalGreaterThan(0, this->getParameter(), 0, true, randomGenerator, NULL);
        T proposal = currentValue + randomNumber;
        if(proposal < 0){
            proposal = -proposal;
        }
        if(proposal > 1){
             proposal = abs(2-proposal);
        }
        return proposal;
    };
    
   long double  logKernel(T& currentValue, T& proposedValue) override{
         long double result= 1.0;
         return(result);
    };
};
template <typename T> class ProportionalScaling : public MCMCKernel<T> {
public:
    ProportionalScaling(T& parameter):MCMCKernel<T>(parameter){};
    T proposedValue(T& currentValue, const gsl_rng *randomGenerator) override{
        long double randomNumber= Random::randomUniformFromGsl2(randomGenerator);
        long double m = exp(2 * log(this->getParameter()) *(randomNumber -0.5));
        T result = m * currentValue ;
        return(result);
    };
   long double  logKernel(T& currentValue, T& proposedValue) override{
         long double result= log(proposedValue / currentValue);
         return(result);
    };
};
template <typename T> class ProportionalScalingWithBounds : public MCMCKernel<T> {
public:
    T from;
    T to;
    ProportionalScalingWithBounds(T& parameter, T& fromPar, T& toPar):MCMCKernel<T>(parameter){
        from = fromPar;
        to = toPar;
    };
    T proposedValue(T& currentValue,const gsl_rng *randomGenerator) override{
        long double randomNumber= Random::randomUniformFromGsl2(randomGenerator);
        long double m = exp(2 * log(this->getParameter()) *(randomNumber -0.5));
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
class DirichletKernel : public MCMCKernel<std::vector<long double>> {
public:
    std::vector<long double> parameter;
    DirichletKernel(std::vector<long double>& parameter):MCMCKernel<std::vector<long double>>(parameter){};
    std::vector<long double> proposedValue(std::vector<long double>& currentValue, const gsl_rng *randomGenerator) override{
        std::vector< long double> newVector= currentValue;
        Random::randomDirichletFromVector (parameter, newVector,true, randomGenerator, NULL);
        return(newVector);
    };
    
    long double  logKernel(std::vector<long double>& currentValue, std::vector<long double>& proposedValue) override{
        std::vector<long double> dummy;
        long double numerator= Distributions::LogDirichletDensity(currentValue, parameter, dummy );
        long double  denominator = Distributions::LogDirichletDensity(proposedValue, parameter, dummy);
        long double result = numerator - denominator;
        return(result);
    };
    
};
template<class T>
using LogPriorDensity = long double(*)(T &, T &, T &);
template<typename T>
class IMCMCParameter{
public:
   
    virtual void  saveCurrentValue()=0;
    virtual T  getCurrentValue()=0;
    virtual void  setValue(T& newValue )=0;
    virtual void  rollbackValue()=0;
    virtual void resetValues()=0;
    virtual void  saveCurrentValueToList()=0;
    virtual long double  getLogPriorCurrentValue()=0;
    virtual void setLogPriorCurrentValue(long double prior_value )=0;
    virtual long double logPriorDensityCurrentvalue(T& second, T& third)=0;
    virtual void setPrior( LogPriorDensity<T> lPriorDensity) =0;
    virtual long double  getOldLogPriorValue()=0;
    virtual void setSecondParLogPrior (T& secondPar ) =0;
    virtual void setThirdParLogPrior (T& thirdPar ) =0;
    virtual void setCurrentMean (T& mean) =0;
    virtual void setCurrentVariance (T& variance) =0;
    virtual std::vector<T>  getChainValues()=0;
    virtual T getCurrentMean()=0;
    virtual T getCurrentVariance()=0;
    virtual std::string getName()=0;
};
template <typename T>
 T operator+(const T& a, const T& b)
{
  return a + b;
}
template <typename T>
bool operator <( const T &lhs, const long double  &rhs )
{
   return ( lhs < rhs );
}
template<typename T>
class MCMCParameter:public IMCMCParameter<T>
{
public:
    int indexEigenMatrix;
    std::string name;
    T currentValue;
    long double logPriorCurrentValue;
    T oldValue;
    long double logPriorOldValue;
    std::vector<T> chainValues;
    T currentCumSum;
    T currentMean;
    T currentVariance;
    bool hasPrior;
    LogPriorDensity<T> logPriorDensity;
    T secondParPrior;
    T thirdParPrior;
    bool affectsLikelihood;
    MCMCParameter(std::string namePar, T& initialValue, bool hasPriorPar,LogPriorDensity<T> lPriorDensity ,  T& second, T& third){
        name =namePar;
        currentValue = initialValue;
        oldValue = initialValue;
        chainValues.clear();
        currentCumSum = initialValue;
        currentMean = initialValue;
        currentVariance =T();//0.0 or 0 depending of the type T
        logPriorOldValue = 0.0;
        logPriorDensity = lPriorDensity;
        hasPrior= hasPriorPar;
        secondParPrior = second;
        thirdParPrior = third;
        if (hasPrior){
            logPriorCurrentValue = logPriorDensity(initialValue, second, third);
            affectsLikelihood = true;
        }
        else{
            logPriorCurrentValue =0.0;
            affectsLikelihood = false;
        }
    };
    MCMCParameter(std::string namePar, T& initialValue,  long double initialLogPriorDensity){
        name =namePar;
        currentValue = initialValue;
        oldValue = initialValue;
        chainValues.clear();
        currentCumSum = T();// initialValue;
        currentMean = T();//initialValue;
        currentVariance =T();//0.0 or 0 depending of the type T
        secondParPrior = T();
        thirdParPrior = T();
        logPriorOldValue = 0.0;
        hasPrior= false;
        logPriorDensity = nullptr;
        logPriorCurrentValue =initialLogPriorDensity;
        affectsLikelihood = false;
        //logPriorCurrentValue = logPriorCurrentValue;

    };
    void setPrior( LogPriorDensity<T> lPriorDensity) override{
          logPriorDensity = lPriorDensity;
          hasPrior = true;
          affectsLikelihood = true;
      }
    void saveCurrentValue() override{
        oldValue=currentValue;
        logPriorOldValue = logPriorCurrentValue;
    };
    void saveCurrentValueToList() override{
        chainValues.push_back(currentValue);
    }
    void setPriorSecondParameter( T& par){
        secondParPrior = par;
        
    }
    void setPriorThirdParameter( T& par){
        thirdParPrior = par;
        
    }
    void setValue( T& newValue) override{
        currentValue=newValue;
        if (hasPrior)
                   logPriorCurrentValue = logPriorDensity(currentValue, secondParPrior, thirdParPrior);
               else
                   logPriorCurrentValue =0.0;
    };
    T getCurrentValue( ) override{
        return(currentValue);
    };
    void rollbackValue() override{
        currentValue=oldValue;
        logPriorCurrentValue = logPriorOldValue ;
    };
    void resetValues() override{
        chainValues.clear();
        currentMean=T();
        currentVariance=T();
    }
    T getCurrentMean() override{
        return currentMean;
    }
    T getCurrentVariance()override{
        return currentVariance;
    }
    T getCumSum(){
        return currentCumSum;
    }
    
    long double logPriorDensityCurrentvalue(T& second, T& third)override{
        
        if (hasPrior)
                logPriorCurrentValue = logPriorDensity(currentValue, second, third);
        else
            logPriorCurrentValue=0.0;
        
        return logPriorCurrentValue;
    }
    long double  getLogPriorCurrentValue() override{
           return logPriorCurrentValue;
    }
    long double  getOldLogPriorValue() override{
              return logPriorOldValue;
       }
    void setLogPriorCurrentValue (long double prior_value ) override{
          logPriorCurrentValue = prior_value;
    }
    void setSecondParLogPrior (T& secondPar ) override{
             secondParPrior = secondPar;
       }
    void setThirdParLogPrior (T& thirdPar ) override{
          thirdParPrior = thirdPar;
    }
    std::vector<T>  getChainValues()override{
        
        return chainValues;
    }
    void setCurrentMean (T& mean)override{
        
        currentMean= mean;
    }
    void setCurrentVariance (T& variance)override{
        
        currentVariance = variance;
    }
    std::string getName()override{return name;}
    ~MCMCParameter();
};
template<typename T>
class IMCMCParameterWithKernel{
public:
    virtual void initKernel( T &kernelParameter)=0;
    virtual void initParameter(std::string name,T &initialValue,  long double logPriorDensity)=0;
    virtual void changeKernel(MCMCKernel<T>& newKernel)=0;
    virtual T getKernelParameter()=0;

    virtual void modifyKernelParameter(T &newParameter)=0;
    virtual T stepKernel(const gsl_rng *randomGenerator)=0;
    virtual void resetParameterValues()=0;
    virtual void saveCurrentValue()=0;
    virtual void  saveCurrentValueToList()=0;
    virtual void rollbackValue()=0;
    virtual void setParameterValue(T &newValue)=0;
    virtual std::vector<T> autoCorrelation()=0;
    virtual std::vector<T> autoCovariance()=0;
    virtual long double ESS()=0;
    
};

typedef MCMCParameter<long double>  DoubleMCMCParameter;
typedef MCMCKernel<long double>  DoubleKernel;
class MCMCParameterWithKernel: public  IMCMCParameterWithKernel<long double>{
    //DoubleMCMCParameter *MCMCparameter;
    //DoubleMCMCKernel *kernel;
    std::vector<linalgebra::complex> chainComplexValues;
    IMCMCParameter<long double> *MCMCparameter;
    IMCMCKernel<long double > *kernel;
public:
    MCMCParameterWithKernel();
    MCMCParameterWithKernel(std::string name,long  double initialValue, long double kernelParameter,
                            long double logPriorDensity);
    MCMCParameterWithKernel(std::string name,long  double & initialValue, long double kernelParameter,
                            LogPriorDensity<long double> lPriorDensity,  long double& second, long double& third);
    void initKernel(long double &kernelParameter);
    void initParameter(std::string name,long double &initialValue,  long double logPriorDensity);
    void changeKernel(MCMCKernel<long double>& newKernel);
    long double getKernelParameter();
    IMCMCKernel<long double > *getKernel();
    IMCMCParameter<long double> *getParameter();
    void modifyKernelParameter(long double &newParameter);
    long double stepKernel(const gsl_rng *randomGenerator);
    void resetParameterValues();
    void saveCurrentValue();
    void saveCurrentValueToList();
    void rollbackValue();
    void setParameterValue(long double & newValue);
    std::vector<long double> autoCorrelation();
    std::vector<long double> autoCovariance();
    long double ESS();
};
typedef MCMCParameter<std::vector< long double>>  MCMCVector;
typedef MCMCKernel<std::vector< long double>> VectorKernel;
class MCMCVectorParameterWithKernel:public  IMCMCParameterWithKernel<std::vector<long double>>{
    IMCMCParameter<std::vector<long double>> *MCMCparameter;
    IMCMCKernel<std::vector<long double> > *kernel;
public:
    MCMCVectorParameterWithKernel();
    MCMCVectorParameterWithKernel( std::string name, std::vector<long double> &nitialValue,
                                   std::vector<long double> &kernelParameter,  long double logPriorDensity);
    
    void initKernel(std::vector<long double> &kernelParameter);
    void initParameter(std::string name,std::vector<long double> &nitialValue,  long double logPriorDensity);
    void changeKernel(MCMCKernel<std::vector<long double>>& newKernel);
    std::vector<long double> getKernelParameter();
    IMCMCKernel<std::vector< long double>> *getKernel();
    IMCMCParameter<std::vector<long double>> *getParameter();
    void modifyKernelParameter(std::vector<long double> &newParameter);
    std::vector<long double> stepKernel(const gsl_rng *randomGenerator);
    void resetParameterValues();
    void saveCurrentValue();
    void saveCurrentValueToList();
    void rollbackValue();
    void setParameterValue(std::vector<long double> &newValue);
    std::vector<std::vector<long double>> autoCorrelation();
    std::vector<std::vector<long double>> autoCovariance();
    long double ESS();
  
};
#endif /* mcmc_parameter_hpp */

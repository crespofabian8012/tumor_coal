//
//  mcmc_parameter.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/7/20.
//
//#include  <algorithm>

#include "autocovariance.hpp"
#include "autocorrelation.hpp"
#include "mcmc_parameter.hpp"
#include "random.h"
#include "utils.hpp"


MCMCParameterWithKernel::MCMCParameterWithKernel(std::string name, long double initialValue,long double kernelParameter,  long double logPriorDensity){
    MCMCparameter = new  DoubleMCMCParameter(name, initialValue,  logPriorDensity);
    kernel = new DoubleProportionalScaling(kernelParameter);
    kernel->setParameter(kernelParameter);
};
MCMCParameterWithKernel::MCMCParameterWithKernel(std::string name, long  double & initialValue, long double kernelParameter,  LogPriorDensity<long double> lPriorDensity,  long double& second, long double& third){
    
    long double initialLogPrior = lPriorDensity(initialValue, second,  third);
    MCMCparameter = new  DoubleMCMCParameter(name, initialValue, initialLogPrior);
    MCMCparameter->setPrior(  lPriorDensity );
    kernel = new DoubleProportionalScaling(kernelParameter);
};
void MCMCParameterWithKernel::changeKernel(DoubleMCMCKernel& newKernel){
    kernel = &newKernel;
};
long double MCMCParameterWithKernel::getKernelParameter()
{
    return kernel->getParameter();
};
void MCMCParameterWithKernel::modifyKernelParameter(long double &newParameter)
{
    return kernel->setParameter(newParameter);
};
long double MCMCParameterWithKernel::stepKernel(const gsl_rng *randomGenerator){
    long double currentValue= MCMCparameter->getCurrentValue();
    long double newValue= kernel->proposedValue(currentValue, randomGenerator);
    MCMCparameter->setValue(newValue );
    return(newValue);
};

void MCMCParameterWithKernel::resetParameterValues(){
    MCMCparameter->resetValues();
};

void MCMCParameterWithKernel::initKernel(long double& kernelParameter){
    DoubleProportionalScaling* mcmckernel = new DoubleProportionalScaling(kernelParameter);
    kernel = mcmckernel;
};
void MCMCParameterWithKernel::initParameter(std::string name, long double& initialValue, long double logPriorDensity){
    MCMCparameter = new  DoubleMCMCParameter(name, initialValue, logPriorDensity);
};
IMCMCKernel<long double> * MCMCParameterWithKernel::getKernel(){
    return kernel;
};
IMCMCParameter<long double> * MCMCParameterWithKernel::getParameter(){
    return MCMCparameter;
};
void MCMCParameterWithKernel::saveCurrentValue(){
    MCMCparameter->saveCurrentValue();
};
void MCMCParameterWithKernel::saveCurrentValueToList(){
    
    MCMCparameter->saveCurrentValueToList();
    
    int path = MCMCparameter->getChainValues().size() ;
    long double currentMean = MCMCparameter->getCurrentMean();;
    long double currentValue = MCMCparameter->getCurrentValue();
    long double currentVariance;
    
    if (path > 1){
        currentVariance = MCMCparameter->getCurrentVariance() ;
        
        long double diff = currentMean -currentValue;
        long double term1 = diff * diff;
        long double term2 = currentVariance + term1  / path;
        long double factor = (path - 1.0) / path ;
        currentVariance = factor * (currentVariance + diff * diff / path);
        currentVariance =  (path - 1.0) / path * (term2) ;
        MCMCparameter->setCurrentVariance(currentVariance);
        
    }
    currentMean = ((path - 1.0) * currentMean + currentValue) / path;
    MCMCparameter->setCurrentMean(currentMean);
};

void MCMCParameterWithKernel::rollbackValue(){
    MCMCparameter->rollbackValue();
};
void MCMCParameterWithKernel::setParameterValue(long double& newValue){
    MCMCparameter->setValue(newValue);
}
MCMCVectorParameterWithKernel::MCMCVectorParameterWithKernel(std::string name,std::vector<long double> &initialValue, std::vector<long double> &kernelParameter,  long double logPriorDensity){
    MCMCparameter = new  MCMCVector(name, initialValue,   logPriorDensity);
    kernel = new DirichletKernel(kernelParameter);
};

void MCMCVectorParameterWithKernel::initKernel(std::vector<long double> &kernelParameter){
    kernel = new DirichletKernel(kernelParameter);
};
void MCMCVectorParameterWithKernel::initParameter(std::string name, std::vector<long double> &initialValue,  long double logPriorDensity){
    MCMCparameter = new  MCMCVector(name, initialValue, logPriorDensity);
};
void MCMCVectorParameterWithKernel::changeKernel(MCMCKernel<std::vector<long double>>& newKernel){
    
    kernel = &newKernel;
};
std::vector<long double> MCMCVectorParameterWithKernel::getKernelParameter(){
    std::vector<long double> result =kernel->getParameter();
    return result;
};
IMCMCKernel<std::vector< long double>> *MCMCVectorParameterWithKernel::getKernel(){
    return kernel;
};
IMCMCParameter<std::vector<long double>>* MCMCVectorParameterWithKernel::getParameter(){
    return MCMCparameter;
};
void MCMCVectorParameterWithKernel::modifyKernelParameter(std::vector<long double> &newParameter){};
std::vector<long double> MCMCVectorParameterWithKernel::stepKernel(const gsl_rng *randomGenerator){
    std::vector< long double> currentVector = MCMCparameter->getCurrentValue();
    std::vector< long double> newVector = kernel->proposedValue(currentVector, randomGenerator);
    MCMCparameter->setValue(newVector);
    return newVector;
};
void MCMCVectorParameterWithKernel::resetParameterValues(){
    MCMCparameter->resetValues();
};
void MCMCVectorParameterWithKernel::saveCurrentValue(){
    MCMCparameter->saveCurrentValue();
};
void MCMCVectorParameterWithKernel::saveCurrentValueToList(){
    // MCMCparameter->saveCurrentValueToList();
    
};
void MCMCVectorParameterWithKernel::rollbackValue(){
    MCMCparameter->rollbackValue();
};
void MCMCVectorParameterWithKernel::setParameterValue(std::vector<long double> &newValue){
    MCMCparameter->setValue(newValue);
};
std::vector<long double> MCMCParameterWithKernel::autoCorrelation(){
    
    int n= MCMCparameter->getChainValues().size();
    
    int N = Utils::next2Power( n);
    
    std::vector<long double> padded = MCMCparameter->getChainValues();
    
    padded.resize(padded.size() + (N - padded.size() % N) % N, 0.0);
    
    long double mean = MCMCparameter->getCurrentMean();
    for (int i=0;i<n;i++) {
        padded[i]=padded[i]-mean;
    }
    
    std::vector<long double> ans(2* N);
    
    Utils::correl(padded,padded,ans);
    //long double variance =  MCMCparameter->getCurrentVariance();
    for (int i = 0; i < n; ++i) {
        ans[i]= ans[i]/(N*N*2);
        ans[i] = ans[i] / ans[0];
    }
    
    return ans ;
}
std::vector<long double> MCMCParameterWithKernel::autoCovariance(){
    
    int n= MCMCparameter->getChainValues().size();
    std::vector<long double> values = MCMCparameter->getChainValues();
    std::vector<long double> autoCorrelations=autoCorrelation();
    
    std::vector<long double> autoCov(n);
    
    //    using boost::accumulators::accumulator_set;
    //    using boost::accumulators::stats;
    //    using boost::accumulators::tag::variance;
    
    long double variance =  MCMCparameter->getCurrentVariance();
    
    for (int i = 0; i < n; ++i) {
        autoCov[i]= autoCorrelations[i] *  variance;
    }
    return autoCov ;
}
long double MCMCParameterWithKernel::ESS(){
    
    long double resultCorr =0.0;
    // long double resultCov =0.0;
    int i;
    std::vector<long double > values=MCMCparameter->getChainValues();
    if ( std::adjacent_find( values.begin(), values.end(), std::not_equal_to<>() ) == values.end() )
    {
        return 1.0;
    }
    //std::vector<long double> autoCov=autoCovariance();
    std::vector<long double> autoCov;
    std::vector<long double> autoCorr;
    
    //mcmc_utils::autocovariance<long double>(MCMCparameter->getChainValues(), autoCov);
    mcmc_utils::autocorrelationVector<long double>(MCMCparameter->getChainValues(), autoCorr);
    //autoCov[0] is equal to MCMCparameter->getCurrentVariance()
    
    // long double varianceChain = autoCov[0]* MCMCparameter->getChainValues().size() / (MCMCparameter->getChainValues().size() -1);
    
    long double consecutiveCorr = 0.0;
    long double sumConsecutiveCorr = 0.0;
    for (i = 0; i < autoCorr.size()-1; i+=2){

        consecutiveCorr = autoCorr[i] +autoCorr[i+1];
        if (consecutiveCorr <= 0)
            break;
        sumConsecutiveCorr+=consecutiveCorr;
    }
    resultCorr = MCMCparameter->getChainValues().size() * (1.0 / (-autoCorr[0] + 2.0 * sumConsecutiveCorr) );
    assert(resultCorr>=0);
    return resultCorr;
    
}
std::vector<std::vector<long double>> MCMCVectorParameterWithKernel::autoCorrelation(){
    std::vector<std::vector<long double>> result;
    return result;
}
std::vector<std::vector<long double>> MCMCVectorParameterWithKernel::autoCovariance(){
    std::vector<std::vector<long double>> result;
    return result;
    
}
long double MCMCVectorParameterWithKernel::ESS(){
    
    long double result=0.0;
    return result;
}

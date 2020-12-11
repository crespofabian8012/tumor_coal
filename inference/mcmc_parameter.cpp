//
//  mcmc_parameter.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/7/20.
//

#include "mcmc_parameter.hpp"
#include "random.h"
MCMCParameterWithKernel::MCMCParameterWithKernel(long double initialValue,long double kernelParameter){
    MCMCparameter = new  DoubleMCMCParameter(initialValue);
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
long double MCMCParameterWithKernel::stepKernel(gsl_rng *randomGenerator){
    long double currentValue= MCMCparameter->getCurrentValue();
    long double newValue= kernel->proposedValue(currentValue, randomGenerator);
    MCMCparameter->setValue(newValue);
    return(newValue);
};

void MCMCParameterWithKernel::resetParameterValues(){
    MCMCparameter->resetValues();
};

void MCMCParameterWithKernel::initKernel(long double& kernelParameter){
    DoubleProportionalScaling* mcmckernel = new DoubleProportionalScaling(kernelParameter);
    kernel = mcmckernel;
};
void MCMCParameterWithKernel::initParameter(long double& initialValue){
    MCMCparameter = new  DoubleMCMCParameter(initialValue);
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

void MCMCParameterWithKernel::rollbackValue(){
    MCMCparameter->rollbackValue();
};
void MCMCParameterWithKernel::setParameterValue(long double & newValue){
    MCMCparameter->setValue(newValue);
}
MCMCVectorParameterWithKernel::MCMCVectorParameterWithKernel(vector<long double> &initialValue, vector<long double> &kernelParameter){
    MCMCparameter = new  MCMCVector(initialValue);
    kernel = new DirichletKernel(kernelParameter);
};

void MCMCVectorParameterWithKernel::initKernel(vector<long double> &kernelParameter){
    kernel = new DirichletKernel(kernelParameter);
};
void MCMCVectorParameterWithKernel::initParameter(vector<long double> &initialValue){
     MCMCparameter = new  MCMCVector(initialValue);
};
void MCMCVectorParameterWithKernel::changeKernel(MCMCKernel<vector<long double>>& newKernel){
    
    kernel = &newKernel;
};
vector<long double> MCMCVectorParameterWithKernel::getKernelParameter(){
    vector<long double> result =kernel->getParameter();
    return result;
};
IMCMCKernel<vector< long double>> *MCMCVectorParameterWithKernel::getKernel(){
    return kernel;
};
IMCMCParameter<vector<long double>>* MCMCVectorParameterWithKernel::getParameter(){
    return MCMCparameter;
};
void MCMCVectorParameterWithKernel::modifyKernelParameter(vector<long double> &newParameter){};
vector<long double> MCMCVectorParameterWithKernel::stepKernel(gsl_rng *randomGenerator){
    vector< long double> currentVector = MCMCparameter->getCurrentValue();
    vector< long double> newVector = kernel->proposedValue(currentVector, randomGenerator);
    MCMCparameter->setValue(newVector);
    return newVector;
};
void MCMCVectorParameterWithKernel::resetParameterValues(){
    MCMCparameter->resetValues();
};
void MCMCVectorParameterWithKernel::saveCurrentValue(){
        MCMCparameter->saveCurrentValue();
};
void MCMCVectorParameterWithKernel::rollbackValue(){
    MCMCparameter->rollbackValue();
};
void MCMCVectorParameterWithKernel::setParameterValue(vector<long double> &newValue){
    MCMCparameter->setValue(newValue);
};

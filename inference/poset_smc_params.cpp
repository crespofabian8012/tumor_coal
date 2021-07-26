//
//  poset_smc_params.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 17/03/2021.
//

#include "poset_smc_params.hpp"

PosetSMCParams::PosetSMCParams(int numberClones,
                  int sampleSize,
                  pll_msa_t *msa,
                  std::vector<int> &positions,
                 ProgramOptions &programOptions):
numberClones(numberClones),sampleSize(sampleSize)
{
    this->msa= msa;
    this->positions= positions;
    this->programOptions = &programOptions;
    
    assert(positions.size()==sampleSize);
        if (positions.size()==0)
              std::cout << "\n The particle cannot be initialized without tips. \n"<<std::endl;
        
}
ProgramOptions& PosetSMCParams::getProgramOptions(){
    
    return *programOptions;
    
}
void PosetSMCParams::set(std::shared_ptr<MCMCParameterWithKernel> thetaPar,
                std::shared_ptr<MCMCParameterWithKernel> seqErrorPar,
                std::shared_ptr<MCMCParameterWithKernel> dropoutErrorPar,
                std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationDeltaTsPar,
                std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationToriginSTDsPar,
                std::shared_ptr<MCMCVectorParameterWithKernel>proportionsPar
                )
{
    

    this->proportions = proportionsPar;
    this->theta = thetaPar;
    populationDeltaTs=   populationDeltaTsPar;
    populationToriginSTDs=   populationToriginSTDsPar;
    

}
std::shared_ptr<MCMCParameterWithKernel> PosetSMCParams::getTheta(){
    
    return theta;
}
std::shared_ptr<MCMCParameterWithKernel> PosetSMCParams::getPopulationDeltaT(int i){
    
    return populationDeltaTs[i];
}
std::shared_ptr<MCMCParameterWithKernel> PosetSMCParams::getPopulationToriginSTD(int i){
    
    return populationToriginSTDs[i];
}
std::shared_ptr<MCMCVectorParameterWithKernel> PosetSMCParams::getProportionVector(){
    
    return proportions;
}
int  PosetSMCParams::getSampleSize() const{
    
    return sampleSize;
}
int  PosetSMCParams::getNumClones() const{
    
    return numberClones;
}

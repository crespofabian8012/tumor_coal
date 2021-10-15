//
//  poset_smc_params.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 17/03/2021.
//

#include "poset_smc_params.hpp"


PosetSMCParams::PosetSMCParams(int numberClones,
                               int sampleSize,
                               std::vector<int> &sampleSizes,
                               unsigned int num_sites,
                               MSA *msa,
                               //pll_msa_t *msa,
                               const Partition *partition,
                               PLLBufferManager *const pll_buffer_manager,
                               std::vector<int> &positions,
                               ProgramOptions &programOptions,GenotypeErrorModel *gtErrorModel,
                               double theta,
                               std::vector<double> populationDeltaTs,
                               std::vector<double> populationToriginSTDs,
                               std::vector<double>proportions,
                               std::vector<std::vector<double>> coalTimesModelTimePerPopulation):
numberClones(numberClones),sampleSize(sampleSize), msa(msa), partition(partition),
pll_buffer_manager(pll_buffer_manager), gtErrorModel(gtErrorModel), sampleSizes(sampleSizes)
{
    
    this->positions= positions;
    this->programOptions = &programOptions;
    this->coalTimesModelTimePerPopulation = coalTimesModelTimePerPopulation;
    this->populationDeltaTs = populationDeltaTs;
    this->populationToriginSTDs = populationToriginSTDs;
    this->proportions= proportions;
    this->theta = theta;
    this->verbose = 1;
}

//PosetSMCParams::PosetSMCParams(const PosetSMCParams &original ):
//numberClones(original.numberClones),sampleSize(original.sampleSize), msa(original.msa), partition(original.partition),
//pll_buffer_manager(original.pll_buffer_manager), gtErrorModel(original.gtErrorModel), sampleSizes(original.sampleSizes)
//{
//    
//    this->positions= original.positions;
//    this->programOptions = original.programOptions;
//    this->coalTimesModelTimePerPopulation = original.coalTimesModelTimePerPopulation;
//    this->populationDeltaTs = original.populationDeltaTs;
//    this->populationToriginSTDs = original.populationToriginSTDs;
//    this->proportions= original.proportions;
//    this->theta = original.theta;
//    this->verbose = original.verbose;
//}
ProgramOptions& PosetSMCParams::getProgramOptions(){
    
    return *programOptions;
    
}
//void PosetSMCParams::set(std::shared_ptr<MCMCParameterWithKernel> thetaPar,
//                         std::shared_ptr<MCMCParameterWithKernel> seqErrorPar,
//                         std::shared_ptr<MCMCParameterWithKernel> dropoutErrorPar,
//                         std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationDeltaTsPar,
//                         std::vector<std::shared_ptr<MCMCParameterWithKernel>> populationToriginSTDsPar,
//                         std::shared_ptr<MCMCVectorParameterWithKernel>proportionsPar
//                         )
//{
//    this->proportions = proportionsPar;
//    this->theta = thetaPar;
//    populationDeltaTs=   populationDeltaTsPar;
//    populationToriginSTDs=   populationToriginSTDsPar;
//
//
//}
//std::shared_ptr<MCMCParameterWithKernel> PosetSMCParams::getTheta(){
//
//    return theta;
//}
//std::shared_ptr<MCMCParameterWithKernel> PosetSMCParams::getPopulationDeltaT(int i){
//
//    return populationDeltaTs[i];
//}
//std::shared_ptr<MCMCParameterWithKernel> PosetSMCParams::getPopulationToriginSTD(int i){
//
//    return populationToriginSTDs[i];
//}
//std::shared_ptr<double> PosetSMCParams::getProportionVector(){
//
//    return proportions;
//}
int  PosetSMCParams::getSampleSize() const{
    
    return sampleSize;
}
int  PosetSMCParams::getNumClones() const{
    
    return numberClones;
}
double PosetSMCParams::getPopulationEvent(int idx_population, int idx_event){
    
    assert(idx_population< numberClones);
    assert(idx_event < coalTimesModelTimePerPopulation.at(idx_population).size());
    return coalTimesModelTimePerPopulation.at(idx_population)[idx_event];
}
void PosetSMCParams::buildListEventTimesPerPopulation(){
    
 //   long double waitingTime, timeNextEvent;
    
//    std::vector<std::vector<double>> coalTimesModelTimePerPopulation(numberClones );
//    bool isCoalescentEvent=false;
//    double currentTime=0.0;
//    double currentTimeKingman=0.0;
//    Population *pop;
//    for(size_t i=0; i <  numberClones; i++){
//
//        pop= result->getPopulationByIndex( i);
//        if (pop->numActiveGametes >1){
//
//            currentTime = result->getHeightModelTime() ;
//            currentTimeKingman = Population::FmodelTstandard (currentTime , pop->timeOriginSTD, pop->delta,   params.getProgramOptions().K);//this is current time in Kingman coalescent
//
//            waitingTime= pop->proposeTimeNextCoalEvent(random, pop->numActiveGametes, params.getProgramOptions().K);
//
//            currentTimeKingman = currentTimeKingman + waitingTime;
//
//            timeNextEvent =   Population::GstandardTmodel(currentTimeKingman, pop->timeOriginSTD, pop->delta, params.getProgramOptions().K);;
//            isCoalescentEvent=true;
//
//        }
//        else{//last event
//
//            if (result->getNumberPopulations()==1 || pop == oldestPop){
//                timeNextEvent = pop->timeOriginSTD;
//            }
//            else{
//
//                timeNextEvent= pop->timeOriginSTD * pop->x / oldestPop->x;
//            }
//            isCoalescentEvent=false;
//        }
//
//
//        long double logLikNextCoal = 0.0;
//
//        logLikNextCoal = chosenPop->logConditionalLikelihoodNextCoalescentTime(minTimeNextEvent, result->getHeightModelTime(), chosenPop->numActiveGametes, params.getProgramOptions().K);
//        log_w = log_w +logLikNextCoal;
//        long double logLikNoCoal;
//        long double logDensityTorigin;
//        //do we  need to add loglik of no coal events in  other populations?
//
//        for(size_t i=0; i <  result->getNumberPopulations(); i++){
//            pop= result->getPopulationByIndex( i);
//            // do we need to include this in the node->ln_likelihood?
//            //logDensityTorigin = pop->LogDensityTime(pop->timeOriginSTD);
//
//            //node->ln_likelihood = node->ln_likelihood +logDensityTorigin;
//
//            if (pop!=chosenPop){
//                logLikNoCoal=pop->LogProbNoCoalescentEventBetweenTimes(result->getHeightModelTime(), minTimeNextEvent,  pop->numActiveGametes, pop->timeOriginSTD, pop->delta, params.getProgramOptions().K);
//
//                log_w = log_w +logLikNoCoal;
//            }
//        }
//
//
    }

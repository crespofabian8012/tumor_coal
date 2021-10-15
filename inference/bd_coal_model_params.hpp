//
//  bd_coal_model_params.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#ifndef bd_coal_model_params_hpp
#define bd_coal_model_params_hpp

#include <vector>

using namespace  std;
class BDCoalPriorParams{
public:
    double lambdaExpPriorScaledGrowthRate;
    double lambdaExpPriorTheta;
    
    double meanBetaPriorSeqError;
    double varBetaPriorSeqError;
    
    double meanBetaPriorADOError;
    double varBetaPriorADOError;
    
    double paramMultiplierProposalScaledGrowthRate;
    double paramMultiplierProposalTheta;
    
    std::vector<double> paramDirichletProportions;

    
public:
    BDCoalPriorParams( double lambdaExpPriorScaledGrowthRate, double                 lambdaExpPriorTheta,
                       double meanBetaPriorSeqError,
                       double varBetaPriorSeqError,
                       double meanBetaPriorADOError,
                       double varBetaPriorADOError,
                      double paramMultiplierProposalScaledGrowthRate,
                      double paramMultiplierProposalTheta);
    
//    getters
  
    ~BDCoalPriorParams(){};
};

#endif /* bd_coal_model_params_hpp */

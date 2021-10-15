//
//  bd_coal_model_params.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#include "bd_coal_model_params.hpp"
BDCoalPriorParams::BDCoalPriorParams(double lambdaExpPriorScaledGrowthRate, double                 lambdaExpPriorTheta,
                    double meanBetaPriorSeqError,
                    double varBetaPriorSeqError,
                    double meanBetaPriorADOError,
                    double varBetaPriorADOError,
                    double paramMultiplierProposalScaledGrowthRate,
                    double paramMultiplierProposalTheta)
:lambdaExpPriorScaledGrowthRate(lambdaExpPriorScaledGrowthRate), lambdaExpPriorTheta(lambdaExpPriorTheta),
meanBetaPriorSeqError(meanBetaPriorSeqError),
varBetaPriorSeqError(varBetaPriorSeqError),
meanBetaPriorADOError(meanBetaPriorADOError),
varBetaPriorADOError(varBetaPriorADOError),
paramMultiplierProposalScaledGrowthRate(paramMultiplierProposalScaledGrowthRate),
paramMultiplierProposalTheta(paramMultiplierProposalTheta)
{
    
    
};



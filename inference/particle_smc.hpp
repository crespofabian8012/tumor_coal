//
//  particle_smc.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 12/03/2021.
//

#ifndef particle_smc_hpp
#define particle_smc_hpp

#include <stdio.h>


#include "data_types.hpp"
#include "tree_node.hpp"
#include <boost/unordered_set.hpp>

extern "C"
{
#include <gsl/gsl_rng.h>
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}
#include <memory>

template <class S>
class ParticleSMC{
    
     S *s;//the State class
    long double logWeight;// log of weight update function

public:
    ParticleSMC(S *s);

    S *getState();
    long double getLogWeight();
};
template <class S>
ParticleSMC<S>::ParticleSMC(S *s )
{
     this->s = s;
    logWeight=0.0;
}

template <class S>
S *ParticleSMC<S>::getState()
{
    return s;
}
template <class S>
long double ParticleSMC<S>::getLogWeight()
{
    return logWeight;
}
long double ParticleSMC<S>::setLogWeight(double newlogWeight)
{
    logWeight =newlogWeight;
}

#endif /* particle_smc_hpp */



/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/

/*
 * random  functions
 */


#ifndef random_h
#define random_h


extern "C"
    {
#include <gsl/gsl_rng.h>
    }
#include <vector>

#include <boost/random.hpp>
namespace Random
    {
    long double    RandomUniform (long int *seed);
    long double  RandomExponential (long double  mean, long int *seed, bool useGSlgenerator,const gsl_rng *rngGsl,  boost::mt19937* rngBoost);
    long double     RandomGamma (long double  shape, long int *seed,  bool useGSlgenerator,const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double  RandomGamma1 (long double  s, long int *seed, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double  RandomGamma2 (long double  s, long int *seed, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    int RandomNegativeBinomial (long double  mean, long double  dispersion, long int *seed);
    long double  RandomNormal(long double  mean, long double  stddev, long int *seed);
    long double  RandomLogNormal( long double  mean, long double  dispersion, long int *seed);
    void  RandomDirichlet2 (long double  *vectorGamma, int vectorSize, long double  *outputVector, long int *seed);
    int RandomPoisson (long double  lambda, long int *seed,const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    int RandomBinomial (long double  prob, int numTrials, long int *seed);
    int RandomUniformTo (int max, long int *seed,bool useGSlgenerator,const gsl_rng* rngGsl, boost::mt19937* rngBoost);
    long double  RandomLogUniform( long double  from, long double  to,  bool useGSlgenerator, const gsl_rng * rngGsl,  boost::mt19937* rngBoost);
    void  RandomDirichlet (long double  s, int vectorSize, std::vector<double> &outputVector, long int *seed,const gsl_rng *randomGsl, boost::mt19937* rngBoost);
    int ChooseUniformState ( long double  *prob, long int *seed, bool useGSlgenerator, const gsl_rng *rngGsl, boost::mt19937* rngBoost);
    long double  randomUniformFromGsl();
    void   randomDirichletFromGsl(int vectorSize,  double  alpha[], long double  *theta);
    void  randomDirichletFromVector (std::vector<long double> alpha, std::vector<long double> &outputVector,bool useGSlgenerator, const gsl_rng *rngGsl, boost::mt19937* rngBoost);
    long double  randomGammaShapeScale( long double  shape, long double  scale, boost::mt19937* rng );
    long double  randomGammaMeanVariance( long double  mean , long double  variance, boost::mt19937* rng );
    long double randomNormalBoost(long double mean, long double sigma, boost::mt19937* rng);
    long double  RandomExponentialStartingFrom (long double  lambda, long double from, bool useGSlgenerator,const gsl_rng* rngGsl, boost::mt19937* rngBoost);
    long double  RandomPowerLawDistribution (long double  a, long double from, bool useGSlgenerator,const gsl_rng *randomGsl, boost::mt19937* rngBoost);
    long  double   randomNormalGreaterThan(long double mean, long double sigma, long double from, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double  RandomDensityModelTimeOrigin (long double  lambda, bool useGSlgenerator, long double from,   const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    int  randomUniformIntegerInterval(const gsl_rng * r, int from, int to);
    
    long double RandomDensityModelTimeOriginLambda (long double  lambda, bool useGSlgenerator, long double gamma,double populationSize, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double randomUniformBoost(boost::mt19937* rng);
    void allocateListRandomNumbersGenerators(std::vector<gsl_rng * > &vec);
    std::vector<boost::mt19937*> allocateListRandomNumbersGeneratorsBoost(int numberRngs);
    gsl_rng* generateRandomObject(long seed);
    void freeListRandomNumbersGenerators(std::vector<gsl_rng * > &vec);
    long double  randomUniformFromGsl2(const gsl_rng * r);
    long double  RandomExponentialExponential (long double  lambda, bool useGSlgenerator, long double from, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double randomBeta(long double  alpha, long double beta, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double randomGamma(long double  mean, long double variance, bool useGSlgenerator,const  gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    long double RandomBetaMeanVar(long double mean, long double var, long int* seed, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost);
    int randomDiscreteFromProbabilityVector( const gsl_rng *rngGsl,  double  *prob, int size);
    
    };
namespace Distributions{

long double LogExponentialDensity(long double& value, long double& lambda, long double& from);
long double  logMultinomialProbability(const  size_t size, const  double  *p, const unsigned int *n );
long double LogPowerLawDistibutionDensity(long double& value, long double& a, long double& from);
long double LogUniformDensity(long double& value, long double& from, long double& to);
long double  LogDirichletDensity(std::vector<long double> &proportionsVector,  std::vector<long double> &concentrationVector, std::vector<long double> &dummy);
long double LogBetaDensity(long double& value, long double& alpha, long double& beta);

};
#endif /* random_h */


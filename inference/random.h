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

#include <gsl/gsl_rng.h>
#include <vector>



using namespace std;
struct Random
{
static long double    RandomUniform (long int *seed);
static long double  RandomExponential (long double  mean, long int *seed, bool useGSlgenerator);
static long double     RandomGamma (long double  shape, long int *seed,  bool useGSlgenerator);
static long double  RandomGamma1 (long double  s, long int *seed, bool useGSlgenerator);
static long double  RandomGamma2 (long double  s, long int *seed, bool useGSlgenerator);
static int RandomNegativeBinomial (long double  mean, long double  dispersion, long int *seed);
static long double  RandomNormal(long double  mean, long double  stddev, long int *seed);
static long double  RandomLogNormal( long double  mean, long double  dispersion, long int *seed);
static void  RandomDirichlet2 (long double  *vectorGamma, int vectorSize, long double  *outputVector, long int *seed);
static int RandomPoisson (long double  lambda, long int *seed);

static int RandomBinomial (long double  prob, int numTrials, long int *seed);

static int RandomUniformTo (int max, long int *seed,bool useGSlgenerator);

static  long double  RandomLogUniform( long double  from, long double  to,  bool useGSlgenerator, gsl_rng * randomGenerator);

static void  RandomDirichlet (long double  s, int vectorSize, std::vector<double> &outputVector, long int *seed);


static int ChooseUniformState ( long double  *prob, long int *seed, bool useGSlgenerator);

static long double  randomUniformFromGsl();
static void   randomDirichletFromGsl(int vectorSize,  double  alpha[], long double  *theta);

static void  randomDirichletFromVector (std::vector<long double> alpha, std::vector<long double> &outputVector);
static long double  randomGammaBoost( long double  shape, long double  scale );
static long double  randomGammaBoost2( long double  mean , long double  variance );

static long double randomNormalBoost(long double mean, long double sigma);
static long double  RandomExponentialStartingFrom (long double  lambda, long double from, bool useGSlgenerator, gsl_rng* randomGenerator);
static long double  RandomPowerLawDistribution (long double  a, long double from, bool useGSlgenerator);
static long  double   randomNormalGreaterThan(long double mean, long double sigma, long double from);
static long double  RandomDensityModelTimeOrigin (long double  lambda, bool useGSlgenerator, long double from,  gsl_rng * randomGenerator);

static long double randomUniformBoost();
static void allocateListRandomNumbersGenerators(vector<gsl_rng * > &vec);
static void freeListRandomNumbersGenerators(vector<gsl_rng * > &vec);
static long double  randomUniformFromGsl2(gsl_rng * r);
    static long double  RandomExponentialExponential (long double  lambda, bool useGSlgenerator, long double from);
};
struct Distributions{
static long double LogExponentialDensity(long double lambda, long double value, long double from);
static long double  logMultinomialProbability(const  size_t size, const  double  *p, const unsigned int *n );
static double LogUniformDensity(double value, double from, double to);
static long double LogPowerLawDistibutionDensity(long double a, long double value, long double from);
    static long double LogUniformDensity(long double value, long double from, long double to);
};
#endif /* random_h */


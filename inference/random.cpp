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
 * random functions
 */
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#define BOOST_PENDING_INTEGER_LOG2_HPP

#include "random.h"


#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <sys/time.h>

//#include <chrono>

#include <boost/integer/integer_log2.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/nondet_random.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/function.hpp>
#include <boost/random.hpp>
#include <boost/multiprecision/random.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/number.hpp>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>



using namespace boost::multiprecision;
using namespace boost::random;
using namespace boost::math;

typedef boost::random::independent_bits_engine<boost::random::mt19937, 256, boost::multiprecision::cpp_int> generator_type;

class MCMCoptions;
/***************************** RandomUniform **********************************/
/* It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller. 1988. Random number generators: good
 ones are hard to find. Communications of the ACM, 31(10):1192-1201.
 */

long double  Random::RandomUniform (long int *seed)
{
    long int  lo, hi, test;
    
    hi = (*seed) / 127773;
    lo = (*seed) % 127773;
    test = 16807 * lo - 2836 * hi;
    if (test > 0)
        *seed = test;
    else
        *seed = test + 2147483647;
    return (double)(*seed) / (double)2147483647;
}

/**************************** RandomGamma *************************/
/*    Generates a gamma number using routines in Ziheng's
 Yang tools.h in PAML
 
 Random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
 r^(s-1)*exp(-r)
 
 J. Dagpunar (1988) Principles of random variate generation,
 Clarendon Press, Oxford
 
 Calling rndgamma1() if s<1 or rndgamma2() if s>1 or exponential if s=1
 */

long double     Random::RandomGamma (long double  shape, long int *seed, bool useGSlgenerator,const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    long double  gammaNumber = 0;
    
    if (shape <= 0)
        fprintf (stderr, "ERROR: problems with gamma variable generation, shape < 0");
    else if (shape < 1)
        gammaNumber = Random::RandomGamma1 (shape, seed, useGSlgenerator, rngGsl,   rngBoost );
    else if (shape > 1)
        gammaNumber = Random::RandomGamma2 (shape, seed, useGSlgenerator, rngGsl,   rngBoost);
    else
        gammaNumber = -log (Random::randomUniformFromGsl2(rngGsl));
    return (gammaNumber);
}
/*************** RandomGamma1 ***************/
long double  Random::RandomGamma1 (long double  s, long int *seed, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    /* Random standard gamma for s<1
     switching method
     */
    long double             r, x=0.0, small=1e-37, w;
    static long double    a, p, uf, ss=10.0, d;
    
    if (s!=ss)
    {
        a  = 1.0-s;
        p  = a/(a+s*exp(-a));
        uf = p*pow(small/a,s);
        d  = a*log(a);
        ss = s;
    }
    for (;;)
    {
        if (useGSlgenerator)
            r = randomUniformFromGsl2(rngGsl);
        else
            r = randomUniformBoost(rngBoost);
        if (r > p)
        {
            x = a-log((1.0-r)/(1.0-p));
            w=a*log(x)-d;  /* this was with comma in line above before 270917*/
        }
        else if (r>uf)
        {
            x = a*pow(r/p,1/s);
            w=x; /* this was with comma in line above before 270917*/
        }
        else
            return (0.0);
        if (useGSlgenerator)
            r = randomUniformFromGsl2(rngGsl);
        else
            r = randomUniformBoost(rngBoost);
        if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
        break;
    }
    return (x);
}
/*************** RandomGamma2 ***************/
long double  Random::RandomGamma2 (long double  s, long int *seed, bool useGSlgenerator, const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
    /* Random standard gamma for s>1
     Best's (1978) t distribution method
     */
    long double             r ,d, f, g, x;
    static long double     b, h, ss=0;
    
    if (s!=ss)
    {
        b  = s-1.0;
        h  = sqrt(3.0*s-0.75);
        ss = s;
    }
    for (;;)
    {
        if (useGSlgenerator)
            r = randomUniformFromGsl2(randomGsl);
        else
            r = randomUniformBoost(rngBoost);
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        if (useGSlgenerator)
            r = randomUniformFromGsl2(randomGsl);
        else
            r = randomUniformBoost(rngBoost);
        d = 64*r*r*g*g*g;
        if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
            break;
    }
    return (x);
}

/********************* RandomPoisson ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

int Random::RandomPoisson (long double  lambda, long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    int        poissonNumber;
    long double     sum;
    
    sum = 0;
    poissonNumber = -1;
    
    while (sum <= 1.0)
    {
        sum += RandomExponential (lambda, seed, true, rngGsl,  rngBoost);
        poissonNumber++;
    }
    
    return poissonNumber;
}
/************************************************************/
/********************* RandomExponential ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda.
 */

long double  Random::RandomExponential (long double  lambda, long int *seed, bool useGSlgenerator,const gsl_rng *randomGsl,  boost::mt19937* rngBoost)
{
    long double   exponentialNumber, U;
    
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl2(randomGsl);
        else
           U = randomUniformBoost(rngBoost);
       }
    while (U == 0 || U == 1);
    
    exponentialNumber = -log (U) / lambda;
    
    return exponentialNumber;
}
/********************* RandomExponentialStartingFrom ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda and with the 0 swifted to the  value from.
 */
long double  Random::RandomExponentialStartingFrom (long double  lambda, long double from, bool useGSlgenerator,const gsl_rng *randomGsl,  boost::mt19937* rngBoost)
{
    //from  have to be non negative
    assert(from>=0);
    long double   exponentialNumber, U;
    
    do{
         if (useGSlgenerator)
           U = randomUniformFromGsl2(randomGsl);
        else
            U = randomUniformBoost(rngBoost);
    }
    while (U == 1 || U == 0);
    
 //  exponentialNumber = -log ( (1 - U)/ exp(lambda* from)) / lambda;
    exponentialNumber = from - log (1-U) / lambda;

    return exponentialNumber;
}
/********************* RandomPowerLawDistribution ********************/
/* Generates a random number from a RandomPowerLaw distibution with
 parametrr a>0 and with the 0 swifted to the  value from.
 */
long double  Random::RandomPowerLawDistribution (long double  a, long double from, bool useGSlgenerator,const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
    //from value have to be positive
    long double   output, U;
    
    if (useGSlgenerator)
        U = randomUniformFromGsl2(randomGsl);
    else
        U = randomUniformBoost(rngBoost);
   
    output =from + U / (a*(1-U));
    return output;
}
/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

int Random::RandomUniformTo (int max, long int *seed,  bool useGSlgenerator, const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
    long double     U;
    //rd = randomUniformFromGsl();
    if (useGSlgenerator)
        U = randomUniformFromGsl2(randomGsl);
    else
        U = randomUniformBoost(rngBoost);
    return (floor(U*max));
}

/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */
int Random::ChooseUniformState (long double  *prob, long int *seed,  bool useGSlgenerator, const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
    int            chosenState;
    long double         U, cumProb;
    
    chosenState = 0;
    cumProb = prob[chosenState];
    // ran = RandomUniform(seed);
    
    if (useGSlgenerator)
        U = randomUniformFromGsl2(randomGsl);
    else
        U = randomUniformBoost(rngBoost);
    
    while (U > cumProb)
        cumProb += prob[++chosenState];
    
    return chosenState;
}


/*************** RandomDirichlet ***************/
//first generates random samples from a gamma and then divide each value by the sum
void  Random::RandomDirichlet (long double  s, int vectorSize, std::vector<double> &outputVector, long int *seed,const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{   int i;
    long double  sum=0.0;
    long double  current;
    // *outputVector = malloc(vectorSize * sizeof(double));
    // if (*outputVector == NULL)
    //     return;
    for (i=0; i < vectorSize; i++){
        current = RandomGamma(s, seed, true, randomGsl, rngBoost);
        outputVector[i] = current;
        //outputVector[i] = current;
        //*(*outputVector + i)=current;
        sum=sum+current;
    }
    for (i=0; i < vectorSize; i++){
        outputVector[i]=  outputVector[i]/ sum;
        // *(*outputVector + i)= *(*outputVector + i)/sum;
    }
}
/*************** RandomDirichlet ***************/
//first generates random samples from a gamma and then divide each value by the sum
void  Random::randomDirichletFromVector (std::vector<long double> alpha, std::vector<long double> &outputVector,bool useGSlgenerator,const gsl_rng *rngGsl, boost::mt19937* rngBoost)
{   int i;
    long double sum=0.0;
    long double current;
    //long int seed =0;
    // *outputVector = malloc(vectorSize * sizeof(double));
    // if (*outputVector == NULL)
    //     return;
    outputVector.clear();
    
    for (i=0; i < alpha.size(); i++){
        
        if (useGSlgenerator)
            current = RandomGamma (alpha.at(i), NULL, useGSlgenerator, rngGsl,  rngBoost);
            //current = RandomGamma(alpha.at(i), &seed);
        else
           current = randomGammaShapeScale(alpha.at(i), 1, rngBoost);
        outputVector.push_back(current);
        sum=sum+current;
    }
    for (i=0; i < alpha.size(); i++){
        outputVector.at(i)=  outputVector.at(i)/ sum;
        // *(*outputVector + i)= *(*outputVector + i)/sum;
    }
}
void   Random::randomDirichletFromGsl(int vectorSize,  double  alpha[], long double  *theta)
{
    //long double  thetaDouble[vectorSize];
    double *thetaDouble;
    thetaDouble = new double[vectorSize];
    for (unsigned int i = 0; i < vectorSize; ++i){
        thetaDouble[i]=0;
        theta[i]=0;
    }
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
   
    //T = gsl_rng_ranlux389; // Generator setup
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, mySeed);
    //unsigned long mySeed2 =gsl_rng_uniform_int( r, mySeed);
    //unsigned long mySeed2 = ((int)randomUniformFromGsl() )* mySeed;
    //gsl_rng_set(r, mySeed2);
    
    //gsl_ran_dirichlet( r,  vectorSize, alpha, &thetaDouble[0]);
    gsl_ran_dirichlet( r,  vectorSize, alpha, thetaDouble);
    for (unsigned int i = 0; i < vectorSize; ++i){
        theta[i]=thetaDouble[i];
       }
    gsl_rng_free (r);
}

long double  Random::RandomLogUniform( long double  from, long double  to,  bool useGSlgenerator, const gsl_rng * rngGsl,  boost::mt19937* rngBoost){
    long double U;
    
       if (useGSlgenerator)
            U = randomUniformFromGsl2(rngGsl);
        else
            U = randomUniformBoost(rngBoost);
        
    return(exp(from + U*(to -from)));
}
// from https://stackoverflow.com/questions/9768519/gsl-uniform-random-number-generator
//long double  Random::randomUniformFromGsl(){
//    const gsl_rng_type * T;
//   gsl_rng * r;
//    gsl_rng_env_setup();
//    struct timeval tv; // Seed generation based on time
//   gettimeofday(&tv,0);
//    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
//    T= gsl_rng_ranmar;
//     r = gsl_rng_alloc (T);
//   // printf ("seed = %lu\n", gsl_rng_default_seed);
//    gsl_rng_set(r, mySeed);
//    long double  u = gsl_rng_uniform(r); // Generate it!
//    gsl_rng_free (r);
//    return u;
//    
//}
long double  Random::randomUniformFromGsl2(const gsl_rng * r){

    long double  u = gsl_rng_uniform(r);
    return u;
    
}
int  Random::randomUniformIntegerInterval(const gsl_rng * r, int from, int to){
   
    assert(to>=from);
    long int  u;
    if (to>from)
     u = gsl_rng_uniform_int(r, to-from);
    else
        u=0;
    return (from +u);
    
}
void Random::allocateListRandomNumbersGenerators(std::vector< gsl_rng * > &vec)
{
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    unsigned long mySeed;
    for(unsigned int i=0; i < vec.size(); i++)
    {
        gettimeofday(&tv,0);
        mySeed = tv.tv_sec + tv.tv_usec;
        //vec.at(i)= generateRandomObject(mySeed);
        T= gsl_rng_ranmar;
        vec.at(i) = gsl_rng_alloc (T);
        gsl_rng_set(vec.at(i), mySeed);
        
    }
}
gsl_rng* Random::generateRandomObject(long seed)
{
    //const gsl_rng_type * T = gsl_rng_ranmar;
    const gsl_rng_type * T = gsl_rng_taus;
    gsl_rng *random;

    gsl_rng_env_setup();

    random = gsl_rng_alloc(T);
    gsl_rng_set(random, seed);
    return random;
    
}
std::vector<boost::mt19937*> Random::allocateListRandomNumbersGeneratorsBoost(int numberRngs)
{
    std::vector<boost::mt19937 *> vec;
    struct timeval tv; // Seed generation based on time
    unsigned long mySeed;
    for(unsigned int i=0; i < numberRngs; i++)
    {
        gettimeofday(&tv,0);
        mySeed = tv.tv_sec + tv.tv_usec;
        boost::mt19937 rng;
        rng.seed(mySeed);
        vec.push_back( &rng);
        
    }
    return vec;
}
void Random::freeListRandomNumbersGenerators(std::vector<gsl_rng * > &vec)
{
    for(unsigned int i=0; i < vec.size(); i++)
    {
      gsl_rng_free (vec.at(i));
        
    }
}

long double  Random::randomGammaShapeScale( long double  shape, long double  scale, boost::mt19937* rng ) {
    //    boost::random_device rd;
    //    std::uint64_t value = rd();
    //value = (value << 32) | rd();
    //    boost::gamma_distribution<std::uint64_t> dis;
    //    boost::function<std::uint64_t()> gen = boost::bind(dis, boost::ref(rd));
    //    std::uint64_t value = gen();
    //boost::mt19937 rng(56) ;
    boost::gamma_distribution<> gd( shape );
    boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( *rng, gd );
    return scale*var_gamma();
}
long double  Random::randomGammaMeanVariance( long double  mean , long double  variance, boost::mt19937* rng )
{
    const long double  shape = ( mean * mean )/variance;
    long double  scale = variance/mean;
    //boost::mt19937 rng(56) ;
    
    boost::gamma_distribution<> gd( shape );
    boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( *rng, gd );
    return scale*var_gamma();
}
long double  Distributions::logMultinomialProbability(const  size_t size, const  double  *p, const unsigned int *n )
{
    return gsl_ran_multinomial_lnpdf(size, p, n);
}
long  double   Random::randomNormalBoost(long double mean, long double sigma, boost::mt19937* rng){

    boost::normal_distribution< boost::multiprecision::cpp_dec_float_50> normal_dist(0, sigma);
    boost::random::independent_bits_engine<boost::mt19937, 50L * 1000L / 301L, boost::multiprecision::number<boost::multiprecision::cpp_int::backend_type, boost::multiprecision::et_off> > gen;
    
    long double result = mean +(long double)normal_dist(gen);
    return result;
}
long  double   Random::randomNormalGreaterThan(long double mean, long double sigma, long double from,bool useGSlgenerator,const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    long  double result;
    if (useGSlgenerator){
        result = mean + gsl_ran_gaussian(rngGsl, sigma);
    }
    else{
        result= randomNormalBoost(mean, sigma, rngBoost);
        }
    
    
    while(result < from)
    {
        result = from + (from-result);
    }
    return result;
}
/********************* RandomDensityModelTimeOrigin ********************/
long double  Random::RandomDensityModelTimeOrigin (long double  lambda, bool useGSlgenerator, long double from, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    //from value have to be positive
    long double   result, U;
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl2(rngGsl);
        else
            U = randomUniformBoost(rngBoost);
    }
    while (U == 0);
    
 
    result = from +  (1.0 / lambda) * log(1- lambda / (log(U)));

    return result;
}
/********************* RandomDensityModelTimeOriginLambda ********************/
long double  Random::RandomDensityModelTimeOriginLambda (long double  lambda, bool useGSlgenerator, long double gamma,double populationSize, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    //from value have to be positive
    long double   result, U;
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl2(rngGsl);
        else
            U = randomUniformBoost(rngBoost);
    }
    while (U == 0);
    
    if (useGSlgenerator)
        U= gsl_ran_gamma(rngGsl, populationSize, 1.0);
    else
        U= randomGammaShapeScale(populationSize, 1.0, rngBoost );
    //result = from +  (1.0 / lambda) * log(1- lambda / (log(U)));
    result =    -(1.0 / (gamma)) * log(1- gamma / (U+gamma));
    return result;
}
/***************************** LogExponentialDensity *******************************/
/* Log Exponential distribution density */
long double Distributions::LogExponentialDensity(long double& value, long double& from, long double& lambda)
{
    double result=0.0;
    
    if (value >= from )
    {
        result = log(lambda)-lambda*(value -from);
        
        auto const lambda1 = lambda;
        const auto d = boost::math::exponential_distribution<>{static_cast<double>(lambda1)};
        result=  log(boost::math::pdf(d, value-from));
      
        return result;
    }
    else
    {
        result=log(0);
        return result;
    }
    
}
/***************************** LogBetaDensity *******************************/
/* Log of Beta distribution density */
long double Distributions::LogBetaDensity(long double& value, long double& alpha, long double& beta)
{
    double result=0.0;
    
    if (value >= 0 && value <= 1 )
    {
        auto const alpha1 = alpha;
        auto const beta1 = beta;
        const auto d = boost::math::beta_distribution<>(alpha1, beta1);
        result=  log(boost::math::pdf(d, value));
      
        return result;
    }
    else
    {
        result=log(0);
        return result;
    }
    
}
long double Random::randomBeta(long double  alpha, long double beta, bool useGSlgenerator,const gsl_rng *rngGsl, boost::random::mt19937 * rngBoost)
{
    double result=0.0;
    assert(alpha>0);
    assert(beta>0);
    long double    U;
       
       do{
            if (useGSlgenerator)
              U = randomUniformFromGsl2(rngGsl);
           else
               U = randomUniformBoost(rngBoost);
       }
       while (U == 1 || U == 0);
       
      auto const alpha1 = alpha;
      auto const beta1 = beta;
      const auto dist = boost::math::beta_distribution<>(alpha1, beta1);
     
      result = quantile(dist, U);

       return result;
    
}
long double Random::randomGamma(long double  mean, long double variance, bool useGSlgenerator,const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    double result=0.0;
    const double shape = ( mean*mean )/variance;
    double scale = variance/mean;
    long double    U;
       
       do{
            if (useGSlgenerator)
              U = randomUniformFromGsl2(rngGsl);
           else
               U = randomUniformBoost(rngBoost);
       }
       while (U == 1 || U == 0);
       
   
      const auto dist = boost::gamma_distribution<>(shape);
      boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma(*rngBoost, dist );
     
      //result = quantile(dist, U);
      result=  scale*var_gamma();

       return result;
}
long double Random::randomUniformBoost(boost::random::mt19937 * rngBoost)
{
    
 boost::random::independent_bits_engine<boost::random::mt19937, 50L*1000L/301L, cpp_int> gen;
     std::setprecision(50);
    
    long double result = (long double)boost::random::generate_canonical<cpp_bin_float_50, std::numeric_limits<cpp_bin_float_50>::digits>(gen);
    
    std:: cout << "\n random boost: "   <<  result << std::endl;
    return result;
    
}
/********************* RandomExponentialExponential********************/
long double  Random::RandomExponentialExponential (long double  lambda, bool useGSlgenerator, long double from, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
{
    //from parameter have to be positive
    long double   result, U;
    
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl2(rngGsl);
        else
            U = randomUniformBoost(rngBoost);
    }
    while (U == 1);
    
    result = from +  lambda* U /(1-U);
    return result;
}
long double Random::RandomBetaMeanVar(long double mean, long double var, long int* seed, bool useGSlgenerator, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost){
    
    long double sample_size, shape1, shape2, gamma1, gamma2, randBeta;
    
    /* assuming variance < mean (1-mean) */
    sample_size = (mean * (1.0 - mean) / var) - 1.0;
    
    shape1 = mean * sample_size;
    shape2 = (1.0 - mean) * sample_size;
    
    gamma1 = RandomGamma(shape1, seed,  useGSlgenerator,  rngGsl,  rngBoost);
    gamma2 = RandomGamma(shape2, seed, useGSlgenerator,  rngGsl,  rngBoost);
    randBeta = gamma1 / (gamma1 + gamma2);
    
    return randBeta;
    
    
}
/***************************** LogUniformDensity *******************************/
/* Log Uniform distribution density */
long double Distributions::LogUniformDensity(long double& value,  long double& from, long double& to)
{
    double result;
    if (value >= exp(from) && value<= exp(to) )
    {
        result = -log(value) - log(to-from);
        return result;
    }
    else
    {
        result=log(0);
        return result;
    }
}
/***************************** LogPowerLawDistibutionDensity *******************************/
/* Log Power Law distribution density */
long double Distributions::LogPowerLawDistibutionDensity( long double& value, long double& a,  long double& from){
    
    if (value < from)
        return log(0);
    
    long double firstTerm=log(a);
    long double secondTerm=2*log(1+a*(value-from));
    //long double fullTerm= a / (1 + a*(value-from))*(1 + a*(value-from));
    long double logfullTerm= firstTerm -secondTerm;
    long double result=logfullTerm;
    return result;
}

long double  Distributions::LogDirichletDensity(std::vector<long double> &proportionsVector,  std::vector<long double> &concentrationVector,  std::vector<long double> &dummy)
{
    int i;
    long double  sum=0;
    long double  logResult=0;
    for( i = 0 ; i < proportionsVector.size(); i++)
    {
        sum = sum +concentrationVector[i];
        logResult= logResult+(concentrationVector[i]-1)*log(proportionsVector[i]);
        logResult= logResult+ lgamma(concentrationVector[i]);
    }
    logResult = logResult-lgamma(sum);
    return logResult;
}

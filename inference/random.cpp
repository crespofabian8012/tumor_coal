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
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <sys/time.h>
#include <chrono>
#include <ctime>
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
#include <boost/multiprecision/cpp_bin_float.hpp>



using namespace boost::multiprecision;
using namespace boost::random;
using namespace boost::math;

typedef boost::random::independent_bits_engine<boost::random::mt19937, 256, boost::multiprecision::cpp_int> generator_type;


/***************************** RandomUniform **********************************/
/* It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller. 1988. Random number generators: good
 ones are hard to find. Communications of the ACM, 31(10):1192-1201.
 */

long double  RandomUniform (long int *seed)
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

long double     RandomGamma (long double  shape, long int *seed, bool useGSlgenerator)
{
    long double  gammaNumber = 0;
    
    if (shape <= 0)
        fprintf (stderr, "ERROR: problems with gamma variable generation, shape < 0");
    else if (shape < 1)
        gammaNumber = RandomGamma1 (shape, seed, useGSlgenerator );
    else if (shape > 1)
        gammaNumber = RandomGamma2 (shape, seed, useGSlgenerator);
    else
        gammaNumber = -log (randomUniformFromGsl());
    return (gammaNumber);
}
/*************** RandomGamma1 ***************/
long double  RandomGamma1 (long double  s, long int *seed, bool useGSlgenerator)
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
            r = randomUniformFromGsl();
        else
            r = randomUniformBoost();
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
            r = randomUniformFromGsl();
        else
            r = randomUniformBoost();
        if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
        break;
    }
    return (x);
}
/*************** RandomGamma2 ***************/
long double  RandomGamma2 (long double  s, long int *seed, bool useGSlgenerator)
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
            r = randomUniformFromGsl();
        else
            r = randomUniformBoost();
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        if (useGSlgenerator)
            r = randomUniformFromGsl();
        else
            r = randomUniformBoost();
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

int RandomPoisson (long double  lambda, long int *seed)
{
    int        poissonNumber;
    long double     sum;
    
    sum = 0;
    poissonNumber = -1;
    
    while (sum <= 1.0)
    {
        sum += RandomExponential (lambda, seed, true);
        poissonNumber++;
    }
    
    return poissonNumber;
}
/************************************************************/
/********************* RandomExponential ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda.
 */

long double  RandomExponential (long double  lambda, long int *seed, bool useGSlgenerator)
{
    long double   exponentialNumber, U;
    
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl();
        else
           U = randomUniformBoost();
       }
    while (U == 0);
    
    exponentialNumber = -log (U) / lambda;
    
    return exponentialNumber;
}
/********************* RandomExponentialStartingFrom ********************/
/* Generates a random number from a Exponential distibution with
 mean lambda and with the 0 swifted to the  value from.
 */
long double  RandomExponentialStartingFrom (long double  lambda, long double from, bool useGSlgenerator)
{
    //from value have to be positive
    long double   exponentialNumber, U;
    
    do{
         if (useGSlgenerator)
           U = randomUniformFromGsl();
        else
            U = randomUniformBoost();
    }
    while (U == 0);
    
 //  exponentialNumber = -log ( (1 - U)/ exp(lambda* from)) / lambda;
    exponentialNumber = from - log (U) / lambda;

    return exponentialNumber;
}
/********************* RandomPowerLawDistribution ********************/
/* Generates a random number from a RandomPowerLaw distibution with
 parametrr a>0 and with the 0 swifted to the  value from.
 */
long double  RandomPowerLawDistribution (long double  a, long double from, bool useGSlgenerator)
{
    //from value have to be positive
    long double   output, U;
    
    if (useGSlgenerator)
        U = randomUniformFromGsl();
    else
        U = randomUniformBoost();
   
    output =from + U / (a*(1-U));
    return output;
}
/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

int RandomUniformTo (int max, long int *seed,  bool useGSlgenerator)
{
    long double     U;
    //rd = randomUniformFromGsl();
    if (useGSlgenerator)
        U = randomUniformFromGsl();
    else
        U = randomUniformBoost();
    return (floor(U*max));
}

/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */
int ChooseUniformState (long double  *prob, long int *seed,  bool useGSlgenerator)
{
    int            chosenState;
    long double         U, cumProb;
    
    chosenState = 0;
    cumProb = prob[chosenState];
    // ran = RandomUniform(seed);
    
    if (useGSlgenerator)
        U = randomUniformFromGsl();
    else
        U = randomUniformBoost();
    
    while (U > cumProb)
        cumProb += prob[++chosenState];
    
    return chosenState;
}


/*************** RandomDirichlet ***************/
//first generates random samples from a gamma and then divide each value by the sum
void  RandomDirichlet (long double  s, int vectorSize, vector<double> &outputVector, long int *seed)
{   int i;
    long double  sum=0.0;
    long double  current;
    // *outputVector = malloc(vectorSize * sizeof(double));
    // if (*outputVector == NULL)
    //     return;
    for (i=0; i < vectorSize; i++){
        current = RandomGamma(s, seed, true);
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
void  randomDirichletFromVector (vector<long double> alpha, vector<long double> &outputVector)
{   int i;
    long double sum=0.0;
    long double current;
    //long int seed =0;
    // *outputVector = malloc(vectorSize * sizeof(double));
    // if (*outputVector == NULL)
    //     return;
    outputVector.clear();
    
    for (i=0; i < alpha.size(); i++){
        
        //current = RandomGamma(alpha.at(i), &seed);
        current = randomGammaBoost(alpha.at(i), 1);
        outputVector.push_back(current);
        sum=sum+current;
    }
    for (i=0; i < alpha.size(); i++){
        outputVector.at(i)=  outputVector.at(i)/ sum;
        // *(*outputVector + i)= *(*outputVector + i)/sum;
    }
}
void   randomDirichletFromGsl(int vectorSize,  double  alpha[], long double  *theta)
{
    long double  thetaDouble[vectorSize];
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
    //mySeed * randomUniformFromGsl();
    T = gsl_rng_ranlux389; // Generator setup
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, mySeed);
    //unsigned long mySeed2 =gsl_rng_uniform_int( r, mySeed);
    //unsigned long mySeed2 = ((int)randomUniformFromGsl() )* mySeed;
    //gsl_rng_set(r, mySeed2);
    gsl_ran_dirichlet( r,  vectorSize, alpha, &thetaDouble[0]); // Generate it!
    for (unsigned int i = 0; i < vectorSize; ++i){
        theta[i]=thetaDouble[i];
       }
    gsl_rng_free (r);
}
long double  RandomLogUniform( long double  from, long double  to,  bool useGSlgenerator){
    long double U;
    
       if (useGSlgenerator)
            U = randomUniformFromGsl();
        else
            U = randomUniformBoost();
        
    return(exp(from + U*(to -from)));
}
// from https://stackoverflow.com/questions/9768519/gsl-uniform-random-number-generator
long double  randomUniformFromGsl(){
    const gsl_rng_type * T;
   gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
   gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
//    // mySeed =125673899;
    T= gsl_rng_ranmar;
     r = gsl_rng_alloc (T);
   // printf ("seed = %lu\n", gsl_rng_default_seed);
    gsl_rng_set(r, mySeed);
    long double  u = gsl_rng_uniform(r); // Generate it!
    gsl_rng_free (r);
    //fprintf (stderr, "\n random gsl:  %.10Lf \n", u);
    return u;
    
}
long double  randomUniformFromGsl2(gsl_rng * r){
    //    const gsl_rng_type * T;
    //    gsl_rng * r;
    //    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
     gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    //    // mySeed =125673899;
    //    T= gsl_rng_ranmar;
    //    r = gsl_rng_alloc (T);
    // printf ("seed = %lu\n", gsl_rng_default_seed);
    gsl_rng_set(r, mySeed);
    long double  u = gsl_rng_uniform(r); // Generate it!
    gsl_rng_free (r);
    //fprintf (stderr, "\n random gsl:  %.10Lf \n", u);
    return u;
    
}
long double  randomGammaBoost( long double  shape, long double  scale ) {
    //    boost::random_device rd;
    //    std::uint64_t value = rd();
    //value = (value << 32) | rd();
    //    boost::gamma_distribution<std::uint64_t> dis;
    //    boost::function<std::uint64_t()> gen = boost::bind(dis, boost::ref(rd));
    //    std::uint64_t value = gen();
    boost::mt19937 rng(56) ;
    boost::gamma_distribution<> gd( shape );
    boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( rng, gd );
    return scale*var_gamma();
}
long double  randomGammaBoost2( long double  mean , long double  variance )
{
    const long double  shape = ( mean * mean )/variance;
    long double  scale = variance/mean;
    boost::mt19937 rng(56) ;
    
    boost::gamma_distribution<> gd( shape );
    boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma( rng, gd );
    return scale*var_gamma();
}
long double  logMultinomialProbability(const  size_t size, const  double  *p, const unsigned int *n )
{
    return gsl_ran_multinomial_lnpdf(size, p, n);
}
long  double   randomNormalBoost(long double mean, long double sigma){
    
    boost::mt19937 rng(57);
    
    boost::normal_distribution< boost::multiprecision::cpp_dec_float_50> normal_dist(0, sigma);
    
boost::random::independent_bits_engine<boost::mt19937, 50L * 1000L / 301L, boost::multiprecision::number<boost::multiprecision::cpp_int::backend_type, boost::multiprecision::et_off> > gen;
    
    long double result = mean +(long double)normal_dist(gen);
    return result;
}
long  double   randomNormalGreaterThan(long double mean, long double sigma, long double from)
{
    long  double result= randomNormalBoost(mean, sigma);
    
    while(result < from)
    {
        result = from + (from-result);
    }
    return result;
}
/********************* RandomDensityModelTimeOrigin ********************/
long double  RandomDensityModelTimeOrigin (long double  lambda, bool useGSlgenerator, long double from)
{
    //from value have to be positive
    long double   result, U;
    
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl();
        else
            U = randomUniformBoost();
    }
    while (U == 0);

    result = from +  (1.0 / lambda) * log(1- lambda / (log(U)));
    return result;
}
/***************************** LogExponentialDensity *******************************/
/* Log Exponential distribution density */
long double LogExponentialDensity(long double lambda, long double value, long double from)
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
long double randomUniformBoost()
{
    
 boost::random::independent_bits_engine<boost::random::mt19937, 50L*1000L/301L, cpp_int> gen;
     std::setprecision(50);
    
    long double result = (long double)boost::random::generate_canonical<cpp_bin_float_50, std::numeric_limits<cpp_bin_float_50>::digits>(gen);
    
    fprintf (stderr, "\n random boost:  %.10Lf \n", result);
    return result;
    
}
/********************* RandomExponentialExponential********************/
long double  RandomExponentialExponential (long double  lambda, bool useGSlgenerator, long double from)
{
    //from value have to be positive
    long double   result, U;
    
    do{
        if (useGSlgenerator)
            U = randomUniformFromGsl();
        else
            U = randomUniformBoost();
    }
    while (U == 1);
    
    result = from +  lambda* U /(1-U);
    return result;
}

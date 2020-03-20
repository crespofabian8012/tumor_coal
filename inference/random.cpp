//
//  random.cpp
//  simul
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//
#include "random.h"
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <sys/time.h>
#include <chrono>
#include <ctime>
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

long double     RandomGamma (long double  shape, long int *seed)
{
    long double  gammaNumber = 0;
    
    if (shape <= 0)
        fprintf (stderr, "ERROR: problems with gamma variable generation, shape < 0");
    else if (shape < 1)
        gammaNumber = RandomGamma1 (shape, seed);
    else if (shape > 1)
        gammaNumber = RandomGamma2 (shape, seed);
    else
        gammaNumber = -log (randomUniformFromGsl());
    return (gammaNumber);
}
/*************** RandomGamma1 ***************/
long double  RandomGamma1 (long double  s, long int *seed)
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
        r = randomUniformFromGsl();
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
        r = randomUniformFromGsl();
        if (1.0-r <= w && r > 0.0)
            if (r*(w+1.0) >= 1.0 || -log(r) <= w)
                continue;
        break;
    }
    return (x);
}
/*************** RandomGamma2 ***************/
long double  RandomGamma2 (long double  s, long int *seed)
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
        r = randomUniformFromGsl();
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        r = randomUniformFromGsl();
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
        sum += RandomExponential (lambda, seed);
        poissonNumber++;
    }
    
    return poissonNumber;
}
/************************************************************/
/********************* RandomExponential ********************/
/* Generates a random number from a Poisson distibution with
 mean lambda.
 */

long double  RandomExponential (long double  lambda, long int *seed)
{
    long double   exponentialNumber, U;
    
    do
        U = randomUniformFromGsl();
    while (U == 0);
    
    exponentialNumber = -log (U) / lambda;
    
    return exponentialNumber;
}
long double  RandomExponentialStartingFrom (long double  lambda, long double from)
{
    //from value have to be positive
    long double   exponentialNumber, U;
    
    do
        U = randomUniformFromGsl();
    while (U == 0);
    
    exponentialNumber = -log (1 - U / exp(lambda* from)) / lambda;
    
    return exponentialNumber;
}

/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

int RandomUniformTo (int max, long int *seed)
{
    long double     rd;
    rd = randomUniformFromGsl();
    return (floor(rd*max));
}

/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */
int ChooseUniformState (long double  *prob, long int *seed)
{
    int            chosenState;
    long double         ran, cumProb;
    
    chosenState = 0;
    cumProb = prob[chosenState];
    ran = RandomUniform(seed);
    
    while (ran > cumProb)
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
        current = RandomGamma(s, seed);
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
    long int seed =0;
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
    //long double  theta[vectorSize];
    for (unsigned int i = 0; i < vectorSize; ++i){
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
    //gsl_rng_set(r, mySeed);
    //unsigned long mySeed2 =gsl_rng_uniform_int( r, mySeed);
    //unsigned long mySeed2 = ((int)randomUniformFromGsl() )* mySeed;
    //gsl_rng_set(r, mySeed2);
    gsl_ran_dirichlet( r,  vectorSize, alpha, theta); // Generate it!
    gsl_rng_free (r);
}
long double  RandomLogUniform( long double  from, long double  to){
    return(exp(from + randomUniformFromGsl()*(to -from)));
}
// from https://stackoverflow.com/questions/9768519/gsl-uniform-random-number-generator
long double  randomUniformFromGsl(){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_ranlux389; // Generator setup
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, mySeed);
    long double  u = gsl_rng_uniform(r); // Generate it!
    gsl_rng_free (r);
    return (float)u;
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

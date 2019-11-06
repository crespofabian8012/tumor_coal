//
//  random.cpp
//  simul
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//

#include <cmath>
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "random.h"
#include <sys/time.h>
/***************************** RandomUniform **********************************/
/* It returns a random uniform variate in range 0..1. It is described in
 Park, S. K. and K. W. Miller. 1988. Random number generators: good
 ones are hard to find. Communications of the ACM, 31(10):1192-1201.
 */

double RandomUniform (long int *seed)
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

double    RandomGamma (double shape, long int *seed)
{
    double gammaNumber = 0;
    
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
double RandomGamma1 (double s, long int *seed)
{
    /* Random standard gamma for s<1
     switching method
     */
    double            r, x=0.0, small=1e-37, w;
    static double   a, p, uf, ss=10.0, d;
    
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
double RandomGamma2 (double s, long int *seed)
{
    /* Random standard gamma for s>1
     Best's (1978) t distribution method
     */
    double            r ,d, f, g, x;
    static double    b, h, ss=0;
    
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

int RandomPoisson (double lambda, long int *seed)
{
    int        poissonNumber;
    double    sum;
    
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

double RandomExponential (double lambda, long int *seed)
{
    double  exponentialNumber, U;
    
    do
        U = randomUniformFromGsl();
    while (U == 0);
    
    exponentialNumber = -log (U) / lambda;
    
    return exponentialNumber;
}

/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

int RandomUniformTo (int max, long int *seed)
{
    double    rd;
    rd = randomUniformFromGsl();
    return (floor(rd*max));
}

/************************************* ChooseUniformState **********************************************/
/* Chooses uniformy a random state according to a vector of state probabilities */
int ChooseUniformState (double *prob, long int *seed)
{
    int            chosenState;
    double        ran, cumProb;
    
    chosenState = 0;
    cumProb = prob[chosenState];
    ran = RandomUniform(seed);
    
    while (ran > cumProb)
        cumProb += prob[++chosenState];
    
    return chosenState;
}


/*************** RandomDirichlet ***************/
//first generates random samples from a gamma and then divide each value by the sum
void  RandomDirichlet (double s, int vectorSize, vector<double> &outputVector, long int *seed)
{   int i;
    double sum=0.0;
    double current;
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
void   randomDirichletFromGsl(int vectorSize, double alpha[], double *theta)
{
    //double theta[vectorSize];
    for (unsigned int i = 0; i < vectorSize; ++i){
        theta[i]=0;
    }
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_ranlux389; // Generator setup
    r = gsl_rng_alloc (T);
    //gsl_rng_set(r, mySeed);
     gsl_ran_dirichlet( r,  vectorSize, alpha, theta); // Generate it!
    gsl_rng_free (r);

}
double RandomLogUniform( double from, double to){
  
    return(exp(from + randomUniformFromGsl()*(to -from)));
}
// from https://stackoverflow.com/questions/9768519/gsl-uniform-random-number-generator
double randomUniformFromGsl(){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_ranlux389; // Generator setup
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, mySeed);
    double u = gsl_rng_uniform(r); // Generate it!
    gsl_rng_free (r);
    return (float)u;
}

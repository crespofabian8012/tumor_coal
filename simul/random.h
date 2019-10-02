//
//  random.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef random_h
#define random_h
static double   RandomUniform (long int *seed);
static double RandomExponential (double mean, long int *seed);
static int    bbinClones (double dat, double *v, int n);
double    RandomGamma (double shape, long int *seed);
double RandomGamma1 (double s, long int *seed);
double RandomGamma2 (double s, long int *seed);
int RandomNegativeBinomial (double mean, double dispersion, long int *seed);
double RandomNormal(double mean, double stddev, long int *seed);
double RandomLogNormal( double mean, double dispersion, long int *seed);
void  RandomDirichlet2 (double *vectorGamma, int vectorSize, double *outputVector, long int *seed);
static int RandomPoisson (double lambda, long int *seed);

int RandomBinomial (double prob, int numTrials, long int *seed);

static int RandomUniformTo (int max, long int *seed);

static double RandomLogUniform( double from, double to, long int *seed);

void  RandomDirichlet (double s, int vectorSize, double **outputVector, long int *seed);
double LogUniformDensity(double value, double from, double to);
#endif /* random_h */

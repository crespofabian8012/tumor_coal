//
//  data_utils.cpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 1/10/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#include "data_utils.hpp"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/timeb.h>
#include <sys/stat.h>
#include <float.h>

#include <zlib.h>
#include <unistd.h>


/*********************************** JC **************************************/
/*    JC performs Jukes-Cantor 69 correction */
/* Note that beta was set such that mean substitution rate will be 1.
 for the JC model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */
void JCmodel (double Pij[4][4], double branchLength, double beta )
{
    int i, j;
    
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if (i == j)
                Pij[i][j] = 0.25 + 0.75*exp(beta*-branchLength);
            else
                Pij[i][j] = 0.25 - 0.25*exp(beta*-branchLength);
        }
    }
}

/*********************************** HKY **************************************/
/*    HKY performs Hasegawa-Kishino-Yano 85 correction */

void HKYmodel (double Pij[4][4], double branchLength, double kappa, double freqR, double freqY, double beta, double freq[4])
{
    int            i, j;
    double        AA1, t, PIj;
    
    t = branchLength;
    
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if (j == A || j == G)    /* purine */
                PIj = freqR;
            else
                PIj = freqY; /* pyrimidine */
            
            AA1 = 1 + PIj*(kappa-1.0);
            
            if (i==j)
                Pij[i][j] = freq[j] + freq[j]*(1/PIj - 1)*exp(-beta*t) + ((PIj-freq[j])/PIj)*exp(-beta*t*AA);
            else if ((i==A && j==G) || (i==C && j==T) || (i==G && j==A) || (i==T && j==C)) /* transition */
                Pij[i][j] = freq[j] + freq[j]*(1/PIj - 1)*exp(-beta*t) - (freq[j]/PIj)*exp(-beta*t*AA);
            else /* transversion */
                Pij[i][j] = freq[j]*(1-exp(-beta*t));
        }
    }
}


/*************** GTR **********************/
void GTRmodel (double Pij[4][4], double branchLength, double Root[], double Cijk[])
{
    int     i, j, k;
    double    t, expt[4];
    
    t = branchLength;
    
    /* P(t)ij = SUM Cijk * exp{Root*t} */
    for (k=1; k<4; k++)
        expt[k]=exp(t*Root[k]);
    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
        {
            Pij[i][j]=Cijk[i*4*4+j*4+0];
            for (k=1; k<4; k++)
                Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
        }
}

/********************* FillSubstitutionMatrix **********************/
/* Sets the apropriate model of nucleotide substitution   */
void FillSubstitutionMatrix (double ch_prob[4][4], double branchLength, int doJC, int doHKY, int doGTR, int doGTnR, double beta, double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[])
{
    int i, j;
    
    if (branchLength<1e-6)
    {
        for (i=0; i<4; i++)
        {
            for (j=0; j<4; j++)
            {
                if (i == j)
                    ch_prob[i][j] = 1.0;
                else
                    ch_prob[i][j] = 0.0;
            }
        }
    }
    else if (doJC == YES)
        JCmodel (ch_prob, branchLength,  beta);
    else if (doHKY == YES)
        HKYmodel (ch_prob, branchLength, kappa,  freqR,  freqY,  beta, freq);
    else if (doGTR == YES)
        GTRmodel (ch_prob, branchLength,  Root,  Cijk);
    else if (doGTnR == YES)
        GTRmodel (ch_prob, branchLength, Root,  Cijk);
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

/************************************* SimulateFiniteDNAforSite **********************************************/
/* Simulates JC, HKY, GTR or GTRnr for a given site */
void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,SiteStr* allSites,  long int *seed, int rateVarAmongSites, double altModelMutationRate, int *numMU, int doJC, int doHKY, int doGTR, int doGTnR, double beta,  double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[])
{
    double    branchLength, Pij[4][4];
    int     cell, anccell, ancstate, newstate;
    
    if (p != NULL)
    {
        // if (p->isHealthyRoot == NO)
        if (p->isOutgroup== NO)
        {
            cell = p->label;
            if (p->anc1 !=NULL)
            {

                anccell = p->anc1->label;
                if (genome == MATERNAL)
                    ancstate = p->anc1->maternalSequence[site];
                else
                    ancstate = p->anc1->paternalSequence[site];
                // }
                // if (doUserTree == YES)
                //branchLength = p->branchLength;
                
                
                if (rateVarAmongSites == YES){
                    // branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
                    branchLength =  p->length * allSites[site].rateMultiplier;
                }
                else{
                    branchLength =  p->length;
                    // branchLength = altModelMutationRate * p->length;
                    
                }
                
                
                FillSubstitutionMatrix (Pij, branchLength,
                                        doJC, doHKY,  doGTR,  doGTnR, beta,  kappa,  freqR,  freqY,  freq,  Root,  Cijk);
                if (genome ==MATERNAL)
                    newstate =p->maternalSequence[site]=ChooseUniformState (Pij[ancstate], seed);
                else// paternal;
                    newstate =p->paternalSequence[site]=ChooseUniformState (Pij[ancstate], seed);
                //newstate = data[genome][cell][site] = ChooseUniformState (Pij[ancstate], seed);
                
                if (newstate != ancstate)
                {
                    if (genome == MATERNAL)
                    {
                        allSites[site].numMutationsMaternal++;
                        
                    }
                    else if (genome == PATERNAL){
                        
                        allSites[site].numMutationsPaternal++;
                    }
                    
                    allSites[site].numMutations++;
                    (*numMU)++;
                }
            }
        }
        SimulateFiniteDNAforSite (p->left,  genome, site,allSites, seed,  rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY, doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk);
        SimulateFiniteDNAforSite (p->right, genome, site,allSites, seed, rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq ,Root,  Cijk);
    }
}

/************************************************************/
/********************* CalculateH ********************/
/* Calculate H for ProbabilityCloneiFromClonej function.
 */

double CalculateH (double t, double TOrigin, double delta)
{
    double  H, AboveTerm, BelowTerm, firstTerm, secondTerm;
    double  a, b;
    
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    a = 0.0;
    b = 0.0;
    H = 0.0;
    
    //printf ("\nInput(H), t=%lf T=%lf delta=%lf ", t, T, delta);
    
    a = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    firstTerm = a * a;
    secondTerm = exp(-1.0 * delta * t);
    AboveTerm = firstTerm * secondTerm;
    //printf ("\nAboveTerm(H) = %lf (%lf %lf %lf) / delta = %lf, T = %lf, t = %lf", AboveTerm, a, firstTerm, secondTerm, delta, T, t);
    
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    BelowTerm = b * b;
    //printf ("\nBelowTerm(H) = %lf", BelowTerm);
    
    H = AboveTerm / BelowTerm;
    //printf ("\nH = %lf", H);
    
    return H;
}

double FmodelTstandard (double t, double TOrigin, double delta)
{
    double  ModelTimeF, firstTerm, secondTerm, thirdTerm;
    double  a, b, c;
    
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    /*fprintf (stderr, "t = %lf\n", t);
     fprintf (stderr, "T = %lf\n", T);
     fprintf (stderr, "delta = %lf\n", delta);*/
    
    a = 1.0 - exp(-1.0 * delta * TOrigin);
    firstTerm = a * a;
    //fprintf (stderr, "firstTerm = %lf\n", firstTerm);
    
    secondTerm = exp(delta * TOrigin);
    //fprintf (stderr, "secondTerm = %lf\n", secondTerm);
    
    b = 1.0 / (1.0 - exp(-1.0 * delta * (TOrigin - t)));
    c = 1.0 / (1.0 - exp(-1.0 * delta * TOrigin));
    thirdTerm = b - c;
    //fprintf (stderr, "thirdTerm = %lf\n", thirdTerm);
    
    
    ModelTimeF = firstTerm * secondTerm * thirdTerm;
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
    
    // New formula from Carsten!
    ModelTimeF = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    
    a = exp(delta * t) - 1.0;
    b = 1.0 - exp(-1.0 * delta * TOrigin);
    c = 1.0 - exp(-1.0 * delta * (TOrigin - t));
    
    ModelTimeF = a * b / (delta * c);
    //fprintf (stderr, "ModelTimeF = %lf\n", ModelTimeF);
    
    
    return ModelTimeF;
}



/************************************************************/
/********************* GmodelTstandard ********************/
/* Conversion from model time t>=0 to standard time.
 */

double GstandardTmodel (double V, double TOrigin, double delta)
{
    double  StandardTimeG, firstTerm, secondTerm, thirdTerm;
    double  a, b, c, d, e;
    
    StandardTimeG = 0.0;
    firstTerm = 0.0;
    secondTerm = 0.0;
    thirdTerm = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    
    //fprintf (stderr, "\nV = %lf, T = %lf, delta = %lf\n", V, T, delta);
    
    
    firstTerm = TOrigin;
    //fprintf (stderr, "\nfirstTerm = %lf\n", firstTerm);
    
    secondTerm = 1 / delta;
    // fprintf (stderr, "secondTerm = %lf\n", secondTerm);
    
    
    //a = 1.0 - exp(-1.0 * delta * T);
    a =  exp(-1.0 * delta * TOrigin);
    
    //fprintf (stderr, "a = %lf\n", a);
    //b = (a * a) * exp(delta * T);
    b = (1 - a) * (1 - a) * (1.0 / a);
    //fprintf (stderr, "b= %lf\n", b);
    
    //c = 1.0 - exp (-1.0 * delta * T);
    c = 1 - a;
    // fprintf (stderr, "c = %lf\n", c);
    //d = c * exp(delta * T);
    d = (1 - a) * (1.0 / a);
    // fprintf (stderr, "d = %lf\n", d);
    e = V + d;
    //fprintf (stderr, "e = %lf\n", e);
    thirdTerm = log(1 - b / e);
    //fprintf (stderr, "valueOfLog = %lf\n", thirdTerm);
    thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
    //fprintf (stderr, "ArgumentOfLog = %lf\n", 1 - b/e);
    //fprintf (stderr, "valueOfLog = %lf\n", thirdTerm);
    
    
    thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
    
    //StandardTimeG = firstTerm + (secondTerm * thirdTerm);
    StandardTimeG = secondTerm * thirdTerm;
    //fprintf (stderr, "StandardTimeG = %lf\n", StandardTimeG);
    
    if ( (1 + delta * V - a) <= 0 ||   (1 + (delta * V - 1)*a ) <= 0 ) // do approximation if required
    {
        fprintf (stderr, "\nApplying approximation of math formula to avoid log(0)\n");
        StandardTimeG = 0.0;
        firstTerm = 0.0;
        secondTerm = 0.0;
        thirdTerm = 0.0;
        a = 0.0;
        b = 0.0;
        c = 0.0;
        d = 0.0;
        e = 0.0;
        
        a = 1 / delta;
        b = log(1 + delta * V);
        firstTerm = a * b;
        //fprintf (stderr, "\nfirstTerm = %lf\n", firstTerm);
        
        
        d = (V * V * delta * exp(-1.0 * delta * TOrigin)) / (1 + V);
        secondTerm =  d;
        //fprintf (stderr, "secondTerm = %lf\n", secondTerm);
        
        StandardTimeG = firstTerm - secondTerm;
        //fprintf (stderr, "StandardTimeG = %lf\n", StandardTimeG);
    }
    
    return StandardTimeG;
}


/***************************** compare******************************/

int compare (const void * a, const void * b)

{
    
    double *p1 = (double *)a;
    
    double *p2 = (double *)b;
    
    
    
    if (*p1 > *p2) return 1;
    
    else if (*p2 > *p1 )return -1;
    
    else return 0;
    
}

/********************************** SimulateTriNucFreqGenome ***********************************/

/*
 
 Simulate a homozygous diploid genome with the trinucleotide frequencies of the human genome
 
 
 
 We simulate first a random trinucleotide.
 
 Then we go site by site taking into account the last two letters  of the
 
 previous trinucleotide.  if the first trinucleotide is, for example, CAT
 
 then we simulate trinucleotide starting in AT_  (i.e. ATA, ATC, ATG or ATT),
 
 and keep going  taking into account the last two letters of the new trinucleotide,
 
 and so on.
 
 */



void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, SiteStr* allSites, int doGeneticSignatures, double cumfreq[4], double *triNucFreq )

{
    
    int         chosenTriNucleotide, nextNucleotide;
    
    int            k, n1, n2, n3, rest, site;
    
    double         *prob4, sum;
    
    
    
    /* memory allocations */
    
    prob4 = (double *) calloc (4, sizeof(double));
    
    if (!prob4)
        
    {
        
        fprintf (stderr, "Could not allocate the prob4 vector\n");
        
        exit (-1);
        
    }
    
    
    
    /* choose first trinucleotide */
    
    chosenTriNucleotide = ChooseUniformState(triNucFreq, seed);
    
    
    
    /* find bases of the selected trinucleotide */
    
    n1 = chosenTriNucleotide/16;
    
    rest = chosenTriNucleotide%16;
    
    n2 = rest/4;
    
    n3 = rest%4;
    
    //fprintf (stderr, "\n%2d %c%c%c ", chosenTriNucleotide, WhichNuc(n1), WhichNuc(n2), WhichNuc(n3));
    
    
    
    site = 0;
    
    p->maternalSequence[site]=p->paternalSequence[site]=n1;
    
    //    data[MATERNAL][cell][site] = data[PATERNAL][cell][site] = allSites[site].referenceAllele = n1;
    
    site++;
    
    p->maternalSequence[site]=p->paternalSequence[site] = allSites[site].referenceAllele = n2;
    
    site++;
    
    p->maternalSequence[site]=p->paternalSequence[site]= allSites[site].referenceAllele = n3;
    
    
    
    /* fill the rest of the genome */
    
    /* choose next nucleotide given the last two bases of the previous trinucleotide  */
    
    for (site=3; site<numSites; site++)
        
    {
        
        /* normalize frequencies given the last two bases */
        
        sum = 0;
        
        for (k=0; k<4; k++)
            
            sum += triNucFreq[trinuc(n2,n3,k)];
        
        for (k=0; k<4; k++)
            
            prob4[k] = triNucFreq[trinuc(n2,n3,k)] / sum;
        
        
        
        nextNucleotide = ChooseUniformState(prob4, seed);
        
        p->maternalSequence[site]=p->paternalSequence[site]= allSites[site].referenceAllele = nextNucleotide;
        
        
        
        /* move downstream one position */
        
        n1 = n2;
        
        n2 = n3;
        
        n3 = nextNucleotide;
        
    }
    
    
    
    free(prob4);
    
    prob4=NULL;
    
}

/***************************** ReadParametersFromFile *******************************/
/* Reads parameter values from the parameter file */

void ReadParametersFromFile(ProgramOptions *programOptions, FilePaths *filePaths,
                            int **CloneNameBegin,
                            int **CloneSampleSizeBegin,
                            int **ClonePopSizeBegin,
                            double **CloneBirthRateBegin,
                            double **CloneDeathRateBegin,
                            double **CloneTimeOriginInput,
                            double Mij[4][4],
                            double freq[4]
                            )
{
    int   j, z;
    char  ch;
    float   argument;
    double    sumPi, sum;
    double argumentDouble;
    double argumentDouble1;
    int argumentInt;
    long int argumentLongInt;
    int *pInt;
    double  *pDouble;
    long int *pLongInt;
    
    /* Used: N X C R D M O T K Y # */
    
    if (feof(stdin))
    {
        fprintf(stderr, "PARAMETER ERROR: Unable to read parameters from stdin\n");
        exit(0);
    }
    
    ch = fgetc(stdin);
    while (isspace(ch))
        ch = fgetc(stdin);
    while (ch == '[')
    {
        ReadUntil(stdin, ']', "closing bracket");
        ch = fgetc(stdin);
        while (isspace(ch))
            ch = fgetc(stdin);
    }
    
    while (!feof(stdin))
    {
        argument = 0;
        ch = toupper(ch);
        switch (ch)
        {
                
            case 'N':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", (int)*pInt);
                    PrintUsage();
                }
                programOptions->numDataSets =argumentInt;
                break;
            case '#':
                if (fscanf(stdin, "%lu bytes", &argumentLongInt) != 1 || argumentLongInt < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad seed (#) (%d)\n\n", (int)argumentLongInt);
                    PrintUsage();
                }
                programOptions->userSeed =argumentLongInt;
                break;
                
                
            case 'X':
                if (fscanf(stdin, "%f", &argument) != 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad haplid/diploid chosen (x) (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*Nscaling = (int) argument;
                programOptions->Nscaling =(int) argument;
                if (programOptions->Nscaling < 1 || programOptions->Nscaling > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Haploid/diplod option (x) (1-2) (%d)\n\n", programOptions->Nscaling);
                    PrintUsage();
                }
                break;
            case 'C':
                if (fscanf(stdin, "%f", &argument) != 1 || argument <= 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of clones (must be 1 or higher) (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*numClones = (int) argument;
                programOptions->numClones=(int) argument;
                
                *CloneNameBegin =  (int*)malloc( programOptions->numClones * (long) sizeof(int));
                *CloneSampleSizeBegin = (int*)malloc( programOptions->numClones* (long) sizeof(int));
                *ClonePopSizeBegin = (int*)calloc( programOptions->numClones , (long) sizeof(int));
                *CloneBirthRateBegin =  (double*)calloc( programOptions->numClones, (long) sizeof(double));
                *CloneDeathRateBegin =  (double*)calloc( programOptions->numClones, (long) sizeof(double));
                *CloneTimeOriginInput =  (double*)calloc( programOptions->numClones, (long) sizeof(double));
                
                if (*CloneNameBegin == NULL || *CloneSampleSizeBegin == NULL || *ClonePopSizeBegin == NULL || *CloneBirthRateBegin == NULL || *CloneDeathRateBegin == NULL || *CloneTimeOriginInput == NULL)
                {
                    fprintf (stderr, "PARAMETER ERROR: Could not allocate variables for clones\n");
                    exit (1);
                }
                
                for (j = 0; j <  programOptions->numClones; j++)
                {
                    
                    for (z = 1; z <= 6; z++)
                    {
                        if (z == 1)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneNameBegin + j) = (int) argument;
                            
                            
                            if ( *(*CloneNameBegin + j)  <= 0 ||  *(*CloneNameBegin + j)  >  programOptions->numClones)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad number for clone %d (should be higher than 0 and lower than the number of clones %d) (%d)\n\n", j,  programOptions->numClones, *(*CloneNameBegin + j) );
                                PrintUsage();
                            }
                        }
                        if (z == 2)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneSampleSizeBegin + j) = (int) argument;
                            programOptions->numCells= programOptions->numCells + (int) argument;
                            if (*(*CloneSampleSizeBegin + j)  < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad sample size for clone %d (should not be negative) (%d)\n\n", j, *(*CloneSampleSizeBegin + j) );
                                PrintUsage();
                            }
                        }
                        if (z == 3)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*ClonePopSizeBegin + j)  = (int) argument;
                            if (*(*ClonePopSizeBegin + j) < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad population size for clone %d (should be higher than 0) (%d)\n\n", j, *(*ClonePopSizeBegin + j) );
                                PrintUsage();
                            }
                        }
                        if (z == 4)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneBirthRateBegin + j) = (double) argument;
                        }
                        if (z == 5)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneDeathRateBegin + j) = (double) argument;
                        }
                        if (z == 6)
                        {
                            fscanf(stdin, "%f", &argument);
                            *(*CloneTimeOriginInput + j) = (double) argument;
                            if (*(*CloneTimeOriginInput + j)  < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be negative) (%lf)\n\n", j,  *(*CloneTimeOriginInput + j) );
                                PrintUsage();
                            }
                        }
                        
                    }
                    
                }
                break;
                
                
            case 'U':
                
                if (fscanf(stdin, "%lf", &argumentDouble) != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mutation rate (%f) \n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->mutationRate=  (double) argumentDouble;
                break;
                
            case 'B':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0  || argument > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphabet (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*alphabet = (int) argument;
                programOptions->alphabet =(int) argument;
                break;
                
            case 'Y':
                if (fscanf(stdin, "%d",  &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", (int) argumentInt);
                    PrintUsage();
                }
                programOptions->noisy =(int) argumentInt;
                break;
                
            case 'D':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 ||argumentDouble < 0 ||argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic dropout rate (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->ADOrate =(double) argumentDouble;
                break;
            case 'O':
                if (fscanf(stdin, "%d", &argumentInt) < 0 )
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    PrintUsage();
                }
                programOptions->outgroupSelection = argumentInt;
                if ( programOptions->outgroupSelection != 0 &&  programOptions->outgroupSelection != 1 &&  programOptions->outgroupSelection != 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n",  programOptions->outgroupSelection);
                    PrintUsage();
                }
                
                if ( programOptions->outgroupSelection == 0)
                {
                    programOptions->thereisOutgroup = NO;
                    programOptions->outgroupBranchLength_RootSample = 0.0;
                    programOptions->outgroupBranchLength_Root1Root2 = 0.0;
                }
                else if (programOptions->outgroupSelection == 1)
                {
                    //*thereisOutgroup = YES;
                    programOptions->thereisOutgroup=YES;
                    programOptions->outgroupBranchLength_Root1Root2 = 0.0;
                    
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root-Sample) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions->outgroupBranchLength_RootSample=argumentDouble;
                }
                else if (programOptions->outgroupSelection == 2)
                {
                    //*thereisOutgroup = YES;
                    programOptions->thereisOutgroup=YES;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root1-Root2) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions->outgroupBranchLength_Root1Root2=argumentDouble;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root2-Sample) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions->outgroupBranchLength_RootSample=argumentDouble;
                }
                else
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    PrintUsage();
                }
                //outgroupBranchLength_RootSample = 0;
                break;
            case 'I':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic imbalance (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->allelicImbalance=argumentDouble;
                break;
            case 'R':
                if (programOptions->doHKY == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
                    PrintUsage();
                }
                if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf",
                           &Mij[0][0], &Mij[0][1], &Mij[0][2], &Mij[0][3],
                           &Mij[1][0], &Mij[1][1], &Mij[1][2], &Mij[1][3],
                           &Mij[2][0], &Mij[2][1], &Mij[2][2], &Mij[2][3],
                           &Mij[3][0], &Mij[3][1], &Mij[3][2], &Mij[3][3])!=16)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix (-rx x x x x x x x x x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
                    PrintUsage();
                }
                
                if (Mij[0][0] != 0  || Mij[1][1] != 0 || Mij[2][2] != 0 || Mij[3][3] != 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix: diagonals should be 0 \n\n");
                    PrintUsage();
                }
                //*thereIsMij = YES;
                programOptions->thereIsMij=YES;
                if (CheckMatrixSymmetry (Mij) == YES)
                {
                    
                    //programOptions->doJC = NO;
                    programOptions->doJC = YES;
                    programOptions->doHKY = NO;
                    //programOptions->doGTR = YES;
                    programOptions->doGTR = NO;
                    programOptions->doGTnR = NO;
                    
                }
                else
                {
                    
                    programOptions->doJC = NO;
                    programOptions->doHKY = NO;
                    //  programOptions->doGTR = NO;
                    programOptions->doGTR =  YES;
                    // programOptions->doGTnR = YES;
                    programOptions->doGTnR = NO;
                }
                break;
            case 'P':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->propAltModelSites=argumentDouble;
                if (programOptions->propAltModelSites < 0 || programOptions->propAltModelSites > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f). It has to be between 0 and 1\n\n", programOptions->propAltModelSites);
                    PrintUsage();
                }
                if (programOptions->propAltModelSites > 0 && programOptions->altModel != ISMhap && programOptions->doSimulateFixedNumMutations == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a proportion of non-ISM  sites bigger than zero if the number of mutations is fixed\n\n");
                    PrintUsage();
                }
                if (programOptions->alphabet == DNA && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) programOptions->altModel);
                        PrintUsage();
                    }
                }
                else if (programOptions->alphabet == BINARY && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) programOptions->altModel);
                        PrintUsage();
                    }
                }
                break;
            case 'M':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative mutation model (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->altModel = (int) argument;
                if (programOptions->alphabet == DNA && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
                        PrintUsage();
                    }
                }
                else if (programOptions->alphabet == BINARY && programOptions->propAltModelSites > 0)
                {
                    if (programOptions->altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the model (%d) specified are incompatible", (int) argument);
                        PrintUsage();
                    }
                }
                break;
            case 'T':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths->treeFile, "trees");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths->treeFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths->treeFile[j] = '\0';
                }
                //*doPrintTrees = YES;
                programOptions->doPrintTrees = YES;
                break;
            case 'A':
                if (fscanf(stdin, "%lf %lf %d", &argumentDouble, &argumentDouble1, &argumentInt) != 3)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean/var/model amplification error (%f ; %f ; model=%d)\n\n",argumentDouble, argumentDouble1, argumentInt);
                    PrintUsage();
                }
                programOptions->meanAmplificationError= argumentDouble;
                programOptions->varAmplificationError= argumentDouble1;
                programOptions->simulateOnlyTwoTemplates=argumentInt;
                
                if ( programOptions->meanAmplificationError < 0 ||  programOptions->meanAmplificationError > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean amplification error (%f)\n\n",  programOptions->meanAmplificationError);
                    PrintUsage();
                }
                if ( programOptions->varAmplificationError < 0 || ( programOptions->meanAmplificationError > 0 &&  programOptions->varAmplificationError >= ( programOptions->meanAmplificationError * (1.0 -  programOptions->meanAmplificationError))))
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n",  programOptions->meanAmplificationError);
                    PrintUsage();
                }
                if ( programOptions->simulateOnlyTwoTemplates != 0 &&  programOptions->simulateOnlyTwoTemplates != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 4 templates) or 1 (assume 2 templates)",  programOptions->simulateOnlyTwoTemplates);
                    PrintUsage();
                }
                break;
            case 'E':
                if (fscanf(stdin, "%lf",  &argumentDouble) !=1 ||  argumentDouble < 0 ||  argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing error (%f)\n\n",  argumentDouble);
                    PrintUsage();
                }
                programOptions->sequencingError= argumentDouble;
                break;
            case 'S':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative  do simulated data or not (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->doSimulateData= (int)argument;
                break;
            case 'K':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths->timesFile, "times");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths->timesFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths->timesFile[j] = '\0';
                }
                //*doPrintTimes = YES;
                programOptions->doPrintTimes = YES;
                break;
            case 'J':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of mutations (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->numFixedMutations = (int) argument;
                //programOptions->doSimulateFixedNumMutations = YES;
                programOptions->doSimulateFixedNumMutations = YES;
                if (programOptions->propAltModelSites > 0 && programOptions->altModel != ISMhap)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a fixed number of mutations if there is any non-ISM  site. Set the proportion of non-ISM diploid sites to zero\n\n");
                    PrintUsage();
                }
                break;
            case 'F':
                if (fscanf(stdin, "%lf %lf %lf %lf", &freq[0], &freq[1], &freq[2], &freq[3])!=4)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad Base Frequencies\n\n");
                    PrintUsage();
                }
                else if (freq[0] == freq[1] == freq[2] == freq[3])
                    programOptions->equalBaseFreq = YES;
                else
                    programOptions->equalBaseFreq = NO;
                sumPi = freq[0] + freq[1] + freq[2] + freq[3];
                if (sumPi != 1.0)
                {
                    freq[0]/=sumPi;
                    freq[1]/=sumPi;
                    freq[2]/=sumPi;
                    freq[3]/=sumPi;
                }
                break;
            case '?':
                PrintUsage();
                break;
                /*case 'H':
                 PrintUsage();
                 break;*/
            case 'V':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble <= 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad coverage dispersion (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->rateVarCoverage = YES;
                break;
            case 'Z':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions->SNPrate = YES;
                break;
                
            case 'H':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing coverage (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions->coverage = (int) argument;
                if (programOptions->coverage > 0){
                    programOptions->doSimulateReadCounts=YES;
                    // *doSimulateReadCounts = YES;
                    
                }
                if (programOptions->genotypingError > 0 && programOptions->doSimulateReadCounts == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
                    PrintUsage();
                }
                break;
            default :
                fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
                PrintUsage();
                break;
        }
        ch = fgetc(stdin);
        while (isspace(ch) && !feof(stdin))
            ch = fgetc(stdin);
        while (ch == '[')
        {
            ReadUntil(stdin, ']', "closing bracket");
            ch = fgetc(stdin);
            while (isspace(ch))
                ch = fgetc(stdin);
        }
    }
    
}
/***************************** PrintUsage *******************************/
/* Prints a short description of program usage */
void PrintUsage()
{
    fprintf (stderr, "\n\nUsage: %s%s [-n# -x# -c# (# # # # # #) -u# -o# -ttrees.tre -ktimes.txt -## -y# -? -h]", PROGRAM_NAME, VERSION_NUMBER);
    fprintf (stderr, "\n-n: number of replicates (e.g. -n1000)");
    fprintf (stderr, "\n-x: haploid/diploid (1 for haploid, 2 for diploid) (e.g. -x2)");
    fprintf (stderr, "\n-c: number of clones; for each clone: ID number, sample size, population size, birth rate, death rate, time to origin (e.g. -c5\n\t1 4 70000 0.3 0.1 56\n\t2 5 60000 0.4 0.3 110\n\t3 5 95000 0.4 0.3 117\n\t4 4 15000 0.5 0.4 95\n\t5 3 44000 0.3 0.1 53)");
    fprintf (stderr, "\n-u: mutation rate (e.g. -u9.1e-6)");
    fprintf (stderr, "\n-o: branch length to the outgroup (root-outgroup) (e.g. -o0.0325) (default is no outgroup)");
    
    /* outputs */
    fprintf (stderr, "\n-t: tree file name (e.g. -ttrees.tre)");
    fprintf (stderr, "\n-k: times file name (e.g. -ktimes.txt)");
    
    /* other settings */
    fprintf (stderr, "\n-#: seed (e.g. -#37864287)");
    fprintf (stderr, "\n-y: noisy (amount of information printed on the screen) (e.g. -y2)");
    fprintf (stderr, "\n-? -h: Print help\n");
    
    exit(-1);
}

/***************************** ReadUntil *******************************/
/* Reading in between [] from input files */
void ReadUntil(FILE *fv, char stopChar, char *what)
{
    char ch;
    
    ch = fgetc(fv);
    while (!feof(fv) && ch != stopChar)
        ch = fgetc(fv);
    
    if (feof(fv) || ch != stopChar)
    {
        fprintf(stderr, "%s missing", what);
        exit(0);
    }
}

/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int CheckMatrixSymmetry(double matrix[4][4])
{
    int i,j;
    
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            if(matrix[i][j] != matrix[j][i])
                return NO;
    return YES;
}

/***************************** InitListClonesTimes*******************************/
/* InitListClonesTimes*/
void InitListClonesTimes(Population **populations, int numClones,  int *doEstimateTimesOriginClones,
                         double *CloneTimeOriginInput
                         ) {
    int z, j;
    z = 0;
    Population *p;
    for (j = 0; j <= (numClones - 1); j++)
    {
        p = *(populations + j);
        p->timeOriginInput = CloneTimeOriginInput[j ];
        p->timeOriginSTD = CloneTimeOriginInput[j ] / p->effectPopSize;
        p->timeMigrationSTDCurrentPop = p->timeOriginSTD;
        p->delta = (double)p->growthRate * p->effectPopSize;
        if ( p->timeOriginSTD == 0) {
            z++;
            p->doEstimateTimeOrigin = 1;
        }
    }
    // active estimation of times
    if (z == numClones)
        *doEstimateTimesOriginClones = YES;
    for (j = 0; j < numClones; j++)
    {
        p = *(populations + j);
        if ( p->delta <= 0)
        {
            fprintf (stderr, "PARAMETER ERROR: The growth rate cannot be lower than the death rate(Delta parameter negative, Delta=(%10.9lf)) for population %d\n\n", p->delta, j);
            PrintUsage();
        }
        if (p->timeOriginInput == 0 && *doEstimateTimesOriginClones == NO)
        {
            fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be 0; excepting estimation of times where all clones must have a time of 0) (%lf)\n\n", j, p->timeOriginInput);
            PrintUsage();
        }
    }
}
/***************************** InitListClones*******************************/
/* InitListClones*/
void InitListClones(Population **populations, int numClones, int noisy, int *CloneNameBegin, int *CloneSampleSizeBegin, double *CloneBirthRateBegin,  double *CloneDeathRateBegin,  int *ClonePopSizeBegin, int TotalNumSequences  ) {
    int z;
    //struct Population* pops = malloc(numClones * (sizeof(struct Population)+TotalNumSequences * sizeof( int) + numClones * sizeof(double) ));
    struct Population* pops =(Population*) malloc(numClones * (sizeof(struct Population) ));
    for (z = 0; z <= (numClones - 1); z++)
    {
        //CloneGrowthRateBegin[z ] = 0.0;
        // ClonePopSizeMeffectBegin[z ] = 0.0;
        // CloneGrowthRateBegin[z ] = CloneBirthRateBegin[z ] - CloneDeathRateBegin[z ];
        // ClonePopSizeMeffectBegin[z] = ClonePopSizeBegin[z ] / CloneBirthRateBegin[z ];
        pops[z].index = CloneNameBegin[z ];
        pops[z].order=0;
        pops[z].birthRate = CloneBirthRateBegin[z ];
        pops[z].deathRate = CloneDeathRateBegin[z ];
        pops[z].growthRate = CloneBirthRateBegin[z ] - CloneDeathRateBegin[z ];
        pops[z].sampleSize = CloneSampleSizeBegin[z];
        pops[z].popSize = ClonePopSizeBegin[z ];
        pops[z].effectPopSize =  ClonePopSizeBegin[z ] / CloneBirthRateBegin[z ];
        pops[z].delta = pops[z].growthRate * pops[z].effectPopSize ;
        pops[z].numActiveGametes = CloneSampleSizeBegin[z];
        pops[z].isAlive = 1;
        pops[z].nodeIdAncesterMRCA = 0;
        pops[z].numCompletedCoalescences = 0;
        pops[z].nextAvailableIdInmigrant = 0;
        pops[z].numIncomingMigrations = 0;
        pops[z].numPossibleMigrations = 0;
        pops[z].doEstimateTimeOrigin = NO;
        *(populations + z) = &pops[z];
        if (noisy > 1)
            printf("\t%d\t\t", CloneNameBegin[z ]);
        if ((z + 1) != CloneNameBegin[z ])
        {
            fprintf (stderr, "PARAMETER ERROR: Check order of clones. Clone (%d) in order is different to (%d). (d)\n\n", z, CloneNameBegin[z]);
            PrintUsage();
        }
        if (noisy > 1)
        {
            printf("\t\t\t%d\t\t\t", CloneSampleSizeBegin[z ]);
            printf("\t%d\t\t\t", ClonePopSizeBegin[z ]);
            printf("\t%lf\t\t",  pops[z].effectPopSize);
            printf("\t\t\t%lf\t", CloneBirthRateBegin[z ]);
            printf("\t%lf\t", CloneDeathRateBegin[z ]);
            printf("\t%lf\t\t",pops[z].growthRate);
            //printf("\t%d\t\t", z);
            printf("\n");
        }
    }
}
 void InitNumberNodes(double *TotalBirthRate, double *TotalDeathRate, int *TotalN,  Population **populations, ProgramOptions *programOptions) {
    programOptions->TotalNumSequences = 0;
    *TotalN = 0;
    *TotalBirthRate = 0.0;
    *TotalDeathRate = 0.0;
    Population* popI;
    int j;
    for (j = 0; j < programOptions->numClones; j++)
    {   popI = *(populations + j);
        popI->FatherPop =NULL;
        programOptions->TotalNumSequences = programOptions->TotalNumSequences + popI->sampleSize;
        *TotalN = *TotalN + popI->popSize;
        *TotalBirthRate = *TotalBirthRate + popI->birthRate;
        *TotalDeathRate = *TotalDeathRate + popI->deathRate;
    }
    programOptions->numNodes = 2 * programOptions->TotalNumSequences + programOptions->numClones+ 10;
    
    programOptions->numCells =programOptions->TotalNumSequences;
}
/***************************** ListClonesAccordingTimeToOrigin*******************************/
/* ListClonesAccordingTimeToOrigin*/
void ListClonesAccordingTimeToOrigin(Population **populations, int numClones) {
    
    
    qsort(populations, numClones, sizeof(Population* ), comparePopulationsByTimeOrigin);
}
/***************************** comparePopulationsByTimeOrigin*******************************/
int comparePopulationsByTimeOrigin(const void *s1, const void *s2)
{
    struct Population *p1 = *(struct Population **)s1;
    struct Population *p2 = *(struct Population **)s2;
    if (  p1->timeOriginInput  > p2 ->timeOriginInput)
        return 1;
    else if (p1->timeOriginInput < p2->timeOriginInput)
        return -1;
    else
        return 0;
}
void InitFilesPathsOptions( FilePaths *filePaths, ProgramOptions *programOptions)
{
    strcpy(filePaths->resultsDir, "Results");
    strcpy(filePaths->treeDir, "trees_dir");
    strcpy(filePaths->timesDir, "times_dir");
    strcpy(filePaths->SNVgenotypesDir, "snv_genotypes_dir");
    strcpy(filePaths->SNVhaplotypesDir, "snv_haplotypes_dir");
    strcpy(filePaths->trueHaplotypesDir, "true_haplotypes_dir");
    strcpy(filePaths->fullHaplotypesDir, "full_haplotypes_dir");
    strcpy(filePaths->MLhaplotypesDir, "ML_haplotypes_dir");
    strcpy(filePaths->fullGenotypesDir, "full_genotypes_dir");
    strcpy(filePaths->CATGdir, "catg_dir");
    strcpy(filePaths->VCFdir, "vcf_dir");
    
    strcpy(filePaths->SNVgenotypesFile, "snv_gen");
    strcpy(filePaths->SNVhaplotypesFile, "snv_hap");
    strcpy(filePaths->trueHaplotypesFile, "true_hap");
    strcpy(filePaths->fullHaplotypesFile, "full_hap");
    strcpy(filePaths->MLhaplotypesFile, "ML_hap");
    strcpy(filePaths->fullGenotypesFile, "full_gen");
    strcpy(filePaths->treeFile, "trees");
    strcpy(filePaths->timesFile, "times");
    if (strlen(filePaths->userTreeFile) == 0)
        strcpy(filePaths->userTreeFile, "usertree");
    if (strlen(filePaths->userGenomeFile) == 0)
        strcpy(filePaths->userGenomeFile, "usergenome");
    strcpy(filePaths->CATGfile, "catg");
    strcpy(filePaths->VCFfile, "vcf");
    strcpy(filePaths->logFile, "log");
    strcpy(filePaths->settingsFile, "log");
#ifdef MYDEBUG
    strcpy(filePaths.mutationsFile, "mutations");
#endif
    programOptions->doPrintSNVgenotypes=1;
    programOptions->doPrintSNVhaplotypes=1;
    programOptions->doPrintTrueHaplotypes=1;
    programOptions->doPrintFullGenotypes=1;
    programOptions->doPrintFullHaplotypes=1;
    programOptions->doNGS=1;
    programOptions->doPrintCATG=1;
    if (programOptions->doSimulateData == NO)
    {
        programOptions->doPrintSNVgenotypes = NO;
        programOptions->doPrintSNVhaplotypes = NO;
        programOptions->doPrintTrueHaplotypes = NO;
        programOptions->doPrintFullHaplotypes = NO;
        programOptions->doPrintMLhaplotypes = NO;
        programOptions->doPrintFullGenotypes = NO;
        programOptions->doPrintAncestors = NO;
        programOptions->doSimulateReadCounts = NO;
        programOptions->doPrintCATG = NO;
    }
}
int SimulateData(ProgramOptions *programOptions, int *CloneNameBegin, int *CloneSampleSizeBegin, int *ClonePopSizeBegin,
                 Population **populations,
                 FilePaths *filePaths,
                 Files*files,
                 char *ObservedCellNames[]
                 
                 )
{
    int i,j,k,z;
    TreeNode    *nodes;
     char *newickString2;
    double totalTreeLength;
     int    HEALTHY_ROOT, TUMOR_ROOT;
    int    numISMdeletions, numISMCNLOH;
    int cumNumMUperTree;
    int     numAltModelSites, numDefaultModelSites, numISMmutations, altModel;
    double freq[4];
    double Mij[4][4];
    double cumfreq[4];
    double cumMij[4][4];
    double Eij[4][4];
    double cumEij[4][4];
    char **cellNames;
    double    kappa, beta, freqR, freqY, freqAG, freqCT;
    double    Rmat[6], NRmat[12], Cijk[256], Root[4];
    double *triNucFreq;
    double   cumNumSNVs, cumNumMU, cumNumDEL, cumNumCNLOH, cumCountMLgenotypeErrors;
    double cumNumMUSq, cumNumSNVsSq, cumNumDELSq, cumNumCNLOHSq;
    CellStr     *cell;
    
    if ( programOptions->alphabet == DNA)
    {
        /* initialize cumfreq */
        cumfreq[0] = freq[0];
        for (i=1; i<4; i++)
            cumfreq[i] = cumfreq[i-1] + freq[i];
        
        /* initialize cumMij */
        for (i=0; i<4; i++)
        {
            cumMij[i][0] = Mij[i][0];
            for (j=1; j<4; j++)
                cumMij[i][j] = cumMij[i][j-1] + Mij[i][j];
        }
        
        /* initialize cumEij */
        for (i=0; i<4; i++)
        {
            cumEij[i][0] = Eij[i][0];
            for (j=1; j<4; j++)
                cumEij[i][j] = cumEij[i][j-1] + Eij[i][j];
        }
    }
    
    /* allocate memory for site information (equal for maternal and paternal) */
    SiteStr* allSites = (SiteStr*) calloc (programOptions->numSites, sizeof(SiteStr));
    if (!allSites)
    {
        fprintf (stderr, "Could not allocate the allSites structure\n");
        exit (-1);
    }
    for (i=0; i< programOptions->numSites; i++)
    {
        allSites[i].alternateAlleles = (int *) calloc (4, sizeof(int));
        if (!allSites[i].alternateAlleles)
        {
            fprintf (stderr, "Could not allocate the allSites[i].alternateAlleles[] structure\n");
            exit (-1);
        }
    }
    
    /* the arrays below keep the index for different types of sites */
    int* SNVsites = (int*) malloc (programOptions->numSites* sizeof(int));
    if (!SNVsites)
    {
        fprintf (stderr, "Could not allocate the SNVsites structure\n");
        exit (-1);
    }
   int* SFS = (int*) malloc (programOptions->numSites* sizeof(int));
    if (!SFS)
    {
        fprintf (stderr, "Could not allocate the SNVsites structure\n");
        exit (-1);
    }
    /* the arrays below keep the index for different types of sites */
    int* variantSites = (int*) malloc (programOptions->numSites* sizeof(int));
    if (!variantSites)
    {
        fprintf (stderr, "Could not allocate the variantSites structure\n");
        exit (-1);
    }
    
    int* DefaultModelSites = (int*) malloc (programOptions->numSites* sizeof(int));
    if (!DefaultModelSites)
    {
        fprintf (stderr, "Could not allocate the DefaultModelSites structure\n");
        exit (-1);
    }
    
    int* AltModelSites = (int*) malloc (programOptions->numSites* sizeof(int));
    if (!AltModelSites)
    {
        fprintf (stderr, "Could not allocate the AltModelSites structure\n");
        exit (-1);
    }
    TreeNode  **treeTips;
    treeTips = (TreeNode **) malloc (programOptions->TotalNumSequences * sizeof(TreeNode*));
    if (!treeTips)
    {
        fprintf (stderr, "Could not allocate the treeTips array\n");
        exit (-1);
    }
    double    *proportionsVector = (double *) calloc((programOptions->numClones), (long) sizeof(double));
    if (!proportionsVector)
    {
        fprintf (stderr, "Could not allocate proportions vector (%lu bytes)\n", (programOptions->numClones ) * (long) sizeof(double));
        exit (-1);
    }
    /* Variance memories */
    int * varEvent = (int *) malloc(programOptions->numDataSets * (long) sizeof(int));
    if (!varEvent)
    {
        fprintf (stderr, "Could not allocate varEvent (%lu bytes)\n", programOptions->numDataSets * (long) sizeof(int));
        exit (1);
    }
    double * varTimeGMRCA = (double *) malloc(programOptions->numDataSets* (long) sizeof(double));
    if (!varTimeGMRCA)
    {
        fprintf (stderr, "Could not allocate varTimeGMRCA (%lu bytes)\n", programOptions->numDataSets  * (long) sizeof (double));
        exit (1);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    TreeNode *root;
    
    double      cumNumCA, meanNumCA, cumNumMIG, meanNumMIG, numEventsTot, countTMRCA,TMRCA;
    int        numCA, numMIG;
    int dataSetNum;
    long int seedFirst =  programOptions->seed;
    int  numMU, numDEL, numCNLOH, numProposedMU, numSNVs, numFixedMutations, numSNVmaternal;
    programOptions->doUseObservedCellNames=NO;
    int ***data;
    ValidateParameters(programOptions,CloneNameBegin , CloneSampleSizeBegin, ClonePopSizeBegin);

    InitListPossibleMigrations(populations,programOptions->numClones);
    InitPopulationsCoalescentEvents( programOptions->numClones,  populations);
   
    for (dataSetNum = 0; dataSetNum < programOptions->numDataSets; dataSetNum++)// dataSetNum refers to a simulated tree number
    {
        if (programOptions->doPrintSeparateReplicates == YES)
            PrepareSeparateFiles(0, 1, dataSetNum,  filePaths, programOptions, files);
        /* Reorganize seed per replicate */
        //seed = seedFirst+dataSetNum+10;
        numCA = numMU = numDEL = numProposedMU = TMRCA = numSNVmaternal = 0;
        programOptions->seed = seedFirst + dataSetNum * 1000;
        //fprintf(stderr, "\n seed = %lu \n", seed);
        /* reset variables per replicate */
        if (programOptions->noisy > 0)
        {
            fprintf (stderr, "\rReplicate #%3d/%d", dataSetNum + 1, programOptions->numDataSets);
            fflush (stdout);
        }
        varEvent[dataSetNum] = 0;
        varTimeGMRCA[dataSetNum] = 0.0;
   
        numCA = numMIG = 0;
        countTMRCA = 0.0;
        /* coalescent tree */
        
        if (programOptions->noisy > 1)
            fprintf (stderr, "\n>> Start coalescent tree .. \n");
        MakeCoalescenceTree2 (&(programOptions->seed), populations,
                              &(programOptions->numNodes),
                              programOptions->numClones,
                              programOptions,
                              cumNumCA,
                              meanNumCA,
                              cumNumMIG,
                              meanNumMIG,
                              &numMIG,
                              &numCA,
                              &numEventsTot,
                              &nodes,
                              treeTips,
                              &root,
                              ObservedCellNames, NULL
                              ) ;
        if (programOptions->noisy > 1)
            fprintf (stderr, "\n>> Finishing coalescent tree ... DONE");
        cumNumCA += numCA;
        cumNumMIG += numMIG;
        // countTMRCA = treeRootInit[0]->timePUnits;
        countTMRCA = root->timePUnits;
        
        //fprintf ( stderr, "\n countTMRCA = %lf\n", countTMRCA);
        varTimeGMRCA[dataSetNum] = countTMRCA;
        varEvent[dataSetNum] = numCA + numMIG;
        // counterTime = counterTime + counterTimeInit;
       
        
        /*************** output files *************/
        
      
        if (programOptions->doPrintTrees == YES)
        {
   
            PrintTrees(dataSetNum, &root, files->fpTrees, programOptions->mutationRate, programOptions->doUseObservedCellNames);
            PrintTrees2(dataSetNum, &root, files->fpTrees2, programOptions->mutationRate,  ObservedCellNames, programOptions->doUseObservedCellNames);
        }
        if (programOptions->doPrintTimes == YES)
        {
            PrintTimes(dataSetNum, files->fpTimes, programOptions->mutationRate, nodes, programOptions->thereisOutgroup);
            PrintTimes2(dataSetNum, files->fpTimes2, programOptions->mutationRate, nodes, programOptions->thereisOutgroup);
        }
        if (programOptions->noisy > 1)
        {
            fprintf (stderr, "\nData set %d", dataSetNum + 1);
            fprintf (stderr, "\n\tNumber of coalescence events   =   %d", numCA);
            fprintf (stderr, "\n\tNumber of migration events     =   %d", numMIG);
        }
        // totalTreeLength = SumBranches(treeRootInit[0], programOptions.mutationRate);
        totalTreeLength = SumBranches(root, programOptions->mutationRate);
        cumNumMUperTree=0;
        
        free(newickString2);
        newickString2=NULL;
        
      
            for (z = 0; z < programOptions->MutationAssignNum; z++)
            {
                numMU=0;
                if (programOptions->doPrintSeparateReplicates == YES)
                    PrepareSeparateFilesGenotypes(1, dataSetNum, z,
                                                  filePaths, programOptions,files);
                
                //here there was the code
                if (programOptions->doSimulateData == YES)
                {
                    numISMmutations = 0;
                    numISMdeletions = 0;
                    numISMCNLOH = 0;
                }
                for (i=0; i< programOptions->numSites; i++)
                {   allSites[i].numMutations =0;
                    allSites[i].numMutationsMaternal =0;
                    allSites[i].numMutationsPaternal =0;
                }
  
                InitializeGenomes (root, &(programOptions->seed), programOptions->alphabet, programOptions->doUserGenome,programOptions->numSites,  allSites, programOptions->doGeneticSignatures,cumfreq, triNucFreq, cellNames);
                // SNPrate=0.01;
                // HEALTHY_ROOT=treeRootInit[0]->label;
                //TUMOR_ROOT=treeRootInit[0]->label;
                HEALTHY_ROOT=root->label;
                TUMOR_ROOT=root->label;
                //        if (SNPrate > 0)
                //            AddGermlineVariation (treeRootInit[0], &seed,  numSites, SNPrate, allSites, alphabet,  data,   HEALTHY_ROOT, cumMij );
                EvolveSitesOnTree (root, MATERNAL, &(programOptions->seed), programOptions->rateVarAmongSites,  programOptions->numSites,  allSites, programOptions->doGeneticSignatures, programOptions->alphaSites, programOptions->propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , &numISMmutations, programOptions->numFixedMutations, numSNVmaternal,  programOptions->doSimulateFixedNumMutations,  programOptions->alphabet,  data,  &numMU, cumMij,  programOptions->altModel, programOptions->altModelMutationRate, programOptions->doUserTree,  programOptions->doJC,  programOptions->doHKY,  programOptions->doGTR,
                                   programOptions->doGTnR,  freqR,  freqY,
                                   freqAG, freqCT, programOptions->titv, freq, Mij ,   Root,  Cijk);
                EvolveSitesOnTree (root, PATERNAL, &(programOptions->seed), programOptions->rateVarAmongSites,  programOptions->numSites,  allSites, programOptions->doGeneticSignatures, programOptions->alphaSites, programOptions->propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , &numISMmutations, numFixedMutations, numSNVmaternal,  programOptions->doSimulateFixedNumMutations,  programOptions->alphabet,  data,  &numMU, cumMij,  programOptions->altModel, programOptions->altModelMutationRate, programOptions->doUserTree,  programOptions->doJC,  programOptions->doHKY,  programOptions->doGTR,
                                   programOptions->doGTnR,  freqR,  freqY,
                                   freqAG, freqCT, programOptions->titv, freq, Mij,   Root,  Cijk );
                cumNumMU += numMU;
                cumNumMUSq += pow(numMU,2);
              
                if (programOptions->doPrintTrueHaplotypes == YES)
                {
                    if (programOptions->doPrintSeparateReplicates == NO)
                        fprintf (files->fpTrueHaplotypes, "[#%d]\n", z+1);
                    //
                    PrintTrueFullHaplotypes (files->fpTrueHaplotypes,  nodes, root , programOptions->numNodes, programOptions->doPrintIUPAChaplotypes, programOptions->doPrintAncestors, programOptions->numSites,  programOptions->numCells, programOptions->alphabet, programOptions->doUserTree , data,    programOptions->doNGS,   cellNames, cell, HEALTHY_ROOT, TUMOR_ROOT, ObservedCellNames, programOptions->doUseObservedCellNames);
                }
                
                
                if (programOptions->doPrintTrees ==YES && programOptions->doPrintSeparateReplicates == YES)
                {
                    
                    fclose(files->fpTrees);
                    fclose(files->fpTrees2);
                }
                if (programOptions->doPrintTimes ==YES && programOptions->doPrintSeparateReplicates == YES)
                {
                    
                    fclose(files->fpTimes);
                    fclose(files->fpTimes2);
                }
                if (programOptions->doPrintTrueHaplotypes ==YES)
                {
                    if (programOptions->doPrintSeparateReplicates == YES)
                        fclose(files->fpTrueHaplotypes);
                    
                }
           
            
                
            }/* end of mutation simulation process */
    
    }
    return 0;
}
/***************************** ValidateParameters*******************************/
        /* Validate parameters*/
        
void ValidateParameters(ProgramOptions *programOptions,
                                int *CloneNameBegin , int *CloneSampleSizeBegin, int *ClonePopSizeBegin)
    {
            if (programOptions->noisy > 1)
            {
                printf("\n>> Settings ..\n");
                
                printf("Number of replicates = %d\n", programOptions->numDataSets);
                if (programOptions->Nscaling == 1)
                    printf("Haploid data (%d)\n", programOptions->Nscaling);
                else
                    printf("Diploid data (%d)\n", programOptions->Nscaling);
                printf("Number of clones = %d\n", programOptions->numClones);
            }
            for (int j = 0; j < programOptions->numClones; j++)
            {
                // Checking
                //  if (numClones <  CloneNameBegin[j] )
                if (programOptions->numClones <  CloneNameBegin[j] )
                {
                    fprintf (stderr, "PARAMETER ERROR: Clon (%d) is higher than the number of clones (%d). (d)\n\n",  CloneNameBegin[j], programOptions->numClones);
                    PrintUsage();
                }
                if ( CloneSampleSizeBegin[j] > ClonePopSizeBegin[j] )
                {
                    fprintf (stderr, "PARAMETER ERROR: Clone (%d) cannot have sample size (%d) higher than population size (%d). (d)\n\n", j, CloneSampleSizeBegin[j] , ClonePopSizeBegin[j]);
                    PrintUsage();
                }
            }
        }
    /***************** bbinClones *****************/
    /* binary search in the probabilities with clones */
        int bbinClones (double dat, double *v, int n)
        {
            int init, end, middle;
            
            if (dat >= 0 && dat <= v[1])
                return (1); /* first population */
            
            init = 1;
            end = n;
            
            while (init <= end)
            {
                middle = (init + end) / 2;
                
                if (dat > v[middle - 1] && dat <= v[middle])
                    return (middle);
                else if (dat > v[middle])
                    init = middle + 1;
                else
                    end = middle - 1;
            }
            
            fprintf (stderr, "\n Warning in bbinClones function");
            exit (-1);
            return -1;
        }

/********************** InitListPossibleMigrations ************************/
        /* Initialize the list of possible migrations times in all populations */
        void  InitListPossibleMigrations(Population **populations, int numClones)
        {
            qsort(populations, numClones, sizeof(Population*), comparePopulationsByTimeOrigin);
            Population *p;
            int i, j;
            struct Population** pops ;
            double* migrationTimes ;
            double d;
            for (i = 0; i < numClones; i++) {
                p = *(populations + i);
                p->order = i;
                p->numPossibleMigrations = i + 1;
                p->numIncomingMigrations = 1; //the time of origin counts as one migration
                // if (!(p->migrationTimes)){
                p->migrationTimes =(double *)  malloc( (p->order + 1) * sizeof( double));
                if (!(p->migrationTimes)){
                    
                    fprintf (stderr, "Could not allocate p->migrationTimes (%lu bytes)\n", (p->order + 1) * (long) sizeof(double));
                    exit (1);
                    
                }
                // }
                
                d = (double)(p->timeOriginSTD);
                p->migrationTimes[0] = d;
                
                
                if(p->order >0 )//&& !(p->immigrantsPopOrderedModelTime) )
                    p->immigrantsPopOrderedModelTime = (Population**) malloc(p->order * sizeof(struct Population*)); //besides the other possible immigrants we need to add this populations itself
                //        //p->immigrantsPopOrderedModelTime[0] = p;
                
                for (j = 1; j < (p->order +1); j++)
                {
                    //p->immigrantsPopOrderedModelTime[j] = (Population *) malloc(sizeof(Population));
                    // if (!(p->immigrantsPopOrderedModelTime[j-1]))
                    p->immigrantsPopOrderedModelTime[j-1]=NULL;
                    //            p->immigrantsPopOrderedModelTime[j-1] =  malloc(sizeof(Population));
                    //            if (!( p->immigrantsPopOrderedModelTime[j-1])){
                    //
                    //                fprintf (stderr, "Could not allocate  p->immigrantsPopOrderedModelTime[j-1] (%lu bytes)\n", (p->order + 1) * (long) sizeof(double));
                    //                exit (1);
                    //
                    //            }
                    
                    d = 2 * ((double)(p->timeOriginSTD)); //  a value greater than time of origin  standarized by  the population
                    p->migrationTimes[j] = d;
                    //pops[j]->timeMigrationSTDCurrentPop= DBL_MAX;//initialize the rest of migration times to some large values
                }
            }
        }
        /********************** resetMigrationsList ************************/
        /* reset list of migrations*/
        void resetMigrationsList(Population **populations, int numClones){
            Population *p;
            int i, j;
            struct Population** pops ;
            double* migrationTimes ;
            double d;
            for (i = 0; i < numClones; i++) {
                p = *(populations + i);
                p->numIncomingMigrations = 1; //the time of origin counts as one migration
                d = (double)(p->timeOriginSTD);
                p->migrationTimes[0] = d;
                for (j = 1; j < (p->order +1); j++)
                {
                    //p->immigrantsPopOrderedModelTime[j] = (Population *) malloc(sizeof(Population));
                    // if (!(p->immigrantsPopOrderedModelTime[j-1]))
                    p->immigrantsPopOrderedModelTime[j-1]=NULL;
                    //            p->immigrantsPopOrderedModelTime[j-1] =  malloc(sizeof(Population));
                    //            if (!( p->immigrantsPopOrderedModelTime[j-1])){
                    //
                    //                fprintf (stderr, "Could not allocate  p->immigrantsPopOrderedModelTime[j-1] (%lu bytes)\n", (p->order + 1) * (long) sizeof(double));
                    //                exit (1);
                    //
                    //            }
                    d = 2 * ((double)(p->timeOriginSTD)); //  a value greater than time of origin  standarized by  the population
                    p->migrationTimes[j] = d;
                    //pops[j]->timeMigrationSTDCurrentPop= DBL_MAX;//initialize the rest of migration times to some large values
                }
            }
        }
        /********************** UpdateListMigrants************************/
        /* update the migrations times of the target population(that receives population PopOrigin)  */
        void  UpdateListMigrants(Population **populations, int numClones, Population *PopChild, Population *PopFather  ) {
            if (PopChild->index ==  PopFather->index) {
                fprintf (stderr, "\nError. The target population %d for  migration must be different than the population of origin %d \n", PopFather->index, PopChild->index);
                exit (-1);
            }
            if (PopFather->order <= PopChild->order ) {
                fprintf (stderr, "\nError. The target population %d for  migration must be older than the population of origin %d \n", PopFather->index, PopChild->index);
                exit (-1);
            }
            Population *pOrigin, *pTarget, *p;
            double *ptr;
            int i, j;
            int lengthMigrationsArray = (int)(PopFather->order) + 1;
            int updatedNumIncomingMigrations = PopFather->numIncomingMigrations;
            // printf ( "\n lengthMigrationsArray= %d \n", lengthMigrationsArray );
            
            Population **pops  = PopFather->immigrantsPopOrderedModelTime;
            ptr = PopFather->migrationTimes;
            for (j = 0; j < lengthMigrationsArray; j++) {
                if (ptr[j] >  (PopFather->timeOriginSTD) ) {
                    break; //first null pointer
                }
            }
            double updatedMigrationTime = (PopChild->timeOriginSTD) * (PopChild->effectPopSize) / (PopFather->effectPopSize);
            PopFather-> migrationTimes[j] = updatedMigrationTime;
            if(updatedNumIncomingMigrations + 1 <= PopFather->numPossibleMigrations){
                updatedNumIncomingMigrations = updatedNumIncomingMigrations + 1;
                PopFather->numIncomingMigrations = updatedNumIncomingMigrations;
                PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
            }
            //  fprintf (stderr ,"\n updatedNumIncomingMigrations %d \n",PopFather->numIncomingMigrations);
            //PopFather->immigrantsPopOrderedModelTime[j-1] = PopChild;
            //order immigrant population by time of origin
            if(PopFather->numIncomingMigrations > 1 )
                qsort(PopFather->migrationTimes, PopFather->numIncomingMigrations, sizeof(double), compare);
            if(PopFather->numIncomingMigrations -1 > 1 )
                qsort(PopFather->immigrantsPopOrderedModelTime, PopFather->numIncomingMigrations -1,  sizeof(Population*), comparePopulationsByTimeOrigin);
        }
        /********************** ChooseFatherPopulation************************/
        /* choose probabilistically  the father population of a  population  */
        Population* ChooseFatherPopulation(Population **populations, int numClones, Population  *PopChild,  long int *seed, int noisy) {
            
            Population *pOrigin, *pTarget, *p;
            double pij, ran;
            double *ptr;
            int i, j, k;
            double cumProb[numClones  - (int)(PopChild->order)];
            cumProb[0] = 0.0;
            for (j = PopChild->order + 1; j < numClones; j++)
            {
                cumProb[j - PopChild->order] = 0.0;
                p = *(populations + j);
                pij = ProbabilityCloneiFromClonej2(PopChild, p, populations, numClones);
                cumProb[j - PopChild->order] = cumProb[j - 1 - PopChild->order] + pij;
                
            }
            // now selecting the ancestral clone
            ran = RandomUniform(seed);
            //fprintf (stderr, "\n ran = %lf ", ran);
            int w = -1;
            for (k = 1; k < numClones ; k++)
            {
                if (ran <= cumProb[k])
                {
                    w = k;
                    break;
                }
            }
            Population *result =  *(populations + PopChild->order + w);
            if (noisy > 3)
                fprintf (stderr, "\nClone %d derived from clone %d\n", PopChild->index, result->index); // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
            /* Update list of migation times considering that clone i comes from clone j */
            if (noisy > 1)
                fprintf (stderr, "\n*** Updating list of migration times (considering that clone %d comes from clone %d) ..\n", PopChild->index, result->index);
            return (result); //the father population has  order  (PopChild->order) + w
        }
        /********************** AssignCurrentSequencesToPopulation************************/
        /* assign current sequences to  population  */
        void AssignCurrentSequencesToPopulation(Population **populations, TreeNode **nodes, ProgramOptions* programOptions,
                                                int numClones, int numNodes, int noisy,  int TotalNumSequences, int *numActiveGametes, int* nextAvailable,
                                                int *labelNodes, char* ObservedCellNames[], int doUseObservedCellNames, int *sampleSizes)
        {
            Population *pop;
            TreeNode *p;
            int i, j;
            double currentSampleSize;
            *numActiveGametes = 0;
            int indexFirstObservedCellName;
            
            int*  CumSamNodes = (int *) malloc ((numClones + 1)* (long) sizeof(int)); /* cumulative of samples in clones */
            if (!CumSamNodes)
            {
                fprintf (stderr, "Could not allocate CumSamNodes (%lu bytes)\n", (numNodes + 1) * (long) sizeof(int));
                exit (1);
            }
            CumSamNodes[0] = 0;
            
            for (i = 1; i <= numClones; i++)
            {
                pop = *(populations + i-1);
                CumSamNodes[i] = 0;
                currentSampleSize = pop->sampleSize;
                //pop->idsActiveGametes= (int*) malloc( (numNodes) * sizeof( int));
                pop->numCompletedCoalescences=0;
                pop->nodeIdAncesterMRCA =0;
                //pop->idsActiveGametes= (int*) malloc( (TotalNumSequences) * sizeof( int));
                pop->idsActiveGametes=NULL;
                pop->idsActiveGametes= (int*) malloc( (pop->sampleSize + pop->numPossibleMigrations) * sizeof( int));
                
                if (!( pop->idsActiveGametes))
                {
                    fprintf (stderr, "Could not allocate  pop->idsActiveGametes (%lu bytes)\n", (pop->sampleSize + pop->numPossibleMigrations ) * (long) sizeof(int));
                    exit (1);
                }
                pop->idsGametes=NULL;
                pop->idsGametes= (int*) malloc( (2* pop->sampleSize + pop->numPossibleMigrations) * sizeof( int));
                if (!( pop->idsGametes))
                {
                    fprintf (stderr, "Could not allocate  pop->idsGametes (%lu bytes)\n", (2* pop->sampleSize + pop->numPossibleMigrations) * (long) sizeof(int));
                    exit (1);
                }
                pop ->numActiveGametes=0;
                pop ->numGametes=0;
                CumSamNodes[i] = CumSamNodes[i - 1] + currentSampleSize;
            }
            if (noisy > 1)
                fprintf (stderr, "\n Initial relation nodes-clones:");
            int  cumIndivid = 0;
            int currentPopIndex;
            *numActiveGametes = 0;
            
            for (j = 0; j < TotalNumSequences; j++)
            {
                cumIndivid++;
                p = *nodes + j;
                //activeGametes[*numActiveGametes] = j;
                p->index = j;
                p->label = j;
                
                *labelNodes = j;
                p->nodeClass = 1;
                // fprintf (stderr,"\nIn cumIndivid=%d j=%d", cumIndivid, j);
                for (i = 1; i <= numClones; i++)
                {
                    pop = *(populations + i - 1);
                    currentPopIndex = pop->index;
                    indexFirstObservedCellName= pop->indexFirstObservedCellName;
                    // Identify to which clone belongs this node (sample)
                    if (cumIndivid <= CumSamNodes[i] && cumIndivid > CumSamNodes[i - 1])
                    {
                        //fprintf (stderr,"\ncumIndivid=%d <= CumSamNodes[i]=%d, in clone %d\n", cumIndivid, CumSamNodes[i], pop->index);
                        pop->idsActiveGametes[pop->numActiveGametes]=j;
                        
                        
                        pop->idsGametes[pop->numGametes]=j;
                        pop->numGametes=pop->numGametes+1;
                        
                        p->indexOldClone = currentPopIndex;
                        p->indexCurrentClone = currentPopIndex;
                        p->effectPopSize= pop->effectPopSize;
                        p->orderCurrentClone = pop->order;
                        
                        
                        
                        if (programOptions->doUseObservedCellNames)
                            strcpy( p->observedCellName,ObservedCellNames[indexFirstObservedCellName + pop->numActiveGametes ]);
                        pop->numActiveGametes=pop->numActiveGametes+1;
                        
                        break;
                    }
                }
                
                //        if(doUseObservedCellNames == YES)
                //            strcpy( p->observedCellName,ObservedCellNames[j]);
                
                
                if (noisy > 1)
                    fprintf (stderr,"\n > The node %d(%d) belongs to clone %d", p->index, cumIndivid, p->indexOldClone);
                *numActiveGametes = *numActiveGametes + 1;
            }
            //AssignObservedCellNamestoTips(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
            // AssignObservedCellNamestoTips2(nodes, populations, sampleSizes,ObservedCellNames,  programOptions);
            free(CumSamNodes);
            CumSamNodes=NULL;
            *nextAvailable = *numActiveGametes;
            *labelNodes = *labelNodes + 1;
        }
        /********************** ChooseRandomIndividual************************/
        /* ChooseRandomIndividual  */
        void ChooseRandomIndividual(int *firstInd,   int numClones, Population *popI,  int *secondInd, long *seed, int choosePairIndividuals)
        {
            double random;
            int k, w;
            double *cumPopulPart = (double *) malloc((popI->numActiveGametes + 1)* (long) sizeof(double));
            if (!cumPopulPart)
            {
                fprintf (stderr, "Could not allocate cumPopulPart (%lu bytes)\n", (popI->numActiveGametes + 1) * (long) sizeof(double));
                exit (-1);
            }
            cumPopulPart[0] = 0;
            for (k = 1; k <= popI->numActiveGametes; k++)
                cumPopulPart[k] = 0;
            for (k = 1; k <= popI->numActiveGametes; k++)
                cumPopulPart[k] = cumPopulPart[k - 1] + 1.0 / (popI->numActiveGametes);
            
            w = 0;
            //    for (k = 0; k < numClones; k++)
            //    {
            //        pop = *(populations + k);
            //        w = w + (pop->numActiveGametes);
            //        //fprintf (stderr, "\nClone %d with %d gametes. w=%d ", k, pop->numActiveGametes, w);
            //    }
            //    if (w != numActiveGametes)
            //    {
            //        fprintf (stderr, "\nError. The sum of partial active gametes is different to the total number of gametes, w %d != numActiveGametes %d. In Coalescence.", w, numActiveGametes);
            //        exit (-1);
            //    }
            random = RandomUniform(seed);
            //fprintf (stderr, "\nran = %lf ", ran);
            *firstInd = bbinClones(random, cumPopulPart, popI->numActiveGametes)-1;
            w = 0;
            
            if (*firstInd >= popI->numActiveGametes || *firstInd < 0 ) /* checking */
            {
                fprintf (stderr, "\n\nERROR: firstInd out of range!\n");
                exit (-1);
            }
            
            if (choosePairIndividuals== YES && popI->numActiveGametes > 1) {
                
                do//choose randomly another individual to coalesce
                {
                    random = RandomUniform(seed);
                    *secondInd = bbinClones(random, cumPopulPart, popI->numActiveGametes)-1;
                    
                } while (*firstInd == *secondInd  );
            }
            free (cumPopulPart);
            cumPopulPart=NULL;
        }
        
        /********************** MakeCoalescenceEvent************************/
        /*  choose 2  active individuals  to make coalescent  */
        void MakeCoalescenceEvent(Population **populations, Population *popI, TreeNode **nodes, int numClones, long int* seed, int noisy,   int *numActiveGametes, int* nextAvailable,
                                  int*labelNodes, double *currentTime, int *numNodes)
        {
            int k, w;
            double rand, ran;
            TreeNode  *p, *q, *r;
            int firstInd, i, j, secondInd=0, newInd=0;
            int choosePairIndividuals = YES;
            
            ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
            
            newInd = *nextAvailable;
            if (noisy > 1)
                fprintf (stderr, "Coalescence involving %d and %d to create node %d (in clone %d)", popI->idsActiveGametes[firstInd], popI->idsActiveGametes[secondInd], newInd, popI->index);
            /*  set pointers between nodes */
            p = *nodes + popI->idsActiveGametes[firstInd];
            q = *nodes + popI->idsActiveGametes[secondInd];
            r = *nodes + newInd;    /* new ancester */
            r->index = *nextAvailable;
            r->label = *labelNodes;
            *labelNodes=*labelNodes+1;
            r->indexOldClone = r->indexCurrentClone = popI->index;//here the clone number is updated
            // r->indexCurrentClone = p->indexCurrentClone;
            // r->orderCurrentClone = p->orderCurrentClone;
            r->orderCurrentClone = popI->order;
            r->effectPopSize=popI->effectPopSize;
            r->nodeClass = 4;
            // link the nodes
            r->left = p;
            r->right = q;
            p->anc1 = r;
            q->anc1 = r;
            r->time = *currentTime;
            r->timePUnits = *currentTime * (popI->effectPopSize);
            
            //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
            if (noisy > 1)
                fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
            /* readjust active nodes */
            
            popI->idsActiveGametes[firstInd] = newInd;
            popI->idsActiveGametes[secondInd] = popI->idsActiveGametes[popI->numActiveGametes - 1];;
            *numActiveGametes = *numActiveGametes - 1; /* less 1 active node */
            
            //update list ids nodes
            popI->idsGametes[popI->numGametes] = newInd;
            popI->numGametes = popI->numGametes +1;
            
            *nextAvailable=*nextAvailable+1; /* 1 node more is available */
            
            popI->CoalescentEventTimes[ popI->numCompletedCoalescences]=  r->time;
            popI->numActiveGametes = popI->numActiveGametes - 1; /* now this clone
                                                                  has 1 less node */
            
            popI->numCompletedCoalescences= popI->numCompletedCoalescences+1;
            /* memory for number of nodes */
            if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
            {
                /* ReallocNodes(&numNodes, activeGametes); */
                if (noisy == 4)
                    fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                *numNodes += INCREMENT_NODES;
                /* realloc */
                *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
                if (!(*nodes))
                {
                    fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
                    exit (-1);
                }
                popI->idsActiveGametes = (int *) realloc (popI->idsActiveGametes, *numNodes * (long) sizeof(int));
                if (!(popI->idsActiveGametes))
                {
                    fprintf (stderr, "Could not reallocate idsActiveGametes for the current population(%lu bytes)\n", *numNodes * (long) sizeof(int));
                    exit (-1);
                }
            }
        }
        /********************** BuildTree************************/
        /*  build tree */
        void BuildTree(Population **populations,Population *CurrentPop,
                       long int *seed,
                       ProgramOptions *programOptions,
                       TreeNode    **nodes,
                       TreeNode   **treeTips,
                       TreeNode    **treeRootInit,
                       int *nextAvailable,
                       int *newInd,
                       double *currentTime,
                       int *labelNodes
                       )
        {
            int i, j, k;
            int indexCurrentTip;
            int  foundSuperflousNode;
            TreeNode *p, *q, *r;
            /********** BUILDING TREES ***********/
            if (programOptions->noisy > 1)
            {
                fprintf (stderr, "\n>> Building trees ..");
            }
            
            j=0;
            indexCurrentTip=0;
            /* get rid of superflous nodes */
            foundSuperflousNode = YES;
            while (foundSuperflousNode == YES)
            {
                foundSuperflousNode = NO;
                for (i = 0; i < *nextAvailable; i++) // available all the nodes
                {
                    p = *nodes + i;
                    
                    //fprintf (stderr, "\n\np->index = %d", p->index);
                    
                    if (p->left == NULL && p->right == NULL && p->anc1 == NULL)
                    {
                        // nothing to do with this node because it is not connected to anything
                        //fprintf (stderr, "\n * nothing to do with this node because it is not connected to anything");
                    }
                    else if (p->left == NULL && p->right == NULL && p->anc1 != NULL)
                    {
                        // (*treeTips)[indexCurrentTip]=p;
                        if(indexCurrentTip <  programOptions->TotalNumSequences){
                            treeTips[indexCurrentTip]=p;
                            indexCurrentTip++;
                        }
                        connectNodes(NULL, NULL, p);
                        // do not do anything with this node because it is a tip
                        //fprintf (stderr, "\n * do not do anything with this node because it is a tip");
                    }
                    else if (p->left != NULL && p->right == NULL && p->anc1 != NULL)
                    {
                        // this is a superflous node and can be removed(this superfluos nodes are the MRCA nodes of the demes
                        foundSuperflousNode = YES;
                        q = p->left;
                        r = p->anc1;
                        if (p->anc1->left == p)  // p->anc up, p->left down, total: up and down for left
                        {
                            r->left = q;
                            q->anc1 = r;
                            p->left = NULL;
                            p->anc1 = NULL;
                            
                            connectNodes(q, r->right, r);
                        }
                        else
                        {
                            r->right = q;
                            q->anc1 = r;
                            p->left = NULL;
                            p->anc1 = NULL;
                            connectNodes(r->left, q, r);
                        }
                        
                        //fprintf (stderr, "\n - this is a superflous node and can be removed (1)");
                    }
                    else if (p->left == NULL && p->right != NULL && p->anc1 != NULL)
                    {
                        // this is a superflous node and can be removed
                        foundSuperflousNode = YES;
                        q = p->right;
                        r = p->anc1;
                        if (p->anc1->left == p)
                        {
                            r->left = q;
                            q->anc1 = r;
                            p->right = NULL;
                            p->anc1 = NULL;
                            connectNodes(q, r->right, r);
                        }
                        else
                        {
                            r->right = q;
                            q->anc1 = r;
                            p->right = NULL;
                            p->anc1 = NULL;
                            connectNodes(r->left, q, r);
                        }
                        
                        //fprintf (stderr, "\n - this is a superflous node and can be removed (2)");
                    }
                    else if (p->left != NULL && p->right != NULL && p->anc1 != NULL)
                    {
                        connectNodes(p->left, p->right, p);
                        // this is an internal node formed by a coalescence event, do not touch
                        //fprintf (stderr, "\n * this is an internal node formed by a coalescence event, do not touch");
                    }
                    else if (p->left != NULL && p->right != NULL && p->anc1 == NULL)
                    {
                        connectNodes(p->left, p->right, p);
                        // this is the last (coalescence event) in the tree, MRCA
                        //fprintf (stderr, "\n * this is the last (coalescence event) in the tree, MRCA");
                    }
                    else if (p->left != NULL && p->right == NULL && p->anc1 == NULL)
                    {
                        // Seems to be the last coalescent event among sequences with non-ancestral material
                        // it is not superfluous, we just remove it
                        p->left->anc1 = NULL;
                        //fprintf (stderr, "\n - this is a superflous node and can be removed (3)");
                    }
                    else if (p->left == NULL && p->right != NULL && p->anc1 == NULL)
                    {
                        // not clear what this node could be doing, but we will remove it anyway
                        fprintf (stderr, "strange\n");
                        p->left = NULL;
                        p->right->anc1 = NULL;
                        //fprintf (stderr, "\n - this is a superflous node and can be removed (4)");
                    }
                    else
                    {
                        fprintf (stderr, "You should not be here, I think\n");
                        fprintf (stderr, "%d %d-- %d %d %d\n", Index(p), j, Index(p->left), Index(p->right), Index(p->anc1));
                    }
                    if (p->anc1 != NULL)
                    {//update length field
                        
                        //   p->length = p->anc1->time- p->time;
                        p->length = (p->anc1->timePUnits- p->timePUnits);
                        //*mutationRate;
                        p->lengthModelUnits = (p->anc1->time- p->time);
                        //*mutationRate;
                        setLength(p);
                    }
                }
                //fprintf (stderr, "\n");
            }//while
            
            /* about the MRCA */
            *newInd=*nextAvailable-1;
            p = *nodes + *newInd; /* because the last one event is the last coalescence */
            //fprintf (stderr, "\n\n\n>> newInd = %d\n", newInd);
            
            if (programOptions->thereisOutgroup == NO)
            {
                p = *nodes + *newInd;
                p->nodeClass = 5;
                *treeRootInit = p;
                //treeRootInit[0] = p;
                p->anc1 = NULL;
            }
            if (programOptions->thereisOutgroup == YES && programOptions->outgroupSelection > 0)  /*** Root and outgroup ***/
            {
                p = *nodes + *newInd; // MRCA
                p->nodeClass = 4;
                
                if (programOptions->noisy > 1)
                    fprintf (stderr, "\n\n>> Attaching outgroup .. ");
                
                //fprintf (stderr, "\n>> ThisCloneNumber = %d, ListMigrationTimesInitial[ThisCloneNumber] = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf \n", ThisCloneNumber, ListMigrationTimesInitial[ThisCloneNumber], ClonePopSizeMeffectBegin[ThisCloneNumber]);
                
                if (programOptions->outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
                    *currentTime = CurrentPop->timeOriginSTD; // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
                else if (programOptions->outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
                    *currentTime = CurrentPop->timeOriginSTD + (programOptions->outgroupBranchLength_Root1Root2 / CurrentPop->effectPopSize) ; // origin of the clone + time given by the user
                    
                }
                else
                {
                    fprintf (stderr, "\n\nError simulationg the outgroup. Check input settings\n");
                    PrintUsage();
                }
                TreeNode*       healthyRoot = *nodes + *nextAvailable;
                healthyRoot->index = *nextAvailable;
                healthyRoot->label = *labelNodes;
                healthyRoot->effectPopSize= p->effectPopSize;
                *labelNodes=*labelNodes+1;
                healthyRoot->left = p;//coalTreeMRCA;
                //        coalTreeMRCA->anc = healthyRoot;
                p->anc1 = healthyRoot;
                
                healthyRoot->timePUnits = p->timePUnits * healthyRoot->effectPopSize;
                healthyRoot->nodeClass = 5;
                //        coalTreeMRCA->length = transformingBranchLength/mutationRate;
                p->length = 0;
                
                //        coalTreeMRCA->branchLength = transformingBranchLength;
                p->lengthModelUnits = 0;
                
                //        healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
                healthyRoot->time = *currentTime  ;
                
                int transformingBranchLength=1.001;
                // healthyRoot->time = p->time * transformingBranchLength ;
                healthyRoot->timePUnits = *currentTime * healthyRoot->effectPopSize;
                p->length = (p->anc1->timePUnits- p->timePUnits);
                //*mutationRate;
                p->lengthModelUnits = (p->anc1->time- p->time);
                //*mutationRate;
                
                healthyRoot->length = 0;
                //        healthyRoot->length = 0;
                
                //        if (noisy > 2)
                //            fprintf (stderr, "DONE");
                //
                (*nextAvailable)++;
                //
                //        /* connect the healthy ancestral cell with the tip healthy cell*/
                //        if (noisy > 2)
                //            fprintf (stderr, "\n>> Adding healthy tip ... ");
                TreeNode* healthyTip = *nodes + *nextAvailable;
                healthyTip->left = NULL;
                healthyTip->right = NULL;
                healthyTip->effectPopSize= healthyRoot->effectPopSize;
                
                connectNodes(NULL, NULL, healthyTip);
                
                healthyTip->anc1 = healthyRoot;
                healthyRoot->right = healthyTip;
                healthyTip->time = 0;
                healthyTip->timePUnits = 0;
                double  healthyTipBranchLengthRatio =1;
    
                healthyTip->length = (healthyTip->anc1->timePUnits- healthyTip->timePUnits);
     
                healthyTip->lengthModelUnits = (healthyTip->anc1->time- healthyTip->time);
    
                healthyTip->isOutgroup= YES;
                
                connectNodes(p, healthyTip, healthyRoot);
                setLength(p);
                setLength(healthyTip);

                *treeRootInit=healthyRoot;
      
            }
            
 
            int intLabel = 0;
            if (programOptions->noisy > 1)
                fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
            if (programOptions->thereisOutgroup == YES)
                intLabel = programOptions->TotalNumSequences + 1;
            else
                intLabel = programOptions->TotalNumSequences;
        
            RelabelNodes(*treeRootInit, treeRootInit, &intLabel );
            
        }
        /********************************** EvolveSitesOnTree ***********************************/
        /* Evolves all sites (maternal and paternal genomes) on the given tree
         We assume that a site will be ISM or Mk in both maternal and paternal genome
         */
        void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, SiteStr* allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, int* DefaultModelSites, int* AltModelSites,  double totalTreeLength , int *numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  int ***data,  int  *numMU, double cumMij[4][4], int altModel, double altModelMutationRate, int doUserTree,int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4],  double Root[], double Cijk[])
        {
            int i = 0;
            

            
            if (rateVarAmongSites == YES)
                for (i=0; i<numSites; i++)
                    allSites[i].rateMultiplier = RandomGamma (alphaSites, seed) / alphaSites;
      
            if (propAltModelSites == 0)/* only default model (ISM diploid) sites */
            {
                numDefaultModelSites = numSites;
                numAltModelSites = 0;
                for (i=0; i<numSites; i++)
                    DefaultModelSites[i] = i;
                SimulateISM (treeRoot, genome, NO, seed,
                             DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,altModelMutationRate
                             );
                
                
            }
            else if (propAltModelSites == 1)
            {
                numDefaultModelSites = 0;
                numAltModelSites = numSites;
                for (i=0; i<numSites; i++)
                    AltModelSites[i] = i;
                
                if (altModel == ISMhap)
                { fprintf(stderr, "only non ISM sites");
                    SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,    altModelMutationRate);
                }
                else if (altModel == Mk2)
                { SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
                }
                else if (altModel == finiteDNA)
                {
                    SimulateFiniteDNA (treeRoot, genome, seed,doJC,  doHKY, doGTR,  doGTnR,  freqR, freqY,  freqAG,  freqCT, titv,  freq, Mij, numAltModelSites, AltModelSites, allSites,   rateVarAmongSites, altModelMutationRate,  numMU,   Root,  Cijk);
                    
                }
                else
                {
                    fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
                    exit(-1);
                }
                
            }
            else {/* both ISM and non-ISM sites */
                numDefaultModelSites = 0;
                numAltModelSites = 0;
                
                for (i=0; i<numSites; i++)
                {
                    if (RandomUniform (seed) < propAltModelSites)
                        AltModelSites[numAltModelSites++] = i;
                    else
                        DefaultModelSites[numDefaultModelSites++] = i;
                }
                //            /* Evolve ISM sites */
                if (numDefaultModelSites > 0)
                {
                    SimulateISM (treeRoot, genome, NO, seed,
                                 DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,altModelMutationRate
                                 );
                    
                }
                //            /* Evolve non-SIM sites */
                if (numAltModelSites > 0)
                {
                    if (altModel == ISMhap)
                    { fprintf(stderr, "pip");
                        SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,   data, allSites, numMU, cumMij,    altModelMutationRate);
                    }
                    else if (altModel == Mk2)
                    { SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
                    }
                    else if (altModel == finiteDNA)
                    {
                        SimulateFiniteDNA (treeRoot, genome, seed,doJC,  doHKY, doGTR,  doGTnR,  freqR, freqY,  freqAG,  freqCT, titv,  freq, Mij, numAltModelSites,AltModelSites,  allSites,    rateVarAmongSites,  altModelMutationRate, numMU,   Root,  Cijk);
                    }
                    else
                    {
                        fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
                        exit(-1);
                    }
                }
            }
         
        }
        void SimulateISM (TreeNode *treeRoot, int genome, int doISMhaploid, long int *seed,  int *DefaultModelSites, int numDefaultModelSites, int* AltModelSites, int numAltModelSites, double totalTreeLength , int *numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  int ***data, SiteStr* allSites, int  *numMU, double cumMij[4][4], double mutationRate)
        {
            int   i, trials, numMutations, mutationsSoFar;
            double  totalBranchSum;
            int   *modelSites, numModelSites;
            double cumBranchLength =0;
            double uniform =0;
            int mutationAdded;
            double ran=0;
          
            modelSites = DefaultModelSites;
            numModelSites = numDefaultModelSites;
          
            totalBranchSum = totalTreeLength *numModelSites;
        
            int ploidyFactor;
            
            if (doISMhaploid==NO){
                ploidyFactor=1;
            }
            else{
                ploidyFactor=2;
            }
            
            do{
                numMutations = RandomPoisson (totalBranchSum, seed);
            }while( *numISMmutations + numMutations > (ploidyFactor *numModelSites));
            
         
            *numISMmutations = *numISMmutations + numMutations;
            trials = 0;
            mutationsSoFar = 0;
            mutationAdded=NO;
            while (mutationsSoFar < numMutations)
            {
                i = RandomUniformTo(numModelSites, seed);
                
     
                while (((doISMhaploid == NO)  && (allSites[i].numMutations != 0))  ||
                       ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[i].numMutationsMaternal != 0)) ||
                       ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[i].numMutationsPaternal != 0)))
                {
                    //            fprintf (stderr, "\n\n site %d again!", i);
                    i = RandomUniformTo(numModelSites, seed);
                    if (trials++ > 1000*numModelSites)
                    {
                        fprintf (stderr, "\n\n ERROR: after %d trials cannot find an unmuted site",1000*numModelSites);
                        fprintf (stderr, "\nmutations = %d   mutations so far = %d\n\n",numMutations, mutationsSoFar);
                        for (i=0; i<numModelSites; i++)
                            fprintf (stderr, "\nsite %d has %d mutations",i,allSites[i].numMutations);
                        exit(-1);
                    }
                }
                if (alphabet == DNA)
                    SimulateISMDNAforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, data, allSites, numMU,cumMij, mutationRate,&uniform, &cumBranchLength,  &ran);
                else{
                    mutationAdded=NO;
                    SimulateISMforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, data, allSites, numMU,cumMij,mutationRate ,    &cumBranchLength,  &uniform,  &mutationAdded);
                }
                mutationsSoFar++;
#ifdef MYDEBUG
                fprintf (stderr, "\nmutations = %d   mutations so far = %d\n",numMutations, mutationsSoFar);
                fprintf (stderr, "\n position = %d ",i);
                if (allSites[i].numMutations > 1)
                {
                    fprintf (stderr, "\n\n ERROR: %d mutations in site %d",allSites[i].numMutations, i);
                    exit(-1);
                }
                for (i=0; i<numDefaultModelSites; i++)
                {
                    if (allSites[i].numMutations > 0)
                        fprintf (stderr, "%2d[%d,%d,%d ] ",i, allSites[i].numMutations, allSites[i].numMutationsMaternal, allSites[i].numMutationsPaternal);
                }
#endif
                
            }
            }
        /***************************** PrintTrueFullHaplotypes *******************************/
        /* Prints observed/ML haplotypes for all sites (variable + invariable) to a file */
        
        void PrintTrueFullHaplotypes (FILE *fp, TreeNode* nodes,TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int ***data,   int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName)
        {
            int         i, j;
            char *temp;
            TreeNode *p;
          
            TreeNode * healthyTip= getHealthyTip(treeRoot);
            if (alphabet == DNA)
            {
                
                if (doPrintIUPAChaplotypes == YES)
                {
                    fprintf (fp,"%d %d\n",numCells +1, numSites);
                    for (i=0; i<numCells; i++){
                        p = (nodes + i);
                        /* print IUPAC haplotype */
                        if (p !=NULL){
                            
                            if (p->left==NULL && p->right ==NULL){
                                if (doUseObservedCellName == YES)
                                    temp=p->observedCellName;
                                else
                                    temp=p->cellName;
                                fprintf (fp,"%-12s ", temp);
                                for (j=0; j<numSites; j++)
                                    fprintf (fp, "%c", WhichIUPAC(p->maternalSequence[j],p->paternalSequence[j]));
                                fprintf (fp,"\n");
                                
                            }
                        }
                    }
                    
                    fprintf (fp,"%-12s ", healthyTip->observedCellName);
                    for (j=0; j<numSites; j++)
                        fprintf (fp, "%c", WhichIUPAC(healthyTip->maternalSequence[j],healthyTip->paternalSequence[j]));
                    fprintf (fp,"\n");
                }
                else // print maternal and paternal DNA haplotypes
                {
                    fprintf (fp,"%d %d\n",2*(numCells+1), numSites);
                    for (i=0; i<numCells; i++){
                        p = (nodes + i);
                        if (p !=NULL){
                            
                            if (p->left==NULL && p->right ==NULL){
                                if (doUseObservedCellName == YES)
                                    temp=p->observedCellName;
                                else
                                    temp=p->cellName;
                                fprintf (fp,"m%-12s ", temp);
                                for (j=0; j<numSites; j++)
                                    fprintf (fp, "%c", WhichNuc(p->maternalSequence[j]));
                                fprintf (fp,"\n");
                                fprintf (fp,"p%-12s ", temp);
                                for (j=0; j<numSites; j++)
                                    fprintf (fp, "%c", WhichNuc(p->paternalSequence[j]));
                                fprintf (fp,"\n");
                            }
                        }
                    }
                    if (doUseObservedCellName == YES)
                        fprintf (fp,"m%-12s ", healthyTip->observedCellName);
                    else
                        fprintf (fp,"m%-12s ", healthyTip->cellName);
                    
                    for (j=0; j<numSites; j++)
                        fprintf (fp, "%c", WhichNuc(healthyTip->maternalSequence[j]));
                    fprintf (fp,"\n");
                    if (doUseObservedCellName == YES)
                        fprintf (fp,"p%-12s ", healthyTip->observedCellName);
                    else
                        fprintf (fp,"p%-12s ", healthyTip->cellName);
                    
                    for (j=0; j<numSites; j++)
                        fprintf (fp, "%c", WhichNuc(healthyTip->paternalSequence[j]));
                    fprintf (fp,"\n");
                  
                }
            }
            else  //print binary haplotypes
            {
                if (doPrintIUPAChaplotypes == YES) // print binary consensus haplotypes
                {
                    for (i = 0; i < numCells; i++)
                    {
                        p = (nodes + i);
                        
                        if (p->left==NULL && p->right ==NULL){
                            if (doUseObservedCellName == YES)
                                temp=p->observedCellName;
                            else
                                temp=p->cellName;
                            
                            fprintf (fp,"%-12s", temp);
                            for (j=0; j<numSites; j++)
                                fprintf (fp, "%c", WhichConsensusBinary(p->maternalSequence[j],p->paternalSequence[j]));
                            fprintf (fp,"\n");
                        }
                    }
     
                }
                else // print maternal and paternal binary haplotypes
                {
                    int i=0;
                    int numAddedTips=0;
                    fprintf (fp,"%d %d\n",(numCells+1), numSites);
                    for (i=0; i<numCells; i++){
                        p = (nodes + i);
                        if (p !=NULL){
                            
                            if (p->left==NULL && p->right ==NULL){
                                if (doUseObservedCellName == YES)
                                    temp=p->observedCellName;
                                else
                                    temp=p->cellName;
                                numAddedTips++;
                                fprintf (fp,"%-12s ", temp);
                                for (j=0; j<numSites; j++)
                                    fprintf (fp, "%c", WhichMut(p->maternalSequence[j]+p->paternalSequence[j]));
                                fprintf (fp,"\n");
                            }
                        }
                    }
                    
                    //this next part is for printing the root/healthy cell
                    if (doUseObservedCellName == YES)
                        fprintf (fp,"%-12s ", healthyTip->observedCellName);
                    else
                        fprintf (fp,"%-12s ", healthyTip->cellName);
                    
                    for (j=0; j<numSites; j++)
                        fprintf (fp, "%c", WhichMut(healthyTip->maternalSequence[j]+healthyTip->paternalSequence[j]));
                    fprintf (fp,"\n");
                    
                    
                    
                }
            }
        }
        /********************* PrepareSeparateFiles **********************/
        /* Open individual files to output results */
        
        void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths *filePaths, const ProgramOptions *programOptions,Files *files)
        {
            char File[MAX_NAME];
            char dir[MAX_NAME];
            /* contains the simulated tree in Newick format
             */
            mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
            //mkdir("Results",0);
#ifdef MAC
            strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
            strcpy (dir, "Results/");
#endif
            //strcpy (resultsDir, dir);
            if (programOptions->doPrintTrees == YES)
            {
                if (programOptions->doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->treeDir );
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d", filePaths->resultsDir, filePaths->treeDir, ChainNumber );
                    
                }
                mkdir(File,S_IRWXU);
                
                if (programOptions->doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir, filePaths->treeFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir,ChainNumber , filePaths->treeFile, paramSetNumber+1, replicate+1);
                }
                //if ((*fpTrees = fopen(File, "w")) == NULL)
                if (openFile(&files->fpTrees, File) == -1)
                {
                    fprintf (stderr, "Can't open \"%s\"\n", File);
                    exit(-1);
                }
                
                sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->treeDir);
                
                if (programOptions->doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_2_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir, filePaths->treeFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_2_%04d_%04d.tre", filePaths->resultsDir, filePaths->treeDir, ChainNumber, filePaths->treeFile, paramSetNumber+1, replicate+1);
                    
                }
                //sprintf(File,"%s/%s/%s_2_%04d.tre", resultsDir, treeDir, treeFile, replicate+1);
                //if ((*fpTrees2 = fopen(File, "w")) == NULL)
                if (openFile(&files->fpTrees2, File) == -1)
                {
                    fprintf(stderr, "Can't open %s.\n", File);
                    exit(-1);
                }
            }
            
            if (programOptions->doPrintTimes == YES)
            {
                if (programOptions->doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->timesDir);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d", filePaths->resultsDir, filePaths->timesDir, ChainNumber );
                }
                mkdir(File,S_IRWXU);
                if (programOptions->doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, filePaths->timesFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, ChainNumber, filePaths->timesFile, paramSetNumber+1, replicate+1);
                    
                }
                //if ((*fpTimes = fopen(File, "w")) == NULL)
                if (openFile(&files->fpTimes, File) == -1)
                {
                    fprintf (stderr, "Can't open \"%s\"\n", File);
                    exit(-1);
                }
                sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->timesDir);
                if (programOptions->doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_2_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, filePaths->timesFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_2_%04d_%04d.txt", filePaths->resultsDir, filePaths->timesDir, ChainNumber,filePaths->timesFile, paramSetNumber+1, replicate+1);
                    
                }
                // sprintf(File,"%s/%s/%s_2_%04d.txt", resultsDir, timesDir, timesFile, replicate+1);
                
                //if ((*fpTimes2 = fopen(File, "w")) == NULL)
                if (openFile(&files->fpTimes2, File) == -1)
                {
                    fprintf(stderr, "Can't open %s.\n", File);
                    exit(-1);
                }
            }
            
     
        }
        /********************* PrepareSeparateFilesGenotypes **********************/
        /* Open individual genotypes files to output results */
        void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
                                           const FilePaths *filePaths, const ProgramOptions *programOptions,Files *files)
        {
            char File[MAX_NAME];
            char dir[MAX_NAME];
            /* contains the simulated tree in Newick format
             */
            mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
            //mkdir("Results",0);
#ifdef MAC
            strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
            strcpy (dir, "Results/");
#endif
            //strcpy (resultsDir, dir);
            
            if (programOptions->doSimulateData == YES)
            {
                /* contains SNV genotypes for every cell */
                if (programOptions->doPrintSNVgenotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->SNVgenotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->SNVgenotypesDir, filePaths->SNVgenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
                    //if ((*fpSNVgenotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpSNVgenotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains haplotypes for variable sites for every cell */
                if (programOptions->doPrintSNVhaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->SNVhaplotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->SNVhaplotypesDir, filePaths->SNVhaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
                    //if ((*fpSNVhaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpSNVhaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains reference haplotypes (before errors)  for every cell */
                if (programOptions->doPrintTrueHaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->trueHaplotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->trueHaplotypesDir, filePaths->trueHaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    
                    
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
                    //if ((*fpTrueHaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpTrueHaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains ML haplotypes  for every cell */
                if (programOptions->doPrintMLhaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->MLhaplotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->MLhaplotypesDir, filePaths->MLhaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, MLhaplotypesDir, MLhaplotypesFile, replicate+1);
                    // if ((*fpMLhaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpMLhaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains all genotypes (variable or invariable) for every cell */
                if (programOptions->doPrintFullGenotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths
                            ->fullGenotypesDir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->fullGenotypesDir, filePaths->fullGenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullGenotypesDir, fullGenotypesFile, replicate+1);
                    //if ((*fpFullGenotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpFullGenotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                /* contains haplotypes for all sites for every cell */
                if (programOptions->doPrintFullHaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->fullHaplotypesDir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->fullHaplotypesDir, filePaths->fullHaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullHaplotypesDir, fullHaplotypesFile, replicate+1);
                    //if ((*fpFullHaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpFullHaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                /* contains reads counts and log10 normalized genotype likelihoods for every SNV and cell */
                if (programOptions->doSimulateReadCounts == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->VCFdir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->VCFdir, filePaths->VCFfile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile, replicate+1);
                    //if ((*fpVCF = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpVCF, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                /* contains reads counts for every SNV and cell */
                if (programOptions->doPrintCATG == YES)
                {
                    sprintf(File,"%s/%s", filePaths->resultsDir, filePaths->CATGdir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths->resultsDir, filePaths->CATGdir, filePaths->CATGfile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile, replicate+1);
                    //if ((*fpCATG = fopen(File, "w")) == NULL)
                    if (openFile(&files->fpCATG, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
            }
        }
        /***************************** InitPopulationsCoalescentEvents*******************************/
        /* InitPopulationsCoalescentEvents*/
        void InitPopulationsCoalescentEvents( int numClones,  Population **populations) {
            int i;
            Population *popI;
            for( i = 0 ; i < numClones; i++)
            {
                popI=*(populations + i);
                popI->CoalescentEventTimes=(double *) calloc((popI->sampleSize -1)+ numClones -1, (long) sizeof(double));
                if (!(popI->CoalescentEventTimes))
                {
                    fprintf (stderr, "CoalescentEventTimes (%lu bytes)\n", (popI->sampleSize -1 +numClones -1) * (long) sizeof(double));
                    exit (-1);
                }
            }
        }
        /********************************** InitializeGenomes ***********************************/
        /* Initialize all genomes with the reference states  */
        void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, SiteStr* allSites, int doGeneticSignatures, double cumfreq[4],double *triNucFreq, char **cellNames)
        {
            int     i, cell, anccell, site;
            double  ran;
            
            if (p != NULL)
            {
                cell = p->label;
           
                
                if (alphabet == DNA)
                {
                    //            if (p->isOutgroup == YES)
                    if (p->anc1 == NULL)
                    {
                        if (doGeneticSignatures == NO) /* initialize genome with mononucleotide frequencies */
                        {
                            for (site=0; site<numSites; site++)
                            {
                                ran = RandomUniform(seed);
                                for (i=0; i<4; i++)
                                {
                                    if (ran <= cumfreq[i])
                                    {
                                        p->maternalSequence[site]=p->paternalSequence[site]=i;
                                        
                                        allSites[site].referenceAllele =p->maternalSequence[site];//   // then allSites[site].referenceAllele hosts the reference genome
                                        break;
                                    }
                                }
                            }
                        }
                        else /* initialize genome with trinucleotide frequencies */
                        {
                            SimulateTriNucFreqGenome (cell, seed, p, alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq);
                        }
                    }
                    else
                    {
                        anccell = p->anc1->label;
                        for (site=0; site<numSites; site++)
                        {
                            p->maternalSequence[site]=p->anc1->maternalSequence[site];
                            p->paternalSequence[site]=p->anc1->paternalSequence[site];
                        }
                    }
                }
                else{
                    for (site=0; site<numSites; site++){
                        p->maternalSequence[site]=0;
                        p->paternalSequence[site]=0;
                        p->numbersMutationsUnderSubtreePerSite[site]=0;
                        p->numbersMaternalMutationsPerSite[site]=0;
                        p->numbersPaternalMutationsPerSite[site]=0;
                        
                    }
                }
                InitializeGenomes (p->left, seed,   alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, cellNames);
                InitializeGenomes (p->right, seed,  alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, cellNames);
            }
        }
        
        /********************* WhichNucChar ************************/
        /* Returns integer representation for character nucleotudes */
        
        int WhichNucChar (char nucleotide)
        {
            if (nucleotide == 'A')
                return (A);
            else if (nucleotide == 'C')
                return (C);
            else if (nucleotide == 'G')
                return (G);
            else if (nucleotide == 'T')
                return (T);
            else if (nucleotide == '?')
                return (ADO);
            else if (nucleotide == '-')
                return (DELETION);
            else if (nucleotide == 'N')
                return (N);
            else if (nucleotide == 'R')
                return (R);
            else
            {
                fprintf (stderr, "\nERROR in WhichNucChar: nucleotide = %c\n",  nucleotide);
                exit(-1);
            }
        }

        /* Returns the sum of the branch lengths for a given tree */
        double SumBranches (TreeNode *p, double mutationRate)
        {
            static double sum;
            
            if (p != NULL)
            {
                if (p->anc1 == NULL)
                    sum = 0;
                else{
                    //sum += (p->anc1->time- p->time)* mutationRate;//p->lengthModelUnits;
                    sum += p->length;
                    //            sum += p->length;
                }
                //            sum += p->lengthModelUnits;//length;
                SumBranches (p->left,  mutationRate);
                SumBranches (p->right,   mutationRate);
            }
            
            return sum;
        }
        /************************* MakeCoalescenceTree2 ************************/
        /* Builds a genealogy under the structured  ' */ /* this function go by events CA, MIG, CONV */
        
        void MakeCoalescenceTree2 (long int *seed, Population **populations,
                                   int *numNodes,
                                   int numClones,
                                   ProgramOptions *programOptions,
                                   double      cumNumCA,
                                   double meanNumCA,
                                   double cumNumMIG,
                                   double meanNumMIG,
                                   int  *numMIG,
                                   int  *numCA,
                                   double *numEventsTot,
                                   TreeNode** nodes,
                                   TreeNode** treeTips,
                                   TreeNode    **treeRootInit,
                                   char* ObservedCellNames[],
                                   int *sampleSizes
                                   ) {
            
            int      c, d, i, j, w, k, m, cumIndivid, *activeGametes = NULL, isCoalescence, whichInd,
            firstInd, secondInd, newInd, eventNum, numActiveGametes, foundSuperflousNode,
            isMigration, whichClone, currentNumberAliveClones;
            int     labelNodes;
            double    currentTime, eventTime;
            TreeNode  *p, *q, *r;
            double    ran;
            double    *cumPopulPart;
            int     *numParcialActiveGametes, *CumSamNodes;
            int         *fromClone;
            int         ThisCloneNumber, ThisOriginCloneNumber;
            double      minValue;
            double      ThisRateCA;
            double      ThisTimeCA_W;
            double      ThisTimeCA_V1;
            double      ThisTimeCA_V2;
            
            int         doAmigration;
            int         ThisCloneNumberMigrations, ThisM;
            int nextAvailable;
            
            
            /* defaults */
            isCoalescence = NO;
            isMigration = NO;
            newInd = whichClone = labelNodes = cumIndivid = 0;
            numParcialActiveGametes = NULL;
            fromClone = NULL;
            eventTime = 0.0;
            c = m = 0;
            whichInd = 0;
            minValue = 0.0;
            ThisCloneNumber = 0;
            ThisOriginCloneNumber = 0;
            ThisRateCA = 0.0;
            ThisTimeCA_W = 0.0;
            ThisTimeCA_V1 = 0.0;
            ThisTimeCA_V2 = 0.0;
            doAmigration = -1;
            ThisCloneNumberMigrations = -1;
            ThisM = -1;
            
            
            //*numNodes = 2 * TotalNumSequences * numClones+ 1; // (2 * TotalNumSequences) + numClones(superfluos) - 1, but let's allocate some more..
            currentNumberAliveClones = numClones;
            
            
            //InitListPossibleMigrations(populations, numClones);
            resetMigrationsList( populations,  numClones);
            
            //allocate memory for the treenodes
            *nodes = (TreeNode *) malloc ((programOptions->numNodes + 1)* (sizeof(TreeNode)+ 5* programOptions->numSites * sizeof(int) + 2*MAX_NAME * sizeof(char)+ 3 * sizeof(pll_unode_t) + 3*sizeof(pll_tree_edge_t) )); /* nodes */
            if (!(*nodes))
            {
                fprintf (stderr, "Could not allocate nodes (%lu bytes)\n", (programOptions->numNodes+ 1)  * (long) sizeof(TreeNode));
                exit (1);
            }
            
            for (i=0; i< programOptions->TotalNumSequences; i++){
                treeTips[i]=NULL;
            }
            
            //    treeRootInit = (TreeNode **) calloc(1, sizeof(TreeNode *)); /* nodes pointers */
            //    if (!treeRootInit)
            //    {
            //        fprintf (stderr, "Could not allocate treeRootInit (%lu bytes)\n", 1  * (long) sizeof(TreeNode));
            //        exit (1);
            //    }
            /* set everything to null */
            for (i = 0; i < *numNodes; i++)
            {
                p = (*nodes + i);
                p->left = NULL;
                p->right = NULL;
                p->anc1 = NULL;
                p->outgroup = NULL;
                p->nodeRight=NULL;
                p->nodeRight=NULL;
                p->nodeBack=NULL;
                p->time = 0;
                p->timePUnits = 0;
                p->length = 0;
                p->lengthModelUnits = 0;
                p->index = 0;
                p->label = 0;
                p->isOutgroup = NO;
                p->nodeClass = 0;
                p->indexOldClone = 0;
                p->indexCurrentClone = 0;
                p->orderCurrentClone = 0;
                p->effectPopSize=0;
                p->isLeaf = NO;
                //  p->cellName =(char *) malloc(MAX_NAME* sizeof(char));
                //        if (!(p->cellName))
                //        {
                //            fprintf (stderr, "Could not allocate p->cellName (%lu bytes)\n", MAX_NAME * (long) sizeof(char));
                //            exit (1);
                //        }
                if (programOptions->doSimulateData ==YES)
                {
                    
                    p->maternalSequence= (int*) malloc (programOptions->numSites * sizeof(int));
                    if (!(p->maternalSequence))
                    {
                        fprintf (stderr, "Could not allocate p->maternalSequence (%lu bytes)\n", programOptions->numSites* sizeof(int));
                        exit (1);
                    }
                    
                    p->paternalSequence= (int*) malloc (programOptions->numSites* sizeof(int));
                    if (!(p->paternalSequence))
                    {
                        fprintf (stderr, "Could not allocate p->paternalSequence (%lu bytes)\n", programOptions->numSites* sizeof(int));
                        exit (1);
                    }
                    p->numbersMutationsUnderSubtreePerSite=(int*) calloc (programOptions->numSites, sizeof(int));
                    if (!(p->numbersMutationsUnderSubtreePerSite))
                    {
                        fprintf (stderr, "Could not allocate p->NumbersMutationsUnderSubtreePerSite (%lu bytes)\n", programOptions->numSites* sizeof(int));
                        exit (1);
                    }
                    p->numbersMaternalMutationsPerSite=(int*) calloc (programOptions->numSites, sizeof(int));
                    if (!(p->numbersMaternalMutationsPerSite))
                    {
                        fprintf (stderr, "Could not allocate p->NumbersMaternalMutationsPerSite (%lu bytes)\n", programOptions->numSites* sizeof(int));
                        exit (1);
                    }
                    p->numbersPaternalMutationsPerSite=(int*) calloc (programOptions->numSites , sizeof(int));
                    if (!(p->numbersMaternalMutationsPerSite))
                    {
                        fprintf (stderr, "Could not allocate p->NumbersMaternalMutationsPerSite (%lu bytes)\n", programOptions->numSites* sizeof(int));
                        exit (1);
                    }
                    
                }
                
                p->nodeLeft=(pll_unode_t*) calloc (1 , sizeof(pll_unode_t));
                if (!(p->nodeLeft))
                {
                    fprintf (stderr, "Could not allocate p->nodeLeft (%lu bytes)\n", 1* sizeof(pll_unode_t));
                    exit (1);
                }
                p->nodeRight=(pll_unode_t*) calloc (1 , sizeof(pll_unode_t));
                if (!(p->nodeRight))
                {
                    fprintf (stderr, "Could not allocate p->nodeRight (%lu bytes)\n",  sizeof(pll_unode_t));
                    exit (1);
                }
                p->nodeBack=(pll_unode_t*) calloc (1 , sizeof(pll_unode_t));
                if (!(p->nodeBack))
                {
                    fprintf (stderr, "Could not allocate p->nodeBack (%lu bytes)\n",  sizeof(pll_unode_t));
                    exit (1);
                }
                
                p->edgeBack=(pll_tree_edge_t*) calloc (1 , sizeof(pll_tree_edge_t));
                if (!(p->edgeBack))
                {
                    fprintf (stderr, "Could not allocate p->edgeBack (%lu bytes)\n",  sizeof(pll_tree_edge_t));
                    exit (1);
                }
                p->edgeLeft=(pll_tree_edge_t*) calloc (1 , sizeof(pll_tree_edge_t));
                if (!(p->edgeLeft))
                {
                    fprintf (stderr, "Could not allocate p->edgeLeft (%lu bytes)\n",  sizeof(pll_tree_edge_t));
                    exit (1);
                }
                p->edgeRight=(pll_tree_edge_t*) calloc (1 , sizeof(pll_tree_edge_t));
                if (!(p->edgeRight))
                {
                    fprintf (stderr, "Could not allocate p->edgeRight (%lu bytes)\n",  sizeof(pll_tree_edge_t));
                    exit (1);
                }
            }
            AssignCurrentSequencesToPopulation(populations, nodes, programOptions, numClones, *numNodes, programOptions->noisy, programOptions->TotalNumSequences, &numActiveGametes,  &nextAvailable,
                                               &labelNodes, ObservedCellNames, programOptions->doUseObservedCellNames, sampleSizes);
            Population *currentPop;
            Population *fatherPop;
            i=0;
            currentTime=0.0;
            while (i < numClones) {
                currentPop = *(populations + i);
                SimulatePopulation(currentPop, populations,programOptions, seed,
                                   &(programOptions->numNodes),
                                   numClones,
                                   cumNumCA,
                                   meanNumCA,
                                   cumNumMIG,
                                   meanNumMIG,
                                   numMIG,
                                   numCA,
                                   numEventsTot,
                                   nodes,
                                   &nextAvailable ,
                                   &numActiveGametes,
                                   &labelNodes,
                                   &currentTime,
                                   &eventNum);
                if (i< numClones-1)   //if it is not the last one
                {
                    //choose the father population from which the population i came
                    
                    fatherPop= ChooseFatherPopulation(populations, numClones, currentPop, seed,  programOptions->noisy);
                    currentPop->FatherPop = fatherPop;
                    //update list of migrant times
                    UpdateListMigrants(populations, numClones, currentPop, fatherPop);
                }
                i = i + 1;
            }
            //    free (CumSamNodes);
            //   free (activeGametes);
            BuildTree(populations,currentPop,
                      seed,
                      programOptions,
                      nodes,
                      treeTips,
                      treeRootInit,
                      //TreeNode    **treeRootInit,
                      &nextAvailable,
                      &newInd,
                      &currentTime,
                      &labelNodes
                      );
            
            if (programOptions->noisy > 1)
                fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
            //    if (thereisOutgroup == YES)
            //        intLabel = TotalNumSequences + 2;
            //    else
            //        intLabel = TotalNumSequences;
            
            // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
        }
        /************************* SimulatePopulation ************************/
        /* simulates the evolution of population until its MRCA and can receive inmigrants  ' */ /* */
        
        void SimulatePopulation( Population *popI, Population** populations,
                                ProgramOptions *programOptions,
                                long int *seed,
                                int *numNodes,
                                int numClones,
                                double      cumNumCA,
                                double meanNumCA,
                                double cumNumMIG,
                                double meanNumMIG,
                                int  *numMIG,
                                int  *numCA,
                                double *numEventsTot,
                                TreeNode    **nodes,
                                int *nextAvailable,
                                int*  numActiveGametes,
                                int* labelNodes,
                                double *currentTime,
                                int* eventNum
                                )
        {
            int  c, d, i, j, w, k, m, cumIndivid, isCoalescence, whichInd,
            firstInd, secondInd, newInd,  foundSuperflousNode,
            isMigration, whichClone, currentNumberAliveClones;
            double     eventTime;
            TreeNode  *p, *q, *r;
            double    ran;
            double    *cumPopulPart;
            *eventNum = 0;
            int numSimClones = 0;
            double ThisRateCA = 0.0;
            double ThisTimeCA_W = 0.0;
            double ThisTimeCA_V1 = 0.0;
            double ThisTimeCA_V2 = 0.0;
            int numParcialActiveGametes;
            int choosePairIndividuals;
            
            numParcialActiveGametes = popI->numActiveGametes;
            //    double CloneDelta = popI->delta;
            //    double CloneTimeOriginInputSTDPOP = popI->timeOriginSTD;
            //int numLeftMigrations=popI->numIncomingMigrations;
            int numMigrations = (popI->numIncomingMigrations); //taking into account also the    time of origin as a migration
            
            //    int ThisCloneNumber = popI->index;
            double timeNextMigration;
            int indexNextMigration = 0;
            Population *incommingPop;
            //    Population *p2;
            //fprintf (stderr, "\n\n> numMigrations= %d \n", numMigrations);
            if (programOptions->noisy > 1)
                fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %lf)\n", popI->index, popI->numActiveGametes, popI->timeOriginInput);
            if (programOptions->noisy > 1)
                fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", popI->order);
            //fprintf (stderr, "\n\n> Simulating evolutionary history of clone %d ..\n", popI->index);
            *currentTime=0;
            while (indexNextMigration < numMigrations) {
                timeNextMigration = (double)(popI->migrationTimes)[indexNextMigration];
                //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
                if ( popI->numActiveGametes >= 2) {
                    ThisRateCA = (double)  popI->numActiveGametes * ((double)  popI->numActiveGametes - 1) / 2.0;
                    ThisTimeCA_W = RandomExponential (ThisRateCA, seed) ;
                    ThisTimeCA_V1 = FmodelTstandard (*currentTime, popI->timeOriginSTD, popI->delta);
                    ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
                    // from standard time to model time, GstandardTmodel(V, T, delta)
                    ThisTimeCA_V2 = GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD, popI->delta);
                }
                else
                {
                    ThisTimeCA_V2 = timeNextMigration + 1.0; // it reached a "provisional" MRCA
                }
                if ( ThisTimeCA_V2 < timeNextMigration)
                {
                    //choose randomly two lineages to coalesce
                    isCoalescence = YES;
                    isMigration = NO;
                    *numCA = *numCA + 1;
                    *eventNum= *eventNum +1;
                    whichClone = popI->index;
                    *currentTime = ThisTimeCA_V2; // update current time in model time
                    eventTime = *currentTime;
                    
                    if (programOptions->noisy > 1)
                    {
                        fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %lf, currentTime (standard time) = %lf\n", *eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                        fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %lf\n", *eventNum, ThisTimeCA_V2);
                    }
                    if (programOptions->noisy == 4)
                        fprintf (stderr, "* Coalescence *\n");
                    MakeCoalescenceEvent(populations, popI, nodes, numClones, seed, programOptions->noisy, numActiveGametes, nextAvailable,
                                         labelNodes, currentTime,  &(programOptions->numNodes));
                }
                else
                {
                    if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
                    {
                        isCoalescence = NO;
                        isMigration = YES;
                        *numMIG = *numMIG + 1;
                        *eventNum= *eventNum +1;
                        *currentTime = timeNextMigration;
                        eventTime = *currentTime;
                        //update migration times in model time
                        if (programOptions->noisy > 1)
                        {
                            fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %lf\n", *eventNum, ThisTimeCA_V2);
                        }
                        if (programOptions->noisy == 4)
                        {fprintf (stderr, "* Migration *\n");}
                        newInd = *nextAvailable;
                        
                        r = *nodes + newInd;   /* new ancester */
                        r->index = *nextAvailable;
                        r->label = *labelNodes;
                        *labelNodes=*labelNodes+1;
                        
                        r->indexCurrentClone = popI->index;
                        r->indexCurrentClone = popI->index;
                        r->orderCurrentClone = popI->order;
                        r->nodeClass = 4;
                        
                        
                        //    p = *nodes + MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration]; // root of younger clone
                        incommingPop = *((popI->immigrantsPopOrderedModelTime) + indexNextMigration );
                        
                        if (!incommingPop){
                            fprintf (stderr, "\nError. The incoming population to  poulation %d is empty \n", popI->index);
                            exit (-1);
                        }
                        //r->indexOldClone = incommingPop->index;
                        p = *nodes + (incommingPop->nodeIdAncesterMRCA); // root of younger clone
                        indexNextMigration = indexNextMigration + 1;
                        p->indexCurrentClone = popI->index;
                        p->indexOldClone = incommingPop->index;
                        p->orderCurrentClone = popI->order;
                        // link the nodes
                        r->left = p;
                        
                        r->right = NULL;
                        //choosePairIndividuals = NO;
                        
                        //ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
                        //q=*nodes + firstInd;
                        //r->right = q;//choose another random living individual of the population
                        
                        p->anc1 = r;
                        //q->anc1 = r;
                        
                        //connectNodes(p, NULL, r);
                        //p->time = *currentTime;
                        // p->timePUnits = *currentTime * (popI->effectPopSize);
                        
                        r->time = *currentTime;// this is equal to the time of the migration
                        r->timePUnits = *currentTime * (popI->effectPopSize);
                        *nextAvailable=*nextAvailable+1; /* 1 node more is available */
                        
                        k = p->indexCurrentClone;
                        incommingPop->numActiveGametes = incommingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                        // remove node from old clone in list of active gametes and add the new node of the current clone
                        //popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                        popI->idsActiveGametes[popI->numActiveGametes]=p->index;//adding the superfluos node
                        popI->numActiveGametes = popI->numActiveGametes + 1; /* now this clone has 1 more node */
                        //                if (noisy > 1)
                        //                    fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, popI->index, incommingPop->nodeIdAncesterMRCA, k);
                        //fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, ThisCloneNumber, MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration], k);
                        if (programOptions->noisy > 1)
                            fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", p->timePUnits);
                        // fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                        /* memory for number of nodes */
                        if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
                        {
                            /* ReallocNodes(&numNodes, activeGametes); */
                            if (programOptions->noisy == 4)
                                fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                            *numNodes += INCREMENT_NODES;
                            /* realloc */
                            *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                            if (!(*nodes))
                            {
                                fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
                                exit (-1);
                            }
                            //                    activeGametes = (int *) realloc (activeGametes, *numNodes * (long) sizeof(int));
                            //                    if (!activeGametes)
                            //                    {
                            //                        fprintf (stderr, "Could not reallocate activeGametes (%lu bytes)\n", *numNodes * (long) sizeof(int));
                            //                        exit (-1);
                            //                    }
                        }
                        
                    }
                    else {
                        //origin reached
                        *currentTime = timeNextMigration;
                        indexNextMigration = indexNextMigration + 1;
                        isCoalescence = NO;
                        isMigration = NO;
                        eventTime = *currentTime;
                        if (programOptions->noisy > 1)
                            fprintf (stderr, "\n\n*** Clone origin ***\n");
                        if (programOptions->noisy == 4)
                            fprintf (stderr, "Clone origin %d at time (model units) = %lf\n", popI->index, *currentTime);
                        if (popI->order < numClones - 1) // do not do it for the last clone
                        {
                            newInd = *nextAvailable;
                            r = *nodes + newInd;    /* new ancestor */
                            r->index = *nextAvailable;
                            r->label = *labelNodes;
                            *labelNodes=*labelNodes+1;
                            r->indexOldClone =r->indexCurrentClone = popI->index;
                            r->orderCurrentClone = popI->order;
                            r->effectPopSize=popI->effectPopSize;
                            popI->nodeIdAncesterMRCA=newInd;
                            
                            r->nodeClass = 4;
                            
                            firstInd = *nextAvailable - 1;
                            //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                            p = *nodes + firstInd; // descendant node (previously generated node, nextAvailable - 1)
                            // link the nodes
                            r->left = p;
                            r->right = NULL;
                            p->anc1 = r;
                            r->time = *currentTime;
                            r->timePUnits = *currentTime * popI->effectPopSize;
                            popI->MRCA = p;
                            
                            //connectNodes(p, NULL, r);
                            //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                            /* readjust active nodes */
                            *nextAvailable=*nextAvailable+1; /* 1 node more is available */
                            popI->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                            
                            //popI->idsGametes[popI->numGametes] = newInd; r is a superflous node and it will be removed so no need to add it
                            //popI->numGametes = popI->numGametes +1;
                            
                            if (programOptions->noisy > 1)
                                fprintf (stderr, "Creating origin node, it creates node %d derived from node %d", newInd, firstInd);
                            if (programOptions->noisy > 1)
                                fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                            /* memory for number of nodes */
                            if (*nextAvailable >= *numNodes)  /* if there aren't enough nodes it go into and it addition more */
                            {
                                /* ReallocNodes(&numNodes); */
                                if (programOptions->noisy == 4)
                                    fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                                *numNodes += INCREMENT_NODES;
                                /* realloc */
                                *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
                                if (!(*nodes))
                                {
                                    fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
                                    exit (-1);
                                }
                                
                            }
                        }
                        else  {//origin of oldest pop reached
                            popI->nodeIdAncesterMRCA=*nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                            r = *nodes + *nextAvailable-1;//popI->idsActiveGametes[0]
                            r->indexOldClone = r->indexCurrentClone = popI->index;
                            r->orderCurrentClone = popI->order;
                            popI->MRCA= r;
                            
                            
                        }
                    }
                }
                if (programOptions->noisy > 3)
                {
                    fprintf (stderr, "\nActive nodes (%d):",  popI->numActiveGametes);
                    for (i = 0; i < popI->numActiveGametes; i++)
                        fprintf (stderr, " %d", popI->idsActiveGametes[i]);
                    fprintf (stderr, "\t|\tNext node available = %d", *nextAvailable);
                }
            }
            if (programOptions->noisy > 1)
                fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", popI->index);
        }
        /**************** RelabelNodes **************/
        /*  After getting rid of superfluos node, we
         need to relabel those so they are consecutive
         Use the indexes as labels when there
         is recombination */
        void RelabelNodes(TreeNode *p, TreeNode **treeRootInit, int *intLabel)
        {
            if (p != NULL)
            {
                RelabelNodes (p->left, treeRootInit, intLabel);
                RelabelNodes (p->right, treeRootInit, intLabel);
                /*RelabelNodes (p->outgroup);*/
                if (p->left == NULL && p->right == NULL) /* is tip */
                {
                    // p->label = intLabel++;
                    //  p->label = (*intLabel);
                    // *intLabel=*intLabel+1;
                    p->label = p->index ;
                }
                else                  /* all ancester */
                {
                    //p->label = intLabel++;
                    p->label = (*intLabel);
                    *intLabel=*intLabel+1;
                }
            }
        }
        
        
        /***************** Index ***************/
        /* Returns index for a given node */
        int Index (TreeNode *p)
        {
            //return (p == NULL) ? -1 : p->index+1; /* If the node haven't got bond => index = -1, else index = index+1 */
            return (p == NULL) ? -1 : p->index; /* If the node haven't got bond => index = -1, else index = index */
        }
        
        
        /***************** Lab ***************/
        /* Returns label for a given node */
        int Label (TreeNode *p)
        {
            return (p->anc1 == NULL && p->left == NULL && p->right == NULL) ? -1 : p->label + 1; /* If the node haven't got ancester and descendants => label = -1, else label = label+1 */
        }
        
        
        /**************** PrintTrees ***************/
        /*  Print unrooted trees to treefile in Newick format */
        void PrintTrees(int replicate, TreeNode **treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames)
        {
            /* there isn´t recombination */
            /*fprintf(fpTrees,"Tree.%05d = ", replicate+1);*/
            //    fprintf(fpTrees, "(");
            WriteTree (treeRootInit[0], mutationRate, fpTrees, doUseObservedCellNames);
            //    fprintf(fpTrees, ");\n");
            fprintf (fpTrees,");\n");
        }
        
        /**************** PrintTrees2 ***************/
        /*  Print unrooted trees to treefile in Newick format */
        void PrintTrees2(int replicate, TreeNode **treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames)
        {
            int indexCurrentCell=0;
            
            /* there isn´t recombination */
            /*fprintf(fpTrees2,"Tree.%05d = ", replicate+1);*/
            //   fprintf(fpTrees2, "(");
            WriteTree2 (treeRootInit[0], mutationRate, fpTrees2, ObservedCellNames, &indexCurrentCell, doUseObservedCellNames);
            //     fprintf(fpTrees2, ");\n");
            //    long len= strlen(newickString);
            //    char *res = malloc(len  + strlen(");\n"));
            //    if (res){
            //        memcpy(res, newickString, len);
            //        memcpy(res + len, ");\n", strlen(");\n")+1);
            //    }
            fprintf(fpTrees2, ");\n");
            
            //fprintf (fpTrees2,"\n");
        }
        
        
        /******************* WriteTree ****************/
        /* Writes a given (unrooted) tree from PrintTrees */
        void WriteTree (TreeNode *p, double mutationRate, FILE    *fpTrees, int doUseObservedCellNames)
        {
            char buffer[1024];
            
            if (p != NULL)
            {
                if(p->isOutgroup == YES)            /* Outgroup*/
                {
                    /*            fprintf (fpTrees, ",outgroup:%8.6f)",p->length*mutationRate);*/
                    //fprintf (fpTrees, ",outgroup:%8.6f",p->length*mutationRate);
                    strcpy( p->cellName,"healthycell");
                    strcpy( p->observedCellName,"healthycell");
                    //p->cellName[MAX_NAME]=0;
                    //                fprintf (fpTrees, ",outgroup:%10.9lf",p->length);
                    //fprintf (fpTrees, ",outgroup:%10.9lf",(p->anc1->time- p->time)*mutationRate);
                    fprintf (fpTrees, "healthycell:%10.9lf",(p->anc1->time- p->time)*mutationRate);
                }
                else if (p->left == NULL && p->right == NULL)        /* tip of the tree */
                {
                    //fprintf (stderr, "\n\n>> p->index = %d, p->class = %d \n\n", p->index, p->class);
                    //fprintf (fpTrees, "samp%05d_C%dR%d:%8.6f", p->index,p->indexOldClone,p->indexOldRegion,(p->anc1->time-p->time)*mutationRate);
                    //   snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
                    snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
                    strcpy( p->cellName,buffer);
                    //strncpy( p->cellName,buffer, sizeof(p->cellName)-1);
                    // p->cellName[MAX_NAME]=0;
                    //            fprintf (fpTrees, "tip_i%05d_C%d_%d:%10.9lf", p->index,p->indexOldClone,p->indexCurrentClone,(p->anc1->time-p->time)*mutationRate);
                    //fprintf (fpTrees, "tip_i%05d_C%d_%d:%10.9lf", p->index,p->indexOldClone,p->indexCurrentClone,(p->anc1->time- p->time)*mutationRate);
                    fprintf (fpTrees, "tip_i%05d_C%d_%d:%10.9lf", p->index,p->indexOldClone,p->indexCurrentClone,(p->anc1->time- p->time)*mutationRate);
                }
                else                                /* all ancester */
                {
                    fprintf (fpTrees, "(");
                    WriteTree (p->left, mutationRate, fpTrees, doUseObservedCellNames);
                    if (p->right != NULL) // Miguel added this condition to consider an outgroup as this right node that is NULL (see add outgroup)
                    {
                        fprintf (fpTrees, ",");
                        WriteTree (p->right, mutationRate, fpTrees, doUseObservedCellNames);
                    }
                    if (p->anc1 !=NULL)
                    {
                        //fprintf (fpTrees, "):%8.6f",(p->anc1->time-p->time)*mutationRate);
                        //snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                        snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                        strcpy( p->cellName,buffer);
                        //strncpy( p->cellName,buffer, sizeof(p->cellName)-1);
                        //p->cellName[MAX_NAME]=0;
                        //                 fprintf (fpTrees, "):%10.9lf", (p->anc1->time-p->time)*mutationRate);
                        fprintf (fpTrees, "):%10.9lf", (p->anc1->time- p->time)*mutationRate);
                        
                        //                fprintf (fpTrees, ")int_i%05d_C%d_%d:%10.9lf",p->index, p->indexOldClone, p->indexCurrentClone, (p->anc1->time-p->time)*mutationRate);
                    }
                    if (p->anc1 ==NULL)  {
                        
                        //snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                        snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                        strcpy( p->cellName,buffer);
                        //strncpy( p->cellName,buffer, sizeof(p->cellName)-1);
                        //p->cellName[MAX_NAME]=0;
                        //                    fprintf (fpTrees, ")"  );
                        //                  fprintf (fpTrees, "):0.00"  );
                        //               fprintf (fpTrees, ")root_i%05d_C%d_%d:0.00", p->index,p->indexOldClone,p->indexCurrentClone );
                        
                    }
                    WriteTree (p->outgroup, mutationRate, fpTrees, doUseObservedCellNames);
                }
            }
        }
        /******************* WriteTree2 ****************/
        /* Writes a given (unrooted) tree from PrintTrees */
        void WriteTree2 ( TreeNode *p, double mutationRate, FILE    *fpTrees2, char *cellNames[], int *indexCurrentCell, int doUseObservedCellNames)
        {
            char *currentNewick = NULL;
            //asprintf(&currentNewick, *newickString);
            char *temp = NULL;
            if (p != NULL)
            {
                if (p->isOutgroup == YES)     /* Outgroup */
                {
                    /*      fprintf (fpTrees2, ",outgroup:%8.6f)",p->length*mutationRate);*/
                    //fprintf (fpTrees2, ",outgroup:%8.6f",p->length*mutationRate);
                    //fprintf (fpTrees2, ",outgroup:%10.9lf", p->length * mutationRate);
                    //fprintf (fpTrees2, "healthycell:%10.9lf", p->length * mutationRate);
                    //fprintf (fpTrees2, "healthycell:%10.9lf", p->length * mutationRate);
                    fprintf (fpTrees2, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
                }
                else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
                {
                    //fprintf (fpTrees2, "samp%05d_C%dR%d:%8.6f", p->index,p->indexOldClone,p->indexOldRegion,(p->anc1->time-p->time)*mutationRate);
                    // fprintf (fpTrees2, "tip_i%05d_C%d_%d:%10.9lf", p->index, p->indexOldClone,p->indexCurrentClone, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                    //fprintf (fpTrees2, "tip_i%05d_C%d_%d:%10.9lf", p->index, p->indexOldClone,p->indexCurrentClone, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                    if (doUseObservedCellNames == YES)
                    {
                        if (strcmp(cellNames[*indexCurrentCell],"healthycell")==0)
                            *indexCurrentCell =*indexCurrentCell+1;
                    }
                    if (doUseObservedCellNames == YES)
                        fprintf (fpTrees2, "%s:%10.9lf", p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                    else
                        fprintf (fpTrees2, "%s:%10.9lf", p->cellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                    // fprintf (fpTrees2, "%s:%10.9lf", cellNames[*indexCurrentCell], (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                    *indexCurrentCell =*indexCurrentCell+1;
                }
                else                /* all ancester */
                {
                    fprintf (fpTrees2, "(");
                    WriteTree2 (p->left, mutationRate, fpTrees2,  cellNames, indexCurrentCell, doUseObservedCellNames);
                    if (p->right != NULL) // Miguel added this condition to consider an outgroup as this right node that is NULL (see add outgroup)
                    {
                        fprintf (fpTrees2, ",");
                        WriteTree2 (p->right, mutationRate, fpTrees2,  cellNames, indexCurrentCell, doUseObservedCellNames);
                    }
                    if (p->anc1 != NULL)
                    {
                        //                //fprintf (fpTrees2, "):%8.6f",(p->anc1->time-p->time)*mutationRate);
                        //                fprintf (fpTrees2, ")int_i%05d_C%d:%10.9lf", p->index, p->indexCoalClone, (p->anc1->timePUnits - p->timePUnits)*1);
                        fprintf (fpTrees2, "):%10.9lf",  (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                    }
                    if (p->anc1 ==NULL)  {
                        //                  fprintf (fpTrees2, ")root_i%05d_C%d_%d:0.00", p->index,p->indexOldClone,p->indexCurrentClone );
                        //                 fprintf (fpTrees2, "):0.00" );
                    }
                    WriteTree2 (p->outgroup, mutationRate, fpTrees2,  cellNames, indexCurrentCell, doUseObservedCellNames);
                }
            }
        }
        /******************* toNewickString2 ****************/
        /*  Build a string with the Newick representation using the left, right, anc1 pointers  */
        char * toNewickString2 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames)
        {
            char buffer[1024];
            //char *newickString = malloc( size);
            char *newickString =NULL;
            char *left=NULL;
            char *right=NULL;
            char *outgroup =NULL;
            char *temp =NULL;
            if (p != NULL)
            {
                if (p->isOutgroup == YES)     /* Outgroup */
                {
                    strcpy( p->cellName,"healthycell");
                    
                    //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
                    if (asprintf(&newickString,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate)<0)
                        return NULL;
                    // snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
                    return newickString;
                }
                else if (p->left == NULL && p->right == NULL)   /* tip of the tree */
                {
                    snprintf(buffer, sizeof(buffer), "tip_i%05d_C%d_%d", p->index,p->indexOldClone,p->indexCurrentClone);
                    strcpy( p->cellName,buffer);
                    //            if (p->isOutgroup == YES)     /* Outgroup */
                    //            {
                    //                //strcat(newickString, "healthycell:%10.9lf", (p->anc1->timePUnits - p->timePUnits) * mutationRate);
                    //                snprintf(newickString,  size,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate);
                    //                return newickString;
                    //            }
                    //else{
                    if  (doUseObservedCellNames == YES)
                    {
                        if (asprintf(&newickString,   "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate)<0)
                            return NULL;
                        //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                        return newickString;
                    }
                    else{
                        if (asprintf(&newickString,   "%s:%10.9lf",  p->cellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate)<0)
                            return NULL;
                        //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                        
                        return newickString;
                    }
                    //  }
                }
                else
                {
                    // fprintf (fpTrees2, "(");
                    if ( p->left != NULL  )
                    {
                        left = toNewickString2 (p->left, mutationRate,   doUseObservedCellNames);
                        //right = toNewickString2 (p->right, mutationRate,   cellNames);
                        //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                        //                free(left);
                        //                left = NULL;
                        //                free(right);
                        //                right=NULL;
                        //                return newickString;
                    }
                    if ( p->right != NULL  )
                    {
                        right = toNewickString2 (p->right, mutationRate,   doUseObservedCellNames);
                        //                snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                        //                free(left);
                        //                left = NULL;
                        //                free(right);
                        //                right=NULL;
                        //                return newickString;
                    }
                    outgroup =toNewickString2 (p->outgroup, mutationRate,   doUseObservedCellNames);
                    if(left!= NULL && right!= NULL && p->anc1 != NULL)
                    {
                        snprintf(buffer, sizeof(buffer), "int_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                        strcpy( p->cellName,buffer);
                        
                        if (asprintf(&newickString, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate )<0)
                            return NULL;
                        //snprintf(newickString, size, "(%s,%s):%10.9lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate );
                        free(left);
                        left = NULL;
                        free(right);
                        right=NULL;
                        free(outgroup);
                        outgroup=NULL;
                    }
                    else if (left != NULL &&  right!= NULL  && p->anc1 == NULL)
                    {
                        snprintf(buffer, sizeof(buffer), "root_i%05d_C%d_%d",  p->index,p->indexOldClone,p->indexCurrentClone);
                        strcpy( p->cellName,buffer);
                        // left = toNewickString2 (p->left, mutationRate,   cellNames);
                        //right = toNewickString2 (p->right, mutationRate,   cellNames);
                        if (asprintf(&newickString,  "(%s,%s);", left, right)<0)
                            return NULL;
                        //snprintf(newickString, size, "(%s,%s);", left, outgroup);
                        free(left);
                        left = NULL;
                        free(right);
                        right=NULL;
                        free(outgroup);
                        outgroup=NULL;
                        // return newickString;
                    }
                    return newickString;
                }
            }
            return newickString;
        }
        /********************* PrintTimes **********************/
        /* Prints to timesfile a detailed description of
         the tree: nodes, times, branch lengths */
        
        void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, TreeNode *nodes,  int thereisOutgroup)
        {
            /* there isn't recombination */
            fprintf (fpTimes, "\n\nDataset %d", replicate + 1);
            fprintf (fpTimes, "\n              ------------ Nodes -------------");
            fprintf (fpTimes, "\n    class    | label  index  (left right anc) |         time     time length    branch length");
            fprintf (fpTimes, "\n----------------------------------------------------------------------------------------------\n");
            ListTimes (0, mutationRate, nodes, fpTimes, thereisOutgroup);
        }
        
        
        /********************* PrintTimes2 **********************/
        /* Prints to timesfile a detailed description of
         the tree: nodes, times, branch lengths */
        
        void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  TreeNode *nodes,  int thereisOutgroup)
        {
            /* there isn't recombination */
            fprintf (fpTimes2, "\n\nDataset %d", replicate + 1);
            fprintf (fpTimes2, "\n              ------------ Nodes -------------");
            fprintf (fpTimes2, "\n    class    | label  index  (left right anc) |         time     time length    branch length");
            fprintf (fpTimes2, "\n----------------------------------------------------------------------------------------------\n");
            ListTimes2 (0, mutationRate, nodes, fpTimes2, thereisOutgroup);
        }
        
        
        
        /********************** ListTimes ************************/
        /* Writes a given tree description from ListTimes   */
        
        void ListTimes (int j, double mutationRate, TreeNode *nodes, FILE *fpTimes, int thereisOutgroup)
        {
            /* It does not list superfluous nodes */
            TreeNode  *p;
            int     i = 0;
            
            do
            {
                p = nodes + i;
                
                if (p->isOutgroup == YES)     /* Outgroup */
                    fprintf (fpTimes, "%13s   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                             "outgroup", Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, p->anc1->timePUnits - p->timePUnits, (p->anc1->timePUnits - p->timePUnits) * mutationRate);
                
                else if (p->anc1 != NULL && p->left != NULL && p->right != NULL)        /* No MRCA, no tip (internal ancester) */
                    fprintf (fpTimes, "%5s_C%dR%d(f)   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                             "int", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, p->anc1->timePUnits - p->timePUnits, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                
                else if (p->anc1 != NULL && p->left == NULL && p->right == NULL)        /* tip */
                    fprintf (fpTimes, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                             "tip", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, p->anc1->timePUnits - p->timePUnits, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                
                else if (p->nodeClass == 5 || (p->anc1 == NULL && p->left != NULL && p->right != NULL))       /* root, MRCA */
                    fprintf (fpTimes, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.4lf      %10.4lf       %10.9lf\n",
                             "root", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->timePUnits, 0.0, 0.0);
                
                else
                    fprintf (fpTimes, "");
                
                i++;
                
                if (i > 2000)
                    exit(-1);
                
            } while    ((thereisOutgroup == NO  && p->anc1  != NULL)    /* no MRCA */
                        ||  (thereisOutgroup == NO  && p->left == NULL)   /* tip */
                        ||  (thereisOutgroup == YES && p->isOutgroup == NO));
        }
        
        
        
        /********************** ListTimes2 ************************/
        /* Writes a given tree description from ListTimes   */
        
        void ListTimes2 (int j,  double mutationRate, TreeNode *nodes,  FILE *fpTimes2, int thereisOutgroup)
        {
            /* It does not list superfluous nodes */
            TreeNode  *p;
            int     i = 0;
            
            do
            {
                p = nodes + i;
                if (p->isOutgroup == YES)     /* Outgroup */
                    fprintf (fpTimes2, "%13s   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                             "outgroup", Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, p->anc1->time - p->time, (p->anc1->time - p->time) * mutationRate);
                else if (p->anc1 != NULL && p->left != NULL && p->right != NULL)        /* No MRCA, no tip (internal ancester) */
                    fprintf (fpTimes2, "%5s_C%dR%d(f)   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                             "int", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, p->anc1->time - p->time, (p->anc1->time - p->time)*mutationRate);
                else if (p->anc1 != NULL && p->left == NULL && p->right == NULL)        /* tip */
                    fprintf (fpTimes2, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                             "tip", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, p->anc1->time - p->time, (p->anc1->time - p->time)*mutationRate);
                else if (p->nodeClass == 5 || (p->anc1 == NULL && p->left != NULL && p->right != NULL))       /* root, MRCA */
                    fprintf (fpTimes2, "%8s_C%dR%d   %4d   %4d  (%4d %4d %4d) |   %10.9lf      %10.9lf       %10.9lf\n",
                             "root", p->indexOldClone, p->indexCurrentClone, Label(p), p->index, Index(p->left), Index(p->right), Index(p->anc1), p->time, 0.0, 0.0);
                else
                    fprintf (fpTimes2, "");
                i++;
                
                if (i > 2000)
                    exit(-1);
                
            } while    ((thereisOutgroup == NO  && p->anc1  != NULL)    /* no MRCA */
                        ||  (thereisOutgroup == NO  && p->left == NULL)   /* tip */
                        ||  (thereisOutgroup == YES && p->isOutgroup == NO));
        }
/************************************************************/
/********************* ProbabilityCloneiFromClonej2 ********************/
/* Obtain the probability that clone i is originated from clone j
 */
double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, Population **populations, int numClones)
{
    double  ProbabilityIJ, AboveTerm, BelowTerm;
    int     l, j;
    double  h, a, b, c, d, e, t, cum;
    Population *p;
    ProbabilityIJ = 0.0;
    AboveTerm = 0.0;
    BelowTerm = 0.0;
    l = 0;
    h = 0.0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    t = 0.0;
    cum = 0.0;
    // calculate h
    // t = CloneTimeOriginInputSTDPOP[i] * (ClonePopSizeMeffectBegin[i] / ClonePopSizeMeffectBegin[j]);
    t = (PopI->timeOriginSTD ) * (PopI->effectPopSize) / ( PopJ->effectPopSize);
    h = CalculateH(t, PopJ->timeOriginSTD, PopJ->delta);
    AboveTerm = ( PopJ->popSize) * h;
    //fprintf (stderr, "AboveTerm = %lf\n", AboveTerm);
    j=0;
    for (l = PopI->order + 1; l < numClones; l++)
    {    p = *(populations + l);
        //fprintf (stderr, "\ni = %d, j = %d, l = %d\n", i, j, l);
        t = (PopI->timeOriginSTD ) * (PopI->effectPopSize) / ( p->effectPopSize);
        h = CalculateH(t, p->timeOriginSTD, p->delta);
        
        cum = cum + ( ( p->popSize) * h);
    }
    
    BelowTerm = cum;
    //fprintf (stderr, "BelowTerm = %lf\n", BelowTerm);
    
    ProbabilityIJ = AboveTerm / BelowTerm;
    //fprintf (stderr, "ProbabilityIJ = %lf\n", ProbabilityIJ);
    
    return ProbabilityIJ;
}
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

void connectNodelets(TreeNode *node )
{
    if (node != NULL)
    {
        if (node->left == NULL && node->right== NULL)
        {
            char * temp;
            node->isLeaf=YES;
            node->nodeLeft= NULL;
            node->nodeRight= NULL;
            node->nodeBack->next = NULL;
            
            node->nodeBack->node_index= node->index;
            if (asprintf(&temp,  "%d_back",  node->label)<0)
                return;
            node->nodeBack->label=temp;
        }
        else
        {
            char * temp1;
            char * temp2;
            char *temp3;
            node->isLeaf=NO;
            node->nodeBack->next=node->nodeLeft;
            node->nodeLeft->next=node->nodeRight;
            node->nodeRight->next =node->nodeBack;
            
            node->nodeLeft->node_index= node->index;
            
            node->nodeRight->node_index= node->index;
            
            node->nodeBack->node_index= node->index;
            if (asprintf(&temp1,  "%d_back",  node->label)<0)
                return;
            node->nodeBack->label=temp1;
            if (asprintf(&temp2,  "%d_left",  node->label)<0)
                return;
            node->nodeLeft->label= temp2;
            if (asprintf(&temp3,  "%d_right",  node->label)<0)
                return;
            node->nodeRight->label= temp3;
            
        }
    }
}

/********************* connectNodes **********************/
/* connectNodes*/
void connectNodes(TreeNode *left, TreeNode *right, TreeNode *ancester  ){
    if (left!=NULL && right!= NULL && ancester!=NULL )
    {
        connectNodelets(left);
        connectNodelets(right);
        connectNodelets(ancester);
        //connect the child nodes
        left->nodeBack->back =ancester->nodeLeft;
        right->nodeBack->back =ancester->nodeRight;
        
        //connect the ancester node
        ancester->nodeLeft->back =left->nodeBack;
        ancester->nodeRight->back =right->nodeBack;
        
        //connect the edges
        left->edgeBack->edge.utree.parent =ancester->nodeLeft;
        right->edgeBack->edge.utree.parent=ancester->nodeRight;
        
        ancester->edgeLeft->edge.utree.child =left->nodeBack;
        ancester->edgeRight->edge.utree.child=right->nodeBack;
        
        ancester->isLeaf=NO;
    }
    else if(left==NULL && right== NULL && ancester!=NULL )
    { // the ancester node is a leaf
        connectNodelets(ancester);
        //connect the child nodes
        ancester->nodeLeft =NULL;
        ancester->nodeRight =NULL;
        ancester->isLeaf=YES;
        
        ancester->edgeLeft=NULL;
        ancester->edgeRight=NULL;
    }
    else if(left!=NULL && right== NULL && ancester!=NULL )
    {
        connectNodelets(left);
        connectNodelets(ancester);
        //connect the child nodes
        left->nodeBack->back =ancester->nodeLeft;
        //connect the ancester node
        ancester->nodeLeft->back =left->nodeBack;
        
        //connect the edges
        left->edgeBack->edge.utree.parent =ancester->nodeLeft;
        
        ancester->edgeLeft->edge.utree.child =left->nodeBack;
        
        ancester->isLeaf=NO;
    }
    else if(left==NULL && right!= NULL && ancester!=NULL )
    {
        connectNodelets(right);
        connectNodelets(ancester);
        
        //connect the child nodes
        right->nodeBack->back =ancester->nodeRight;
        //connect the ancester node
        ancester->nodeRight->back =right->nodeBack;
        //connect the edges
        right->edgeBack->edge.utree.parent=ancester->nodeRight;
        
        ancester->edgeRight->edge.utree.child=right->nodeBack;
        
        ancester->isLeaf=NO;
    }
}
/********************* setLength **********************/
/* setLength*/
void setLength(TreeNode *node )
{
    double lengthEdge;
    if(node->anc1!=NULL)
    {
        //node->nodeBack->length= node->lengthModelUnits;//this takes into account  the mutation rate
        // lengthEdge = node->anc1->timePUnits - node->timePUnits; // this doesnt take into account the mutation rate
        lengthEdge = node->anc1->timePUnits - node->timePUnits;//not in model time, already includes the effect pop size
        node->nodeBack->length= lengthEdge;
        
        if (node->isLeaf==NO)
        {
            //node->nodeLeft->length=  node->timePUnits - node->left->timePUnits;
            //node->nodeRight->length= node->timePUnits - node->right->timePUnits;
            
            node->nodeLeft->length=  node->timePUnits - node->left->timePUnits;
            node->nodeRight->length= node->timePUnits - node->right->timePUnits;
            
            //node->edgeLeft->length  =node->timePUnits - node->left->timePUnits;
            //node->edgeRight->length  =node->timePUnits - node->right->timePUnits;
            node->edgeLeft->length  =node->timePUnits - node->left->timePUnits;
            node->edgeRight->length  =node->timePUnits - node->right->timePUnits;
        }
        
        if(node->edgeBack!=NULL)
        {
            node->edgeBack->length=lengthEdge;
            
        }
        
    }
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
        gammaNumber = -log (RandomUniform(seed));
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
        r = RandomUniform(seed);
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
        r = RandomUniform(seed);
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
        r = RandomUniform(seed);
        g = r-r*r;
        f = (r-0.5)*h/sqrt(g);
        x = b+f;
        if (x <= 0.0)
            continue;
        r = RandomUniform(seed);
        d = 64*r*r*g*g*g;
        if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))
            break;
    }
    return (x);
}

/************************************* SimulateMk2ForSite ***************************************/
/* Simulates the nucleotide substitution process for a given site under Mk2 model (see Lewis 2001)
 with equal rates. 0 is the reference (healthy) allele */

void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed, int doUserTree, int rateVarAmongSites, double altModelMutationRate, SiteStr* allSites, int ***data, int* numMU )
{
    double    probOfChange, uniform, branchLength;
    int     cell, anccell;
    
    if (p != NULL)
    {
        if (p->isOutgroup == NO)
        {
            cell = p->label;
            anccell = p->anc1->label;
            
            if (doUserTree == YES){
                branchLength = p->lengthModelUnits;//>branchLength;
                //                 branchLength = p->length;//>branchLength;
                
            }
            else
            {
                if (rateVarAmongSites == YES)
                    branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
                else
                    branchLength = altModelMutationRate * p->length;
            }
            
            probOfChange = 0.5 - 0.5 * exp (-2.0 * branchLength);
            
            uniform = RandomUniform(seed);
            if (uniform >= probOfChange) /* => no change */
                data[genome][cell][site] = data[genome][anccell][site];
            else /* => there will be change */
            {
                if (data[genome][anccell][site] == 0)
                    data[genome][cell][site] = 1;
                else
                    data[genome][cell][site] = 0;
                
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                numMU=numMU+1;
            }
        }
        SimulateMk2forSite (p->left,  genome, site, seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
        SimulateMk2forSite (p->right, genome, site, seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
    }
}

/************************************* SimulateMk2 **********************************************/
/* Simulates the nucleotide substitution process under the Mk2 model (see Lewis 2001),
 also called Cavender-Farris-Neyman CFN  model or Jukes-Cantor (1969) model for two alleles */

void SimulateMk2 (TreeNode *p, int genome, long int *seed, int* AltModelSites, int  numAltModelSites,int doUserTree, int rateVarAmongSites, double altModelMutationRate, SiteStr* allSites, int ***data, int* numMU)
{
    int     i;
    
    for (i=0; i<numAltModelSites; i++)
        SimulateMk2forSite (p, genome, AltModelSites[i], seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
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
        U = RandomUniform (seed);
    while (U == 0);
    
    exponentialNumber = -log (U) / lambda;
    
    return exponentialNumber;
}

/********************* WhichIUPAC ************************/
/* Returns the IUPAC representation of the genotype */
/*
 UPAC nucleotide code    Base
 A    Adenine
 C    Cytosine
 G    Guanine
 T (or U)    Thymine (or Uracil)
 R    A or G
 Y    C or T
 S    G or C
 W    A or T
 K    G or T
 M    A or C
 B    C or G or T
 D    A or G or T
 H    A or C or T
 V    A or C or G
 N    unknown state
 . or -    gap
 
 This is what we do:
 
 A/A => A
 A/C => M
 A/G => R
 A/T => W
 A/_ => a
 
 C/A => M
 C/C => C
 C/G => S
 C/T => Y
 C/_ => c
 
 G/A => R
 G/C => S
 G/G => G
 G/T => K
 G/_ => g
 
 T/A => W
 T/C => Y
 T/G => K
 T/T => T
 T/_ => t
 
 _/A => a
 _/C => c
 _/G => g
 _/T => t
 _/_ => -
 
 */

char WhichIUPAC (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //AA
            return ('A');
        else if (allele2 == 1)    //AC
            return ('M');
        else if (allele2 == 2)    //AG
            return ('R');
        else if (allele2 == 3)    //AT
            return ('W');
        else if (allele2 == ADO)    //A?
            return ('a');
        else if (allele2 == DELETION)    //A–
            return ('a');
        else
            return ('N');
    }
    else if (allele1 == 1)
    {
        if (allele2 == 0)        //CA
            return ('M');
        else if (allele2 == 1)    //CC
            return ('C');
        else if (allele2 == 2)    //CG
            return ('S');
        else if (allele2 == 3)    //CT
            return ('Y');
        else if (allele2 == ADO)    //C?
            return ('c');
        else if (allele2 == DELETION)    //C–
            return ('c');
        else
            return ('N');
    }
    else if (allele1 == 2)
    {
        if (allele2 == 0)        //GA
            return ('R');
        else if (allele2 == 1)    //GC
            return ('S');
        else if (allele2 == 2)    //GG
            return ('G');
        else if (allele2 == 3)    //GT
            return ('K');
        else if (allele2 == ADO)    //G?
            return ('g');
        else if (allele2 == DELETION)    //G–
            return ('g');
        else
            return ('N');
    }
    else if (allele1 == 3)
    {
        if (allele2 == 0)        //TA
            return ('W');
        else if (allele2 == 1)    //TC
            return ('Y');
        else if (allele2 == 2)    //TG
            return ('K');
        else if (allele2 == 3)    //TT
            return ('T');
        else if (allele2 == ADO)    //T?
            return ('t');
        else if (allele2 == DELETION)    //T–
            return ('t');
        else
            return ('N');
    }
    else if (allele1 == ADO)
    {
        if (allele2 == 0)        //?A
            return ('a');
        else if (allele2 == 1)    //?C
            return ('c');
        else if (allele2 == 2)    //?G
            return ('g');
        else if (allele2 == 3)    //?T
            return ('t');
        else if (allele2 == ADO)    //??
            return ('-');
        else if (allele2 == DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == DELETION)
    {
        if (allele2 == 0)        //-A
            return ('a');
        else if (allele2 == 1)    //-C
            return ('c');
        else if (allele2 == 2)    //-G
            return ('g');
        else if (allele2 == 3)    //-T
            return ('t');
        else if (allele2 == ADO)    //-?
            return ('-');
        else if (allele2 == DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}
/********************* getHealthyTip **********************/
/* getHealthyTip*/
TreeNode *getHealthyTip(TreeNode *treeRootInit)
{
    if (treeRootInit !=NULL && treeRootInit->right!=NULL)
        return treeRootInit->right;
    else
        return NULL;
}

/********************* WhichMut ************************/
/* Returns character representation for binary data */

char WhichMut (int state)
{
    if (state == 0)
        return ('0');
    else if (state == 1)
        return ('1');
    else if (state == ADO)
        return ('?');
    else if (state == DELETION)
        return ('-');
    else
        return ('N');
}
/********************* WhichConsensusBinary ************************/
/* Returns a consensus representation of the binary genotype */
/*
 0/0 => 0
 0/1 => 1
 1/0 => 1
 1/1 => 2
 
 0/_ => 0
 _/0 => 0
 
 1/_ => 2
 _/1 => 2
 
 _/_ => -
 */

char WhichConsensusBinary (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //00
            return ('0');
        else if (allele2 == 1)    //01
            return ('1');
        else if (allele2 == ADO)    //0?
            return ('0');
        else if (allele2 == DELETION)    //0-
            return ('0');
        else
            return ('N');
    }
    else if (allele1 == 1)
    {
        if (allele2 == 0)        //10
            return ('1');
        else if (allele2 == 1)    //11
            return ('2');
        else if (allele2 == ADO)    //1?
            return ('2');
        else if (allele2 == DELETION)    //0-
            return ('2');
        else
            return ('N');
    }
    else if (allele1 == ADO)
    {
        if (allele2 == 0)        //?0
            return ('0');
        else if (allele2 == 1)    //?1
            return ('2');
        else if (allele2 == ADO)    //??
            return ('-');
        else if (allele2 == DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == DELETION)
    {
        if (allele2 == 0)        //-0
            return ('0');
        else if (allele2 == 1)    //-1
            return ('2');
        else if (allele2 == ADO)    //-?
            return ('-');
        else if (allele2 == DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}
/********************** RandomUniformTo ****************************/
/* it returns random uniform in range 0...max-1          */

int RandomUniformTo (int max, long int *seed)
{
    double    rd;
    rd = RandomUniform (seed);
    return (floor(rd*max));
}
/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char WhichNuc (int nucleotide)
{
    if (nucleotide == A)
        return ('A');
    else if (nucleotide == C)
        return ('C');
    else if (nucleotide == G)
        return ('G');
    else if (nucleotide == T)
        return ('T');
    else if (nucleotide == ADO)
        return ('?');
    else if (nucleotide == DELETION)
        return ('-');
    else
        return ('N');
}
/***************************** openFile*******************************/
/* openFile*/
int openFile(FILE **file, char path[MAX_NAME] )
{
    // if ((*fpTrees = fopen(File, "w")) == NULL)
    if ((*file = fopen(path, "w")) == NULL)
    {
        fprintf(stderr, "Can't open %s.\n", path);
        return -1;
    }
    return 0;
}

/********************************** SimulateISMDNAforSite ***********************************/

/*    Simulates a ACGT mutation under an infinite sites model (ISM) for a given site. The branch
 
 where this mutation is placed is chosen according to its length.
 
 The reference (healthy) allele for each site will be determined by the nucleotide frequencies
 
 A mutation matrix will define the probability of changing to other nucleotides given
 
 the healthy alelle chose
 
 */

void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, int ***data, SiteStr* allSites, int  *numMU, double cumMij[4][4],double mutationRate, double *uniform, double *cumBranchLength, double* ran )

{
    
    //    static double    cumBranchLength, uniform, ran;
    
    int             j, cell, anccell, ancstate;
    
    
    
    if (p != NULL)
        
    {
        
        cell = p->label;
        
        
        
        if ( p->anc1 == NULL)
            
        {
            
            *cumBranchLength = 0;
            
            *uniform = RandomUniform(seed) * totalTreeLength;
            
        }
        
        else
            
        {
            
            anccell = p->anc1->label;
            
            if(genome == MATERNAL )
                
                ancstate = p->anc1->maternalSequence[site];
            
            else
                
                ancstate = p->anc1->paternalSequence[site];
            
            //            cumBranchLength += p->length;// ->branchLength;
            
            *cumBranchLength =*cumBranchLength+ p->length;
            
            
            
            if ((*cumBranchLength < *uniform) || /* => there will be no change */
                
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
                
            {
                
                if(genome == MATERNAL )
                    
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                
                else
                    
                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                
                
                
                p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                
                p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                
                p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site];
                
            }
            
            else /* => there will be change */
                
            {
                
                *ran = RandomUniform(seed) * cumMij[ancstate][3];
                
                for (j=0; j<4; j++)
                    
                {
                    
                    if (*ran <= cumMij[ancstate][j])
                        
                    {
                        
                        //data[genome][cell][site] = j;
                        
                        if(genome == MATERNAL )
                            
                            p->maternalSequence[site]=j;
                        
                        else
                            
                            p->paternalSequence[site]=j;
                        
                        break;
                        
                    }
                    
                }
                
                if (genome == MATERNAL){
                    
                    allSites[site].numMutationsMaternal++;
                    
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site]+1;
                    
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                    
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    
                }
                
                else if (genome == PATERNAL){
                    
                    allSites[site].numMutationsPaternal++;
                    
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                    
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site]+1;
                    
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    
                }
                
                allSites[site].numMutations++;
                
                (*numMU)++;
                
            }
            
        }
        
        SimulateISMDNAforSite (p->left, genome, site, doISMhaploid, seed,  totalTreeLength, data, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran);
        
        SimulateISMDNAforSite (p->right, genome, site, doISMhaploid, seed, totalTreeLength, data, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran);
        
    }
    
}

/********************************** SimulateISMForSite ***********************************/
/*    Simulates a 0/1 mutation under an infinite sites model (ISM) for a given site. The branch
 where this mutation is placed is chosen according to its length.
 0 is the reference (healthy) allele
 */
void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, int ***data, SiteStr* allSites, int  *numMU, double cumMij[4][4], double mutationRate, double*    cumBranchLength, double* uniform, int* mutationAdded)
{
    //    static double    cumBranchLength, uniform;
    //    static double    cumBranchLength, uniform;
    int             cell, anccell;
    //    if (*mutationAdded==YES)
    //        return;
    if (p != NULL)
    {
        cell = p->label;
        
        //        if (p->isOutgroup  == YES || p->anc1 == NULL)
        if ( p->anc1 == NULL)
        {
            *cumBranchLength = 0;
            //            double rUniform=RandomUniform(seed) * totalTreeLength;
            *uniform = RandomUniform(seed) * totalTreeLength;
            
        }
        else
        {
            anccell = p->anc1->label;
            //            *cumBranchLength =*cumBranchLength+ p->length;// ->branchLength;
            // *cumBranchLength =*cumBranchLength+ p->lengthModelUnits;// ->branchLength;
            *cumBranchLength =*cumBranchLength+ p->length;
            //            if ((*cumBranchLength < *uniform) || /* => there will be no change */
            //                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))
            //               ||
            //               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
            //              ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0))
            //                )
            if ((*cumBranchLength < *uniform) ||// (*mutationAdded==YES)/* => there will be no change */
                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))  ||
                ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
                ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0)))
            {
                if(genome == MATERNAL )
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                else
                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                
                p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site];
            }
            else /* => there will be change */
            {
                //if (data[genome][anccell][site] == 0)  /* checking all this might be excessive */
                if(genome == MATERNAL && p->anc1->maternalSequence[site]==0 )
                {
                    p->maternalSequence[site]=1;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site]+1;
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    *numMU=*numMU+1;
                    *mutationAdded=YES;
                }
                else if(genome == PATERNAL && p->anc1->paternalSequence[site]==0 )
                {
                    p->paternalSequence[site]=1;
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site]+1;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    *numMU=*numMU+1;
                    *mutationAdded=YES;
                }
                //                if (data[genome][anccell][site] == 0)
                //                {
                //                    data[genome][cell][site] = 1;
                //                  if (genome == MATERNAL)
                //                    { allSites[site].numMutationsMaternal++;
                //                       }
                //                    else // (genome == PATERNAL)
                //                 { allSites[site].numMutationsPaternal++;
                //                        }
                //                    allSites[site].numMutations++;
                //                    *numMU=*numMU+1;
                //                    *mutationAdded=YES;
                //                    return;
                //
                //
                //
                //                }
                else if (genome == MATERNAL && p->anc1->maternalSequence[site]==1)
                {
                    
                    //                    data[genome][cell][site] = 0;
                    fprintf(stderr,"\n\nERROR: anccell=%d, doISMhaploid:%d , site %d in genome %d of cell %d cannot mutate twice under the ISM model",anccell, doISMhaploid ,site, genome, anccell);
                    fprintf (stderr, "\n\nnumber of mutations %d  maternal: %d and paternal:  %d ", p->numbersMutationsUnderSubtreePerSite[site], p->numbersMaternalMutationsPerSite[site], p->numbersPaternalMutationsPerSite[site]);
                    exit(-1);
                }
                else if (genome == PATERNAL && p->anc1->paternalSequence[site]==1)
                {
                    
                    //                    data[genome][cell][site] = 0;
                    fprintf(stderr,"\n\nERROR: anccell=%d, doISMhaploid:%d , site %d in genome %d of cell %d cannot mutate twice under the ISM model",anccell, doISMhaploid ,site, genome, anccell);
                    fprintf (stderr, "\n\nnumber of mutations %d  maternal: %d and paternal:  %d ", p->numbersMutationsUnderSubtreePerSite[site], p->numbersMaternalMutationsPerSite[site], p->numbersPaternalMutationsPerSite[site]);
                    exit(-1);
                }
                else
                {
                    fprintf (stderr, "\n\nERROR: site %d in genome %d of cell %d has an unknow state %d, %d under the ISM model", site, genome, anccell, p->anc1->maternalSequence[site], p->anc1->paternalSequence[site]);
                    exit(-1);
                }
                if (genome == MATERNAL)
                    allSites[site].numMutationsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numMutationsPaternal++;
                allSites[site].numMutations++;
                
                
            }
        }
        //        double goLeftFirst=RandomUniform(seed);
        //        if(goLeftFirst < 0.5)
        //        {
        SimulateISMforSite (p->left, genome, site, doISMhaploid, seed, totalTreeLength, data, allSites, numMU, cumMij, mutationRate, cumBranchLength, uniform, mutationAdded);
        SimulateISMforSite (p->right, genome, site, doISMhaploid, seed,totalTreeLength, data, allSites, numMU,cumMij, mutationRate, cumBranchLength,  uniform, mutationAdded);
        
        //        }
        //        else{
        //            SimulateISMforSite (p->right, genome, site, doISMhaploid, seed, totalTreeLength, data, allSites, numMU, cumMij, mutationRate, cumBranchLength, uniform, mutationAdded);
        //            SimulateISMforSite (p->left, genome, site, doISMhaploid, seed,totalTreeLength, data, allSites, numMU,cumMij, mutationRate, cumBranchLength,  uniform, mutationAdded);
        //
        //
        //        }
    }
}
/************************************* SimulateFiniteDNA **********************************************/
/* Simulates the nucleotide substitution process under a 4-state Markov model including JC, HKY, GTR and GTRnr */
/* Note that beta is set such that mean substitution rate will be 1.
 E.g., for J-C model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4], int numAltModelSites, int *AltModelSites,SiteStr* allSites,  int rateVarAmongSites, double altModelMutationRate, int *numMU, double Root[], double Cijk[])
{
    int     i, j;
    double beta, kappa;
    double Qij[16];
    double mr;
    
    if (doJC == YES)
    {
        beta = 4./3;
    }
    else if (doHKY == YES)
    {
        freqR = freq[A] + freq[G];
        freqY = freq[C] + freq[T];
        freqAG = freq[A] * freq[G];
        freqCT = freq[C] * freq[T];
        kappa = (titv*freqR*freqY)/(freqAG+freqCT);
        beta = 0.5 / (freqR*freqY + kappa*(freqAG+freqCT));
    }
    else if (doGTR == YES || doGTnR == YES)
    {
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                Qij[i*4+j] = Mij[i][j] / Mij[2][3] * freq[j];
        mr=0;
        for (i=0; i<4; i++)
        {
            Qij[i*4+i]=0;
            Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);
            mr-=freq[i]*Qij[i*4+i];
        }
        EigenREV(mr, Qij, Root, Cijk);
    }
    
    for (i=0; i<numAltModelSites; i++)
        SimulateFiniteDNAforSite (p,  genome, AltModelSites[i], allSites,  seed,  rateVarAmongSites,  altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,    beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk);
    
}

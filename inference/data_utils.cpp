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

#include "constants.hpp"
#include "mutationModel.h"
#include "output_functions.hpp"
#include "random.h"
#include "utils.hpp"

/***************************** ReadParametersFromFile *******************************/
/* Reads parameter values from the parameter file */

void ReadParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths,
                            vector<int> &CloneNameBegin,
                            vector<int> &CloneSampleSizeBegin,
                            vector<int> &ClonePopSizeBegin,
                            vector<double> &CloneBirthRateBegin,
                            vector<double> &CloneDeathRateBegin,
                            vector<double> &CloneTimeOriginInput,
                            double Mij[4][4],
                            double freq[4]
                            )

{
    int   j, z;
    char  ch;
    float   argument;
    double    sumPi;
    double argumentDouble;
    double argumentDouble1;
    int argumentInt;
    long int argumentLongInt;
    
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
                    fprintf(stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", (int)argumentInt);
                    PrintUsage();
                }
                programOptions.numDataSets =argumentInt;
                break;
            case '#':
                if (fscanf(stdin, "%lu bytes", &argumentLongInt) != 1 || argumentLongInt < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad seed (#) (%d)\n\n", (int)argumentLongInt);
                    PrintUsage();
                }
                programOptions.userSeed =argumentLongInt;
                break;
                
                
            case 'X':
                if (fscanf(stdin, "%f", &argument) != 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad haplid/diploid chosen (x) (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*Nscaling = (int) argument;
                programOptions.Nscaling =(int) argument;
                if (programOptions.Nscaling < 1 || programOptions.Nscaling > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Haploid/diplod option (x) (1-2) (%d)\n\n", programOptions.Nscaling);
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
                programOptions.numClones=(int) argument;
                
//                *CloneNameBegin =  (int*)malloc( programOptions.numClones * (long) sizeof(int));
//                *CloneSampleSizeBegin = (int*)malloc( programOptions.numClones* (long) sizeof(int));
//                *ClonePopSizeBegin = (int*)calloc( programOptions.numClones , (long) sizeof(int));
//                *CloneBirthRateBegin =  (double*)calloc( programOptions.numClones, (long) sizeof(double));
//                *CloneDeathRateBegin =  (double*)calloc( programOptions.numClones, (long) sizeof(double));
//                *CloneTimeOriginInput =  (double*)calloc( programOptions.numClones, (long) sizeof(double));
                
//                if (*CloneNameBegin == NULL || *CloneSampleSizeBegin == NULL || *ClonePopSizeBegin == NULL || *CloneBirthRateBegin == NULL || *CloneDeathRateBegin == NULL || *CloneTimeOriginInput == NULL)
//                {
//                    fprintf (stderr, "PARAMETER ERROR: Could not allocate variables for clones\n");
//                    exit (1);
//                }
                
                for (j = 0; j <  programOptions.numClones; j++)
                {
                    
                    for (z = 1; z <= NUM_COLS; z++)
                    {
                        if (z == 1)
                        {
                            fscanf(stdin, "%f", &argument);
                            //*(*CloneNameBegin + j) = (int) argument;
                            CloneNameBegin.push_back((int)argument);
                            if ( CloneNameBegin[j]  <= 0 ||  CloneNameBegin[j]  >  programOptions.numClones)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad number for clone %d (should be higher than 0 and lower than the number of clones %d) (%d)\n\n", j,  programOptions.numClones, CloneNameBegin[j] );
                                PrintUsage();
                            }
                        }
                        if (z == 2)
                        {
                            fscanf(stdin, "%f", &argument);
                            CloneSampleSizeBegin.push_back((int) argument);
                            programOptions.numCells= programOptions.numCells + (int) argument;
                            if (CloneSampleSizeBegin[j] < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad sample size for clone %d (should not be negative) (%d)\n\n", j, CloneSampleSizeBegin[j] );
                                PrintUsage();
                            }
                        }
                        if (z == 3)
                        {
                            fscanf(stdin, "%f", &argument);
                            ClonePopSizeBegin.push_back((int)argument);
                            if (ClonePopSizeBegin[j] < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad population size for clone %d (should be higher than 0) (%d)\n\n", j, ClonePopSizeBegin[j] );
                                PrintUsage();
                            }
                        }
                        if (z == 4)
                        {
                            fscanf(stdin, "%f", &argument);
                            CloneBirthRateBegin.push_back((double) argument);
                        }
                        if (z == 5)
                        {
                            fscanf(stdin, "%f", &argument);
                            CloneDeathRateBegin.push_back((double) argument);
                        }
                        if (z == 6)
                        {
                            fscanf(stdin, "%f", &argument);
                            CloneTimeOriginInput.push_back((double) argument);
                            if (CloneTimeOriginInput[j]  < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be negative) (%lf)\n\n", j,  CloneTimeOriginInput[j] );
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
                programOptions.mutationRate=  (double) argumentDouble;
                break;
                
            case 'B':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0  || argument > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphabet (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                //*alphabet = (int) argument;
                programOptions.alphabet =(int) argument;
                break;
                
            case 'Y':
                if (fscanf(stdin, "%d",  &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", (int) argumentInt);
                    PrintUsage();
                }
                programOptions.noisy =(int) argumentInt;
                break;
                
            case 'D':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 ||argumentDouble < 0 ||argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic dropout rate (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions.ADOrate =(double) argumentDouble;
                break;
            case 'O':
                if (fscanf(stdin, "%d", &argumentInt) < 0 )
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    PrintUsage();
                }
                programOptions.outgroupSelection = argumentInt;
                if ( programOptions.outgroupSelection != 0 &&  programOptions.outgroupSelection != 1 &&  programOptions.outgroupSelection != 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n",  programOptions.outgroupSelection);
                    PrintUsage();
                }
                
                if ( programOptions.outgroupSelection == 0)
                {
                    programOptions.thereisOutgroup = NO;
                    programOptions.outgroupBranchLength_RootSample = 0.0;
                    programOptions.outgroupBranchLength_Root1Root2 = 0.0;
                }
                else if (programOptions.outgroupSelection == 1)
                {
                    //*thereisOutgroup = YES;
                    programOptions.thereisOutgroup=YES;
                    programOptions.outgroupBranchLength_Root1Root2 = 0.0;
                    
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root-Sample) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions.outgroupBranchLength_RootSample=argumentDouble;
                }
                else if (programOptions.outgroupSelection == 2)
                {
                    //*thereisOutgroup = YES;
                    programOptions.thereisOutgroup=YES;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root1-Root2) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions.outgroupBranchLength_Root1Root2=argumentDouble;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root2-Sample) value (%f)\n\n", argumentDouble);
                        PrintUsage();
                    }
                    programOptions.outgroupBranchLength_RootSample=argumentDouble;
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
                programOptions.allelicImbalance=argumentDouble;
                break;
            case 'R':
                if (programOptions.doHKY == YES)
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
                programOptions.thereIsMij=YES;
                if (CheckMatrixSymmetry (Mij) == YES)
                {
                    
                    //programOptions.doJC = NO;
                    programOptions.doJC = YES;
                    programOptions.doHKY = NO;
                    //programOptions.doGTR = YES;
                    programOptions.doGTR = NO;
                    programOptions.doGTnR = NO;
                    
                }
                else
                {
                    
                    programOptions.doJC = NO;
                    programOptions.doHKY = NO;
                    //  programOptions.doGTR = NO;
                    programOptions.doGTR =  YES;
                    // programOptions.doGTnR = YES;
                    programOptions.doGTnR = NO;
                }
                break;
            case 'P':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions.propAltModelSites=argumentDouble;
                if (programOptions.propAltModelSites < 0 || programOptions.propAltModelSites > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f). It has to be between 0 and 1\n\n", programOptions.propAltModelSites);
                    PrintUsage();
                }
                if (programOptions.propAltModelSites > 0 && programOptions.altModel != ISMhap && programOptions.doSimulateFixedNumMutations == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a proportion of non-ISM  sites bigger than zero if the number of mutations is fixed\n\n");
                    PrintUsage();
                }
                if (programOptions.alphabet == DNA && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) programOptions.altModel);
                        PrintUsage();
                    }
                }
                else if (programOptions.alphabet == BINARY && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) programOptions.altModel);
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
                programOptions.altModel = (int) argument;
                if (programOptions.alphabet == DNA && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
                        PrintUsage();
                    }
                }
                else if (programOptions.alphabet == BINARY && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == finiteDNA)
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
                    strcpy(filePaths.treeFile, "trees");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths.treeFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths.treeFile[j] = '\0';
                }
                //*doPrintTrees = YES;
                programOptions.doPrintTrees = YES;
                break;
            case 'A':
                if (fscanf(stdin, "%lf %lf %d", &argumentDouble, &argumentDouble1, &argumentInt) != 3)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean/var/model amplification error (%f ; %f ; model=%d)\n\n",argumentDouble, argumentDouble1, argumentInt);
                    PrintUsage();
                }
                programOptions.meanAmplificationError= argumentDouble;
                programOptions.varAmplificationError= argumentDouble1;
                programOptions.simulateOnlyTwoTemplates=argumentInt;
                
                if ( programOptions.meanAmplificationError < 0 ||  programOptions.meanAmplificationError > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean amplification error (%f)\n\n",  programOptions.meanAmplificationError);
                    PrintUsage();
                }
                if ( programOptions.varAmplificationError < 0 || ( programOptions.meanAmplificationError > 0 &&  programOptions.varAmplificationError >= ( programOptions.meanAmplificationError * (1.0 -  programOptions.meanAmplificationError))))
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n",  programOptions.meanAmplificationError);
                    PrintUsage();
                }
                if ( programOptions.simulateOnlyTwoTemplates != 0 &&  programOptions.simulateOnlyTwoTemplates != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 4 templates) or 1 (assume 2 templates)",  programOptions.simulateOnlyTwoTemplates);
                    PrintUsage();
                }
                break;
            case 'E':
                if (fscanf(stdin, "%lf",  &argumentDouble) !=1 ||  argumentDouble < 0 ||  argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing error (%f)\n\n",  argumentDouble);
                    PrintUsage();
                }
                programOptions.sequencingError= argumentDouble;
                break;
//            case 'S':
//                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 1)
//                {
//                    fprintf (stderr, "PARAMETER ERROR: Bad alternative  do simulated data or not (%d)\n\n", (int) argument);
//                    PrintUsage();
//                }
//                programOptions.doSimulateData= (int)argument;
//                break;
            case 'K':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths.timesFile, "times");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths.timesFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths.timesFile[j] = '\0';
                }
                //*doPrintTimes = YES;
                programOptions.doPrintTimes = YES;
                break;
            case 'J':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of mutations (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions.numFixedMutations = (int) argument;
                //programOptions.doSimulateFixedNumMutations = YES;
                programOptions.doSimulateFixedNumMutations = YES;
                if (programOptions.propAltModelSites > 0 && programOptions.altModel != ISMhap)
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
                    programOptions.equalBaseFreq = YES;
                else
                    programOptions.equalBaseFreq = NO;
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
                programOptions.rateVarCoverage = YES;
                break;
            case 'Z':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", argumentDouble);
                    PrintUsage();
                }
                programOptions.SNPrate = YES;
                break;
                
            case 'H':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing coverage (%d)\n\n", (int) argument);
                    PrintUsage();
                }
                programOptions.coverage = (int) argument;
                if (programOptions.coverage > 0){
                    programOptions.doSimulateReadCounts=YES;
                    // *doSimulateReadCounts = YES;
                    
                }
                if (programOptions.genotypingError > 0 && programOptions.doSimulateReadCounts == YES)
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

/***************************** ReadUntil *******************************/
/* Reading in between [] from input files */
void ReadUntil(FILE *fv, char stopChar, string what)
{
    char ch;
    
    ch = fgetc(fv);
    while (!feof(fv) && ch != stopChar)
        ch = fgetc(fv);
    
    if (feof(fv) || ch != stopChar)
    {
        cerr << what << " missing" << endl;
        exit(0);
    }
}

/***************************** InitListClones*******************************/
/* InitListClones*/
void InitListClones(vector<Population *> &populations, int numClones, int verbose, const vector<int> &CloneNameBegin, const vector<int> &CloneSampleSizeBegin, const vector<double> &CloneBirthRateBegin,  const vector<double> &CloneDeathRateBegin, const vector<int> &ClonePopSizeBegin,
                    const vector<double> &CloneTimeOriginInput,  int TotalNumSequences  )
 {
    int z;
    for (z = 0; z <= (numClones - 1); z++)
    {
        int ind = CloneNameBegin[z];
        int ord = 0;
        double timeOriginInput = CloneTimeOriginInput[z];
        int sampleSize = CloneSampleSizeBegin[z];
        int popSize = ClonePopSizeBegin[z];
        double birthRate = CloneBirthRateBegin[z];
        double deathRate = CloneDeathRateBegin[z];
        
        auto pop = new Population(ind, ord, timeOriginInput, sampleSize, popSize, birthRate, deathRate, NO);
        populations.push_back(pop);

        if (verbose > 1)
            printf("\t%d\t\t", CloneNameBegin[z ]);
        if ((z + 1) != CloneNameBegin[z])
        {
            fprintf (stderr, "PARAMETER ERROR: Check order of clones. Clone (%d) in order is different to (%d). (d)\n\n", z, CloneNameBegin[z]);
            PrintUsage();
        }
        if ( pop->delta <= 0)
        {
            fprintf (stderr, "PARAMETER ERROR: The growth rate cannot be lower than the death rate(Delta parameter negative, Delta=(%10.9lf)) for population %d\n\n", pop->delta, z);
            PrintUsage();
        }

        if (verbose > 1)
        {
            printf("\t\t\t%d\t\t\t", CloneSampleSizeBegin[z ]);
            printf("\t%d\t\t\t", ClonePopSizeBegin[z ]);
            printf("\t%lf\t\t",  populations[z]->effectPopSize);
            printf("\t\t\t%lf\t", CloneBirthRateBegin[z ]);
            printf("\t%lf\t", CloneDeathRateBegin[z ]);
            printf("\t%lf\t\t",populations[z]->growthRate);
            //printf("\t%d\t\t", z);
            printf("\n");
        }
    }
}

 void InitNumberNodes(double &TotalBirthRate, double &TotalDeathRate, int &TotalN,  vector<Population *> &populations, ProgramOptions &programOptions) {
    programOptions.TotalNumSequences = 0;
    TotalN = 0;
    TotalBirthRate = 0.0;
    TotalDeathRate = 0.0;
    Population *popI;
    int j;
    for (j = 0; j < programOptions.numClones; j++)
    {
        popI = populations[j];
        popI->FatherPop = NULL;
        programOptions.TotalNumSequences = programOptions.TotalNumSequences + popI->sampleSize;
        TotalN = TotalN + popI->popSize;
        TotalBirthRate = TotalBirthRate + popI->birthRate;
        TotalDeathRate = TotalDeathRate + popI->deathRate;
    }
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    
    programOptions.numCells = programOptions.TotalNumSequences;
}
/***************************** ListClonesAccordingTimeToOrigin*******************************/
/* ListClonesAccordingTimeToOrigin*/
void ListClonesAccordingTimeToOrigin(vector<Population *> &populations, int numClones)
{
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
//    qsort(populations, numClones, sizeof(Population), comparePopulationsByTimeOrigin);
}
/***************************** comparePopulationsByTimeOrigin*******************************/
bool comparePopulationsByTimeOrigin(const void *s1, const void *s2)
{
    Population *p1 = (Population *)s1;
    Population *p2 = (Population *)s2;
    return (p1->timeOriginInput < p2->timeOriginInput);
}
void InitFilesPathsOptions( FilePaths &filePaths, ProgramOptions &programOptions)
{
    strcpy(filePaths.resultsDir, "Results");
    strcpy(filePaths.treeDir, "trees_dir");
    strcpy(filePaths.timesDir, "times_dir");
    strcpy(filePaths.SNVgenotypesDir, "snv_genotypes_dir");
    strcpy(filePaths.SNVhaplotypesDir, "snv_haplotypes_dir");
    strcpy(filePaths.trueHaplotypesDir, "true_haplotypes_dir");
    strcpy(filePaths.fullHaplotypesDir, "full_haplotypes_dir");
    strcpy(filePaths.MLhaplotypesDir, "ML_haplotypes_dir");
    strcpy(filePaths.fullGenotypesDir, "full_genotypes_dir");
    strcpy(filePaths.CATGdir, "catg_dir");
    strcpy(filePaths.VCFdir, "vcf_dir");
    
    strcpy(filePaths.SNVgenotypesFile, "snv_gen");
    strcpy(filePaths.SNVhaplotypesFile, "snv_hap");
    strcpy(filePaths.trueHaplotypesFile, "true_hap");
    strcpy(filePaths.fullHaplotypesFile, "full_hap");
    strcpy(filePaths.MLhaplotypesFile, "ML_hap");
    strcpy(filePaths.fullGenotypesFile, "full_gen");
    strcpy(filePaths.treeFile, "trees");
    strcpy(filePaths.timesFile, "times");
    if (strlen(filePaths.userTreeFile) == 0)
        strcpy(filePaths.userTreeFile, "usertree");
    if (strlen(filePaths.userGenomeFile) == 0)
        strcpy(filePaths.userGenomeFile, "usergenome");
    strcpy(filePaths.CATGfile, "catg");
    strcpy(filePaths.VCFfile, "vcf");
    strcpy(filePaths.logFile, "log");
    strcpy(filePaths.settingsFile, "log");
#ifdef MYDEBUG
    strcpy(filePaths.mutationsFile, "mutations");
#endif
    programOptions.doPrintSNVgenotypes=1;
    programOptions.doPrintSNVhaplotypes=1;
    programOptions.doPrintTrueHaplotypes=1;
    programOptions.doPrintFullGenotypes=1;
    programOptions.doPrintFullHaplotypes=1;
    programOptions.doNGS=1;
    programOptions.doPrintCATG=1;
    if (programOptions.doSimulateData == NO)
    {
        programOptions.doPrintSNVgenotypes = NO;
        programOptions.doPrintSNVhaplotypes = NO;
        programOptions.doPrintTrueHaplotypes = NO;
        programOptions.doPrintFullHaplotypes = NO;
        programOptions.doPrintMLhaplotypes = NO;
        programOptions.doPrintFullGenotypes = NO;
        programOptions.doPrintAncestors = NO;
        programOptions.doSimulateReadCounts = NO;
        programOptions.doPrintCATG = NO;
    }
}
int SimulateData(ProgramOptions &programOptions, vector<int> &CloneNameBegin, vector<int> &CloneSampleSizeBegin,
                 vector<int> &ClonePopSizeBegin,
                 vector<Population *> &populations,
                 FilePaths &filePaths,
                 Files &files,
                 double freq[4],
                 double Mij[4][4])
{
    int i,j,k,z;
    vector<TreeNode *> nodes;
    char *newickString2;
    double totalTreeLength;
    int    HEALTHY_ROOT, TUMOR_ROOT;
    int    numISMdeletions, numISMCNLOH;
    int cumNumMUperTree;
    int     numAltModelSites = 0, numDefaultModelSites = 0, numISMmutations = 0, altModel = 0;

    double cumfreq[4];
    double cumMij[4][4];
    double Eij[4][4];
    double cumEij[4][4];
    char **cellNames;

    double    kappa = 0.0, beta = 0.0, freqR = 0.0, freqY = 0.0, freqAG = 0.0, freqCT = 0.0;
    double    Rmat[6], NRmat[12], Cijk[256], Root[4];
    double *triNucFreq;
    double   cumNumSNVs = 0, cumNumMU = 0, cumNumDEL = 0, cumNumCNLOH = 0, cumCountMLgenotypeErrors = 0;
    double cumNumMUSq = 0, cumNumSNVsSq = 0, cumNumDELSq = 0, cumNumCNLOHSq = 0;
    CellStr     *cell;
    
    
    if ( programOptions.alphabet == DNA)
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
    //SiteStr* allSites = (SiteStr*) calloc (programOptions.numSites, sizeof(SiteStr));
    vector<SiteStr> allSites(programOptions.numSites);
//    if (!allSites)
//    {
//        fprintf (stderr, "Could not allocate the allSites structure\n");
//        exit (-1);
//    }
    for (i=0; i< programOptions.numSites; i++)
    {
        //allSites[i].alternateAlleles = (int *) calloc (4, sizeof(int));
        allSites[i].alternateAlleles = new int[4];
        if (!allSites[i].alternateAlleles)
        {
            fprintf (stderr, "Could not allocate the allSites[i].alternateAlleles[] structure\n");
            exit (-1);
        }
    }
    
    vector<int> SNVsite(programOptions.numSites);
    vector<int> SFS(programOptions.numSites);
    vector<int> variantSites(programOptions.numSites);
    vector<int> DefaultModelSites(programOptions.numSites);
    vector<int> AltModelSites(programOptions.numSites);
    
    vector<TreeNode *> treeTips;
    vector<double> proportionsVector(programOptions.numClones);
    vector<int> varEvent(programOptions.numDataSets);
    vector<double> varTimeGMRCA(programOptions.numDataSets);

    ////////////////////////////////////////////////////////////////////////////////////////////
    TreeNode *init_root = new TreeNode();
    
    double cumNumCA = 0.0, meanNumCA = 0.0, cumNumMIG = 0.0, meanNumMIG = 0.0;
    double numEventsTot = 0.0, countTMRCA = 0.0,TMRCA = 0.0;
    int        numCA = 0, numMIG = 0;
    int dataSetNum;
    long int seedFirst =  programOptions.seed;
    int  numMU = 0, numDEL = 0, numCNLOH = 0, numProposedMU = 0, numSNVs = 0, numFixedMutations = 0, numSNVmaternal = 0;
    programOptions.doUseObservedCellNames=NO;
    //int ***data;
    ValidateParameters(programOptions,CloneNameBegin , CloneSampleSizeBegin, ClonePopSizeBegin);

    InitListPossibleMigrations(populations, programOptions.numClones);
    InitPopulationsCoalescentEvents( programOptions.numClones,  populations);
   
    for (dataSetNum = 0; dataSetNum < programOptions.numDataSets; dataSetNum++)// dataSetNum refers to a simulated tree number
    {
        if (programOptions.doPrintSeparateReplicates == YES)
            PrepareSeparateFiles(0, 1, dataSetNum,  filePaths, programOptions, files);
        /* Reorganize seed per replicate */
        //seed = seedFirst+dataSetNum+10;
        numCA = numMU = numDEL = numProposedMU = TMRCA = numSNVmaternal = 0;
        programOptions.seed = seedFirst + dataSetNum * 1000;
        //fprintf(stderr, "\n seed = %lu \n", seed);
        /* reset variables per replicate */
        if (programOptions.noisy > 0)
        {
            fprintf (stderr, "\rReplicate #%3d/%d", dataSetNum + 1, programOptions.numDataSets);
            fflush (stdout);
        }
        varEvent[dataSetNum] = 0;
        varTimeGMRCA[dataSetNum] = 0.0;
   
        numCA = numMIG = 0;
        countTMRCA = 0.0;
        /* coalescent tree */
        
        if (programOptions.noisy > 1)
            fprintf (stderr, "\n>> Start coalescent tree .. \n");
        
        TreeNode *root = MakeCoalescenceTree2 (&programOptions.seed, populations,
                              programOptions.numNodes,
                              programOptions.numClones,
                              programOptions,
                              cumNumCA,
                              meanNumCA,
                              cumNumMIG,
                              meanNumMIG,
                              numMIG,
                              numCA,
                              numEventsTot,
                              nodes,
                              treeTips,
                              init_root) ;
        if (programOptions.noisy > 1)
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
        newickString2=NULL;
        newickString2 = toNewickString2 ( root, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
        printf("\n newick = %s  \n", newickString2);
      
        if (programOptions.doPrintTrees == YES)
        {
   
            PrintTrees(dataSetNum, root, files.fpTrees, programOptions.mutationRate, programOptions.doUseObservedCellNames);
            PrintTrees2(dataSetNum, root, files.fpTrees2, programOptions.mutationRate, NULL, NO);
        }
        if (programOptions.doPrintTimes == YES)
        {
            PrintTimes(dataSetNum, files.fpTimes, programOptions.mutationRate, nodes, programOptions.thereisOutgroup);
            PrintTimes2(dataSetNum, files.fpTimes2, programOptions.mutationRate, nodes, programOptions.thereisOutgroup);
        }
        if (programOptions.noisy > 1)
        {
            fprintf (stderr, "\nData set %d", dataSetNum + 1);
            fprintf (stderr, "\n\tNumber of coalescence events   =   %d", numCA);
            fprintf (stderr, "\n\tNumber of migration events     =   %d", numMIG);
        }
        // totalTreeLength = SumBranches(treeRootInit[0], programOptions.mutationRate);
        totalTreeLength = SumBranches(root, programOptions.mutationRate);
        cumNumMUperTree=0;
        
        free(newickString2);
        newickString2=NULL;
      
        for (z = 0; z < programOptions.MutationAssignNum; z++)
        {
            numMU=0;
            if (programOptions.doPrintSeparateReplicates == YES)
                PrepareSeparateFilesGenotypes(1, dataSetNum, z,
                                              filePaths, programOptions,files);
            
            //here there was the code
            if (programOptions.doSimulateData == YES)
            {
                numISMmutations = 0;
                numISMdeletions = 0;
                numISMCNLOH = 0;
            }
            for (i=0; i< programOptions.numSites; i++)
            {   allSites[i].numMutations =0;
                allSites[i].numMutationsMaternal =0;
                allSites[i].numMutationsPaternal =0;
            }

            InitializeGenomes (root, &(programOptions.seed), programOptions.alphabet, programOptions.doUserGenome,programOptions.numSites,  allSites, programOptions.doGeneticSignatures,cumfreq, triNucFreq, NULL);
            // SNPrate=0.01;
            // HEALTHY_ROOT=treeRootInit[0]->label;
            //TUMOR_ROOT=treeRootInit[0]->label;
            HEALTHY_ROOT=root->label;
            TUMOR_ROOT=root->label;
            //        if (SNPrate > 0)
            //            AddGermlineVariation (treeRootInit[0], &seed,  numSites, SNPrate, allSites, alphabet,  data,   HEALTHY_ROOT, cumMij );
            


            EvolveSitesOnTree (root, MATERNAL, &(programOptions.seed), programOptions.rateVarAmongSites,  programOptions.numSites,  allSites, programOptions.doGeneticSignatures, programOptions.alphaSites, programOptions.propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , numISMmutations, programOptions.numFixedMutations, numSNVmaternal,  programOptions.doSimulateFixedNumMutations,  programOptions.alphabet,  numMU, cumMij,  programOptions.altModel, programOptions.altModelMutationRate, programOptions.doUserTree,  programOptions.doJC,  programOptions.doHKY,  programOptions.doGTR,
                               programOptions.doGTnR,  freqR,  freqY,
                               freqAG, freqCT, programOptions.titv, freq, Mij ,   Root,  Cijk);
            EvolveSitesOnTree (root, PATERNAL, &(programOptions.seed), programOptions.rateVarAmongSites,  programOptions.numSites,  allSites, programOptions.doGeneticSignatures, programOptions.alphaSites, programOptions.propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , numISMmutations, numFixedMutations, numSNVmaternal,  programOptions.doSimulateFixedNumMutations,  programOptions.alphabet,  numMU, cumMij,  programOptions.altModel, programOptions.altModelMutationRate, programOptions.doUserTree,  programOptions.doJC,  programOptions.doHKY,  programOptions.doGTR,
                               programOptions.doGTnR,  freqR,  freqY,
                               freqAG, freqCT, programOptions.titv, freq, Mij,   Root,  Cijk );
            cumNumMU += numMU;
            cumNumMUSq += pow(numMU,2);
          
            if (programOptions.doPrintTrueHaplotypes == YES)
            {
                if (programOptions.doPrintSeparateReplicates == NO)
                    fprintf (files.fpTrueHaplotypes, "[#%d]\n", z+1);
                //
                PrintTrueFullHaplotypes (files.fpTrueHaplotypes,  nodes, root , programOptions.numNodes, programOptions.doPrintIUPAChaplotypes, programOptions.doPrintAncestors, programOptions.numSites,  programOptions.numCells, programOptions.alphabet, programOptions.doUserTree,    programOptions.doNGS,   NULL, NULL, HEALTHY_ROOT, TUMOR_ROOT, NULL, NO);
            }
            
            if (programOptions.doPrintTrees ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpTrees);
                fclose(files.fpTrees2);
            }
            if (programOptions.doPrintTimes ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpTimes);
                fclose(files.fpTimes2);
            }
            if (programOptions.doPrintTrueHaplotypes ==YES)
            {
                if (programOptions.doPrintSeparateReplicates == YES)
                    fclose(files.fpTrueHaplotypes);
                
            }
       
            
        }/* end of mutation simulation process */
        
//        TreeNode *p;
//        for(i=0; i<programOptions.numNodes; i++)
//        {
//            p= nodes[i] + i;
//            //             free( p->cellName);
//            //    p->cellName=NULL;
//            free( p->maternalSequence);
//            p->maternalSequence=NULL;
//            free( p->paternalSequence);
//            p->paternalSequence=NULL;
//            free(p->numbersMutationsUnderSubtreePerSite);
//            p->numbersMutationsUnderSubtreePerSite=NULL;
//            free(p->numbersMaternalMutationsPerSite);
//            p->numbersMaternalMutationsPerSite=NULL;
//            free(p->numbersPaternalMutationsPerSite);
//            p->numbersPaternalMutationsPerSite=NULL;
//                        free( p->nodeLeft);
//                        p->nodeLeft=NULL;
//                        free(p->nodeRight);
//                        p->nodeRight=NULL;
//                        free(p->nodeBack);
//                        p->nodeBack=NULL;
//                        free(p->edgeBack);
//                        p->edgeBack=NULL;
//                        free(p->edgeLeft);
//                        p->edgeLeft=NULL;
//                        free(p->edgeRight);
//                        p->edgeRight=NULL;
//
//        }
//        free (nodes);
//        nodes=NULL;
   
    }
    
 
 
    for (i=0; i< programOptions.numSites; i++)
    {
        free( allSites[i].alternateAlleles);
    }
    //free(allSites);
    //allSites=NULL;
//    free(SNVsites);
//    SNVsites=NULL;
//    free(SFS);
//    SFS=NULL;
//    free(variantSites);
//    variantSites=NULL;
//    free(DefaultModelSites);
//    DefaultModelSites=NULL;
//    free(AltModelSites);
//    AltModelSites=NULL;
//    free(treeTips);
//    treeTips=NULL;
//    free(proportionsVector);
//    proportionsVector=NULL;
//    free(varEvent);
//    varEvent=NULL;
//    free(varTimeGMRCA);
//    varTimeGMRCA=NULL;
    
    return 0;
}
/***************************** ValidateParameters*******************************/
        /* Validate parameters*/
        
void ValidateParameters(ProgramOptions &programOptions,
                        vector<int> CloneNameBegin , vector<int> CloneSampleSizeBegin, vector<int> ClonePopSizeBegin)
{
            if (programOptions.noisy > 1)
            {
                printf("\n>> Settings ..\n");
                
                printf("Number of replicates = %d\n", programOptions.numDataSets);
                if (programOptions.Nscaling == 1)
                    printf("Haploid data (%d)\n", programOptions.Nscaling);
                else
                    printf("Diploid data (%d)\n", programOptions.Nscaling);
                printf("Number of clones = %d\n", programOptions.numClones);
            }
            for (int j = 0; j < programOptions.numClones; j++)
            {
                // Checking
                //  if (numClones <  CloneNameBegin[j] )
                if (programOptions.numClones <  CloneNameBegin[j] )
                {
                    fprintf (stderr, "PARAMETER ERROR: Clon (%d) is higher than the number of clones (%d). (d)\n\n",  CloneNameBegin[j], programOptions.numClones);
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
void  InitListPossibleMigrations(vector<Population *> &populations, int numClones)
{
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        p->InitListPossibleMigrations(i);
    }
}
/********************** resetMigrationsList ************************/
/* reset list of migrations*/
void resetMigrationsList(vector<Population*> &populations, int numClones){
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        p->resetMigrationsList();
    }
}
/********************** ChooseFatherPopulation************************/
/* choose probabilistically  the father population of a  population  */
Population* ChooseFatherPopulation(vector<Population*> &populations, int numClones, Population  *PopChild,  long int *seed, int noisy) {
    
    Population *pOrigin, *pTarget, *p;
    double pij, ran;
    double *ptr;
    int i, j, k;
    double cumProb[numClones  - (int)(PopChild->order)];
    cumProb[0] = 0.0;
    for (j = PopChild->order + 1; j < numClones; j++)
    {
        cumProb[j - PopChild->order] = 0.0;
        p = populations[j];
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
    Population *result =  populations[PopChild->order + w];
    if (noisy > 3)
        fprintf (stderr, "\nClone %d derived from clone %d\n", PopChild->index, result->index); // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
    /* Update list of migation times considering that clone i comes from clone j */
    if (noisy > 1)
        fprintf (stderr, "\n*** Updating list of migration times (considering that clone %d comes from clone %d) ..\n", PopChild->index, result->index);
    return (result); //the father population has  order  (PopChild->order) + w
}
/********************** AssignCurrentSequencesToPopulation************************/
/* assign current sequences to  population  */
void AssignCurrentSequencesToPopulation(vector<Population *> &populations, vector<TreeNode*> &nodes,
                                        ProgramOptions &programOptions,
                                        int numClones, int numNodes, int noisy,  int TotalNumSequences,
                                        int &numActiveGametes, int &nextAvailable,
                                        int &labelNodes, char* ObservedCellNames[], int doUseObservedCellNames)
{
    Population *pop;
    TreeNode *p;
    int i, j;
    double currentSampleSize;
    numActiveGametes = 0;
    int indexFirstObservedCellName;
    
    vector<int> CumSumNodes(numClones+1);
    //int*  CumSamNodes = (int *) malloc ((numClones + 1)* (long) sizeof(int)); /* cumulative of samples in clones */
    
//            if (!CumSamNodes)
//            {
//                fprintf (stderr, "Could not allocate CumSamNodes (%lu bytes)\n", (numNodes + 1) * (long) sizeof(int));
//                exit (1);
//            }
    
    CumSumNodes[0] = 0;
    
    for (i = 1; i <= numClones; i++)
    {
        pop = populations[i-1];
        pop->resetActiveGametes();
        currentSampleSize = pop->sampleSize;
        CumSumNodes[i] = CumSumNodes[i - 1] + currentSampleSize;
    }
    if (noisy > 1)
        fprintf (stderr, "\n Initial relation nodes-clones:");
    int  cumIndivid = 0;
    int currentPopIndex;
    numActiveGametes = 0;
    
    for (j = 0; j < TotalNumSequences; j++)
    {
        cumIndivid++;
        p = nodes[j];
        //activeGametes[*numActiveGametes] = j;
        p->index = j;
        p->label = j;
        
        labelNodes = j;
        p->nodeClass = 1;
        // fprintf (stderr,"\nIn cumIndivid=%d j=%d", cumIndivid, j);
        for (i = 1; i <= numClones; i++)
        {
            pop = populations[i-1];
            currentPopIndex = pop->index;
            indexFirstObservedCellName= pop->indexFirstObservedCellName;
            // Identify to which clone belongs this node (sample)
            if (cumIndivid <= CumSumNodes[i] && cumIndivid > CumSumNodes[i - 1])
            {
                //fprintf (stderr,"\ncumIndivid=%d <= CumSamNodes[i]=%d, in clone %d\n", cumIndivid, CumSamNodes[i], pop->index);
                pop->idsActiveGametes.push_back(j);
                
                pop->idsGametes.push_back(j);
                pop->numGametes = pop->numGametes+1;

                p->indexOldClone = currentPopIndex;
                p->indexCurrentClone = currentPopIndex;
                p->effectPopSize= pop->effectPopSize;
                p->orderCurrentClone = pop->order;
                p->isLeaf = YES;
                
                if (programOptions.doUseObservedCellNames)
                    strcpy(p->observedCellName, ObservedCellNames[indexFirstObservedCellName + pop->numActiveGametes ]);
                pop->numActiveGametes=pop->numActiveGametes+1;
                break;
            }
        }
        
        if (noisy > 1)
            fprintf (stderr,"\n > The node %d(%d) belongs to clone %d", p->index, cumIndivid, p->indexOldClone);
        numActiveGametes = numActiveGametes + 1;
    }
    CumSumNodes.clear();
    nextAvailable = numActiveGametes;
    labelNodes = labelNodes + 1;
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
void MakeCoalescenceEvent(vector<Population*> &populations, Population *popI, vector<TreeNode *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                          int &labelNodes, double &currentTime, int &numNodes)
        {
            int k, w;
            double rand, ran;
            TreeNode  *p, *q, *r;
            int firstInd, i, j, secondInd=0, newInd=0;
            int choosePairIndividuals = YES;
            
            ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
            
            newInd = nextAvailable;
            //if (noisy > 1)
                fprintf (stderr, "\n Coalescence involving %d and %d to create node %d (in clone %d)", popI->idsActiveGametes[firstInd], popI->idsActiveGametes[secondInd], newInd, popI->index);
            
            
     
            /*  set pointers between nodes */
            p = nodes[popI->idsActiveGametes[firstInd]];
            q = nodes[popI->idsActiveGametes[secondInd]];
            r = nodes[newInd];    /* new ancester */
            r->index = nextAvailable;
            r->label = labelNodes;
            labelNodes=labelNodes+1;
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
            r->time = currentTime;
            r->timePUnits = currentTime * (popI->effectPopSize);
            
            //fprintf (stderr, "\n r->index = %d, r->time = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf, ThisCloneNumber = %d\n", r->index, r->time, ClonePopSizeMeffectBegin[ThisCloneNumber], ThisCloneNumber);
            if (noisy > 1)
                fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
            /* readjust active nodes */
            
            popI->idsActiveGametes[firstInd] = newInd;
            popI->idsActiveGametes[secondInd] = popI->idsActiveGametes[popI->numActiveGametes - 1];;
            numActiveGametes = numActiveGametes - 1; /* less 1 active node */
            
            
            
            
            
            //update list ids nodes
           // popI->idsGametes[popI->numGametes] = newInd;
            popI->numGametes = popI->numGametes +1;
            
            nextAvailable=nextAvailable+1; /* 1 node more is available */
            
            popI->CoalescentEventTimes[ popI->numCompletedCoalescences]=  r->time;
            popI->numActiveGametes = popI->numActiveGametes - 1; /* now this clone
                                                                  has 1 less node */
            
            popI->numCompletedCoalescences= popI->numCompletedCoalescences+1;
            
            fprintf (stderr, "\n pop of order  %d, number of Active gametes %d ", popI->order, popI->numActiveGametes);
            
            for(int i=0; i < popI->idsActiveGametes.size();i++)
                fprintf (stderr, "\n pop of order  %d Active gamete id %d", popI->order, popI->idsActiveGametes[i]);
            /* memory for number of nodes */
            if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
            {
                /* ReallocNodes(&numNodes, activeGametes); */
                if (noisy == 4)
                    fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                numNodes += INCREMENT_NODES;
                /* realloc */
//                *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode) );
//                if (!(*nodes))
//                {
//                    fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
//                    exit (-1);
//                }
//                popI->idsActiveGametes = (int *) realloc (popI->idsActiveGametes, *numNodes * (long) sizeof(int));
//                if (!(popI->idsActiveGametes))
//                {
//                    fprintf (stderr, "Could not reallocate idsActiveGametes for the current population(%lu bytes)\n", *numNodes * (long) sizeof(int));
//                    exit (-1);
//                }
            }
        }
    /********************** BuildTree************************/
    /*  build tree */
    TreeNode *BuildTree(vector<Population* > &populations,
                   Population *CurrentPop,
                   long int *seed,
                   ProgramOptions &programOptions,
                   vector<TreeNode *> &nodes,
                   vector<TreeNode *> &treeTips,
                   TreeNode *tumour_mrca,
                   int &nextAvailable,
                   int &newInd,
                   double &currentTime,
                   int &labelNodes)
    {
        int i, j;
        int indexCurrentTip;
        int  foundSuperflousNode;
        TreeNode *p, *q, *r;
        TreeNode *treeRootInit;
        /********** BUILDING TREES ***********/
        if (programOptions.noisy > 1)
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
            for (i = 0; i < nextAvailable; i++) // available all the nodes
            {
                p = nodes[i];
                
                //fprintf (stderr, "\n\np->index = %d", p->index);
                
                if (p->left == NULL && p->right == NULL && p->anc1 == NULL)
                {
                    // nothing to do with this node because it is not connected to anything
                    //fprintf (stderr, "\n * nothing to do with this node because it is not connected to anything");
                }
                else if (p->left == NULL && p->right == NULL && p->anc1 != NULL)
                {
                    // (*treeTips)[indexCurrentTip]=p;
                    if(indexCurrentTip <  programOptions.TotalNumSequences){
                        //treeTips[indexCurrentTip]=p;
                        treeTips.push_back(p);
                        indexCurrentTip++;
                    }
                    //connectNodes(NULL, NULL, p);
                    // do not do anything with this node because it is a tip
                    //fprintf (stderr, "\n * do not do anything with this node because it is a tip");
                }
                else if (p->left != NULL && p->right == NULL && p->anc1 != NULL )
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
                        
                        //connectNodes(q, r->right, r);
                    }
                    else
                    {
                        r->right = q;
                        q->anc1 = r;
                        p->left = NULL;
                        p->anc1 = NULL;
                        //connectNodes(r->left, q, r);
                    }
                    
                    //fprintf (stderr, "\n - this is a superflous node and can be removed (1)");
                }
                else if (p->left == NULL && p->right != NULL && p->anc1 != NULL )
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
                        //connectNodes(q, r->right, r);
                    }
                    else
                    {
                        r->right = q;
                        q->anc1 = r;
                        p->right = NULL;
                        p->anc1 = NULL;
                        //connectNodes(r->left, q, r);
                    }
                    
                    //fprintf (stderr, "\n - this is a superflous node and can be removed (2)");
                }
                else if (p->left != NULL && p->right != NULL && p->anc1 != NULL)
                {
                    //connectNodes(p->left, p->right, p);
                    // this is an internal node formed by a coalescence event, do not touch
                    //fprintf (stderr, "\n * this is an internal node formed by a coalescence event, do not touch");
                }
                else if (p->left != NULL && p->right != NULL && p->anc1 == NULL)
                {
                    //connectNodes(p->left, p->right, p);
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
                    //setLength(p);
                }
            }
            //fprintf (stderr, "\n");
        }//while
        
        /* about the MRCA */
        newInd=nextAvailable-1;
        p = nodes[newInd];
        //p = *nodes + *newInd; /* because the last one event is the last coalescence */
        //fprintf (stderr, "\n\n\n>> newInd = %d\n", newInd);
        
        if (programOptions.thereisOutgroup == NO)
        {
            //p = *nodes + *newInd;
            p = nodes[newInd];
            p->nodeClass = 5;
            tumour_mrca = p;
            //*treeRootInit = p;
            //treeRootInit[0] = p;
            p->anc1 = NULL;
            treeRootInit = tumour_mrca;
        }
        if (programOptions.thereisOutgroup == YES && programOptions.outgroupSelection > 0)  /*** Root and outgroup ***/
        {
            p = nodes[newInd]; // MRCA
            p->nodeClass = 4;
            
            if (programOptions.noisy > 1)
                fprintf (stderr, "\n\n>> Attaching outgroup .. ");
            
            //fprintf (stderr, "\n>> ThisCloneNumber = %d, ListMigrationTimesInitial[ThisCloneNumber] = %lf, ClonePopSizeMeffectBegin[ThisCloneNumber] = %lf \n", ThisCloneNumber, ListMigrationTimesInitial[ThisCloneNumber], ClonePopSizeMeffectBegin[ThisCloneNumber]);
            
            if (programOptions.outgroupSelection == 1)  /*** Root 2 times and outgroup ***/
                currentTime = CurrentPop->timeOriginSTD; // origin of the clone; // currentTime + (outgroupBranchLength_Root1Root2 / mutationRate); // set time of the new root (from which the MRCA and outgroup nodes are derived)
            else if (programOptions.outgroupSelection == 2) { /*** Root 2 times and outgroup ***/
                currentTime = CurrentPop->timeOriginSTD + (programOptions.outgroupBranchLength_Root1Root2 / CurrentPop->effectPopSize) ; // origin of the clone + time given by the user
                
            }
            else
            {
                fprintf (stderr, "\n\nError simulationg the outgroup. Check input settings\n");
                PrintUsage();
            }
            //TreeNode* healthyRoot = *nodes + *nextAvailable;
            TreeNode* healthyRoot = nodes[nextAvailable];
            healthyRoot->index = nextAvailable;
            healthyRoot->label = labelNodes;
            healthyRoot->effectPopSize= p->effectPopSize;
            labelNodes=labelNodes+1;
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
            healthyRoot->time = currentTime;
            
            int transformingBranchLength=1.001;
            // healthyRoot->time = p->time * transformingBranchLength ;
            healthyRoot->timePUnits = currentTime * healthyRoot->effectPopSize;
            
            p->length = (p->anc1->timePUnits- p->timePUnits);
            //*mutationRate;
            p->lengthModelUnits = (p->anc1->time- p->time);
            //*mutationRate;
            
            healthyRoot->length = 0;
            //        healthyRoot->length = 0;
            
            //        if (noisy > 2)
            //            fprintf (stderr, "DONE");
            //
            nextAvailable++;
            //
            //        /* connect the healthy ancestral cell with the tip healthy cell*/
            //        if (noisy > 2)
            //            fprintf (stderr, "\n>> Adding healthy tip ... ");
            TreeNode* healthyTip = nodes[nextAvailable];
            healthyTip->left = NULL;
            healthyTip->right = NULL;
            healthyTip->effectPopSize= healthyRoot->effectPopSize;
            
            //connectNodes(NULL, NULL, healthyTip);
            
            healthyTip->anc1 = healthyRoot;
            healthyRoot->right = healthyTip;
            healthyTip->time = 0;
            healthyTip->timePUnits = 0;
            double  healthyTipBranchLengthRatio =1;

            healthyTip->length = (healthyTip->anc1->timePUnits- healthyTip->timePUnits);
 
            healthyTip->lengthModelUnits = (healthyTip->anc1->time- healthyTip->time);

            healthyTip->isOutgroup= YES;
            
            //connectNodes(p, healthyTip, healthyRoot);
            //setLength(p);
            //setLength(healthyTip);

            //*treeRootInit=healthyRoot;
            treeRootInit=healthyRoot;
  
        }

        int intLabel = 0;
        if (programOptions.noisy > 1)
            fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
        if (programOptions.thereisOutgroup == YES)
            intLabel = programOptions.TotalNumSequences + 1;
        else
            intLabel = programOptions.TotalNumSequences;
    
        RelabelNodes(treeRootInit, intLabel);
        
        return treeRootInit;
    }
        /********************* PrepareSeparateFiles **********************/
        /* Open individual files to output results */
        
        void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files)
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
            if (programOptions.doPrintTrees == YES)
            {
                if (programOptions.doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.treeDir );
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d", filePaths.resultsDir, filePaths.treeDir, ChainNumber );
                    
                }
                mkdir(File,S_IRWXU);
                
                if (programOptions.doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir,ChainNumber , filePaths.treeFile, paramSetNumber+1, replicate+1);
                }
                //if ((*fpTrees = fopen(File, "w")) == NULL)
                if (openFile(&files.fpTrees, File) == -1)
                {
                    fprintf (stderr, "Can't open \"%s\"\n", File);
                    exit(-1);
                }
                
                sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.treeDir);
                
                if (programOptions.doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_2_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_2_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, ChainNumber, filePaths.treeFile, paramSetNumber+1, replicate+1);
                    
                }
                //sprintf(File,"%s/%s/%s_2_%04d.tre", resultsDir, treeDir, treeFile, replicate+1);
                //if ((*fpTrees2 = fopen(File, "w")) == NULL)
                if (openFile(&files.fpTrees2, File) == -1)
                {
                    fprintf(stderr, "Can't open %s.\n", File);
                    exit(-1);
                }
            }
            
            if (programOptions.doPrintTimes == YES)
            {
                if (programOptions.doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.timesDir);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d", filePaths.resultsDir, filePaths.timesDir, ChainNumber );
                }
                mkdir(File,S_IRWXU);
                if (programOptions.doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, ChainNumber, filePaths.timesFile, paramSetNumber+1, replicate+1);
                    
                }
                //if ((*fpTimes = fopen(File, "w")) == NULL)
                if (openFile(&files.fpTimes, File) == -1)
                {
                    fprintf (stderr, "Can't open \"%s\"\n", File);
                    exit(-1);
                }
                sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.timesDir);
                if (programOptions.doSimulateData ==YES)
                {
                    sprintf(File,"%s/%s/%s_2_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile, paramSetNumber+1, replicate+1);
                }
                else
                {
                    sprintf(File,"%s/%s_Chain%d/%s_2_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, ChainNumber,filePaths.timesFile, paramSetNumber+1, replicate+1);
                    
                }
                // sprintf(File,"%s/%s/%s_2_%04d.txt", resultsDir, timesDir, timesFile, replicate+1);
                
                //if ((*fpTimes2 = fopen(File, "w")) == NULL)
                if (openFile(&files.fpTimes2, File) == -1)
                {
                    fprintf(stderr, "Can't open %s.\n", File);
                    exit(-1);
                }
            }
            
     
        }
        /********************* PrepareSeparateFilesGenotypes **********************/
        /* Open individual genotypes files to output results */
        void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
                                           const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files)
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
            
            if (programOptions.doSimulateData == YES)
            {
                /* contains SNV genotypes for every cell */
                if (programOptions.doPrintSNVgenotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.SNVgenotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVgenotypesDir, filePaths.SNVgenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
                    //if ((*fpSNVgenotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpSNVgenotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains haplotypes for variable sites for every cell */
                if (programOptions.doPrintSNVhaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.SNVhaplotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVhaplotypesDir, filePaths.SNVhaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
                    //if ((*fpSNVhaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpSNVhaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains reference haplotypes (before errors)  for every cell */
                if (programOptions.doPrintTrueHaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.trueHaplotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.trueHaplotypesDir, filePaths.trueHaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    
                    
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
                    //if ((*fpTrueHaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpTrueHaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains ML haplotypes  for every cell */
                if (programOptions.doPrintMLhaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.MLhaplotypesDir);
                    mkdir(File,S_IRWXU);
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.MLhaplotypesDir, filePaths.MLhaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, MLhaplotypesDir, MLhaplotypesFile, replicate+1);
                    // if ((*fpMLhaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpMLhaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                
                /* contains all genotypes (variable or invariable) for every cell */
                if (programOptions.doPrintFullGenotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.fullGenotypesDir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullGenotypesDir, filePaths.fullGenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullGenotypesDir, fullGenotypesFile, replicate+1);
                    //if ((*fpFullGenotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpFullGenotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                /* contains haplotypes for all sites for every cell */
                if (programOptions.doPrintFullHaplotypes == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.fullHaplotypesDir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullHaplotypesDir, filePaths.fullHaplotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, fullHaplotypesDir, fullHaplotypesFile, replicate+1);
                    //if ((*fpFullHaplotypes = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpFullHaplotypes, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                /* contains reads counts and log10 normalized genotype likelihoods for every SNV and cell */
                if (programOptions.doSimulateReadCounts == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.VCFdir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.VCFdir, filePaths.VCFfile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    //sprintf(File,"%s/%s/%s.%04d", resultsDir, VCFdir, VCFfile, replicate+1);
                    //if ((*fpVCF = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpVCF, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
                /* contains reads counts for every SNV and cell */
                if (programOptions.doPrintCATG == YES)
                {
                    sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.CATGdir);
                    mkdir(File,S_IRWXU);
                    
                    sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.CATGdir, filePaths.CATGfile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
                    // sprintf(File,"%s/%s/%s.%04d", resultsDir, CATGdir, CATGfile, replicate+1);
                    //if ((*fpCATG = fopen(File, "w")) == NULL)
                    if (openFile(&files.fpCATG, File) == -1)
                    {
                        fprintf (stderr, "Can't open \"%s\"\n", File);
                        exit(-1);
                    }
                }
            }
        }
/***************************** InitPopulationsCoalescentEvents*******************************/
/* InitPopulationsCoalescentEvents*/
void InitPopulationsCoalescentEvents( int numClones,  vector<Population*> &populations) {
    int i;
    Population *popI;
    for( i = 0 ; i < numClones; i++)
    {
        popI = populations[i];
        popI->InitCoalescentEvents(numClones);
    }
}
/********************************** InitializeGenomes ***********************************/
/* Initialize all genomes with the reference states  */
void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, double cumfreq[4],double *triNucFreq, char **cellNames)
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
                                //p->maternalSequence[site]=p->paternalSequence[site]=i;
                                p->maternalSequence.push_back(i);
                                p->paternalSequence.push_back(i);
                                
                                allSites[site].referenceAllele = p->maternalSequence[site];//   // then allSites[site].referenceAllele hosts the reference genome
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
//                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
//                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                    
                    p->maternalSequence.push_back(p->anc1->maternalSequence[site]);
                    p->paternalSequence.push_back(p->anc1->paternalSequence[site]);
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

TreeNode *MakeCoalescenceTree2 (long int *seed, vector<Population *> &populations,
                           int &numNodes,
                           int numClones,
                           ProgramOptions &programOptions,
                           double      cumNumCA,
                           double meanNumCA,
                           double cumNumMIG,
                           double meanNumMIG,
                           int  &numMIG,
                           int  &numCA,
                           double &numEventsTot,
                           vector<TreeNode *> &nodes, // empty vector
                           vector<TreeNode *> &treeTips,
                           TreeNode    *treeRootInit
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
    
    resetMigrationsList( populations,  numClones);
    
    //allocate memory for the treenodes
//            *nodes = (TreeNode *) malloc ((programOptions.numNodes + 1)* (sizeof(TreeNode)+ 5* programOptions.numSites * sizeof(int) + 2*MAX_NAME * sizeof(char)+ 3 * sizeof(pll_unode_t) + 3*sizeof(pll_tree_edge_t) )); /* nodes */
//            if (!(*nodes))
//            {
//                fprintf (stderr, "Could not allocate nodes (%lu bytes)\n", (programOptions.numNodes+ 1)  * (long) sizeof(TreeNode));
//                exit (1);
//            }

//            for (i=0; i< programOptions.TotalNumSequences; i++){
//                treeTips[i]=NULL;
//            }
    
    //    treeRootInit = (TreeNode **) calloc(1, sizeof(TreeNode *)); /* nodes pointers */
    //    if (!treeRootInit)
    //    {
    //        fprintf (stderr, "Could not allocate treeRootInit (%lu bytes)\n", 1  * (long) sizeof(TreeNode));
    //        exit (1);
    //    }
    /* set everything to null */
    for (i = 0; i < numNodes; i++)
    {
        p = new TreeNode();
        nodes.push_back(p);
    }
     AssignCurrentSequencesToPopulation(populations, nodes, programOptions, numClones, numNodes, programOptions.noisy, programOptions.TotalNumSequences, numActiveGametes,  nextAvailable,
                                       labelNodes, NULL, NO);
    Population *currentPop;
    Population *fatherPop;
    i=0;
    currentTime=0.0;
    while (i < numClones) {
        currentPop = populations[i];
        SimulatePopulation(currentPop, populations,programOptions, seed,
                           programOptions.numNodes,
                           numClones,
                           cumNumCA,
                           meanNumCA,
                           cumNumMIG,
                           meanNumMIG,
                           numMIG,
                           numCA,
                           numEventsTot,
                           nodes,
                           nextAvailable ,
                           numActiveGametes,
                           labelNodes,
                           currentTime,
                           eventNum);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            
            fatherPop= ChooseFatherPopulation(populations, numClones, currentPop, seed,  programOptions.noisy);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants(numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
    //    free (CumSamNodes);
    //   free (activeGametes);
    TreeNode *root = BuildTree(populations,
                                  currentPop,
                                  seed,
                                  programOptions,
                                  nodes,
                                  treeTips,
                                  treeRootInit,
                                  nextAvailable,
                                  newInd,
                                  currentTime,
                                  labelNodes
                                  );
    
    if (programOptions.noisy > 1)
        fprintf (stderr, "\n\n>> Relabeling nodes on tree... \n\n");
    //    if (thereisOutgroup == YES)
    //        intLabel = TotalNumSequences + 2;
    //    else
    //        intLabel = TotalNumSequences;
    
    // RelabelNodes(treeRootInit[0], treeRootInit, &intLabel );
    
    return root;
}
/************************* SimulatePopulation ************************/
/* simulates the evolution of population until its MRCA and can receive inmigrants  ' */ /* */
void SimulatePopulation( Population *popI, vector<Population*> &populations,
                        ProgramOptions &programOptions,
                        long int *seed,
                        int &numNodes,
                        int numClones,
                        double      cumNumCA,
                        double meanNumCA,
                        double cumNumMIG,
                        double meanNumMIG,
                        int  &numMIG,
                        int  &numCA,
                        double &numEventsTot,
                        vector<TreeNode *> &nodes,
                        int &nextAvailable,
                        int &numActiveGametes,
                        int &labelNodes,
                        double &currentTime,
                        int &eventNum)
{
            int  c, d, i, j, w, k, m, cumIndivid, isCoalescence, whichInd,
            firstInd, secondInd, newInd,  foundSuperflousNode,
            isMigration, whichClone, currentNumberAliveClones;
            double     eventTime;
            TreeNode  *p, *q, *r;
            double    ran;
            double    *cumPopulPart;
            eventNum = 0;
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
            Population *incomingPop;
            //    Population *p2;
            //fprintf (stderr, "\n\n> numMigrations= %d \n", numMigrations);
            if (programOptions.noisy > 1)
                fprintf (stderr, "\n\n>> Simulating evolutionary history of clone %d (number active gametes %d, original time to origin %lf)\n", popI->index, popI->numActiveGametes, popI->timeOriginInput);
            if (programOptions.noisy > 1)
                fprintf (stderr, "\n\n> Simulating evolutionary history of clone  or order  %d ..\n", popI->order);
            //fprintf (stderr, "\n\n> Simulating evolutionary history of clone %d ..\n", popI->index);
            currentTime=0;
            while (indexNextMigration < numMigrations) {
                timeNextMigration = popI->immigrantsPopOrderedByModelTime[indexNextMigration].first;
                //fprintf (stderr, "\n\n> numParcialActiveGametes= %d \n", numParcialActiveGametes);
                if ( popI->numActiveGametes >= 2) {
                    ThisRateCA = (double)  popI->numActiveGametes * ((double)  popI->numActiveGametes - 1) / 2.0;
                    ThisTimeCA_W = RandomExponential (ThisRateCA, seed) ;
                    ThisTimeCA_V1 = Population::FmodelTstandard (currentTime, popI->timeOriginSTD, popI->delta);
                    ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
                    // from standard time to model time, GstandardTmodel(V, T, delta)
                    ThisTimeCA_V2 = Population::GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD, popI->delta);
                }
                else
                {
                    ThisTimeCA_V2 = timeNextMigration + 1.0; // it reached a "provisional" MRCA
                    fprintf (stderr, "\n Only 1 active gamete and %d true migrations out of %d true migrations \n", indexNextMigration, numMigrations - 1);
                }
                if ( ThisTimeCA_V2 < timeNextMigration)
                {
                    //choose randomly two lineages to coalesce
                    isCoalescence = YES;
                    isMigration = NO;
                    numCA = numCA + 1;
                    eventNum= eventNum +1;
                    whichClone = popI->index;
                    currentTime = ThisTimeCA_V2; // update current time in model time
                    eventTime = currentTime;
                    
                    if (programOptions.noisy > 1)
                    {
                        fprintf (stderr, "\n\n*** Event %3d *** currentTime (model time) = %lf, currentTime (standard time) = %lf\n", eventNum, ThisTimeCA_V2, ThisTimeCA_V1 );
                        fprintf (stderr, "\n\n*** Event %3d *** currentTime (input units) = %lf\n", eventNum, ThisTimeCA_V2);
                    }
                    if (programOptions.noisy == 4)
                        fprintf (stderr, "* Coalescence *\n");
                    MakeCoalescenceEvent(populations, popI, nodes, numClones, seed, programOptions.noisy, numActiveGametes, nextAvailable, labelNodes, currentTime,  programOptions.numNodes);
                }
                else
                {
                    if (indexNextMigration < numMigrations - 1) //indexNextMigration corresponds to one of the true migrations
                    {
                        isCoalescence = NO;
                        isMigration = YES;
                        numMIG = numMIG + 1;
                        eventNum= eventNum +1;
                        currentTime = timeNextMigration;
                        eventTime = currentTime;
                        //update migration times in model time
                        if (programOptions.noisy > 1)
                        {
                            fprintf (stderr, "\n\n*** Event %3d *** *currentTime (model units) = %lf\n", eventNum, ThisTimeCA_V2);
                        }
                        if (programOptions.noisy == 4)
                        {fprintf (stderr, "* Migration *\n");}
                       // newInd = nextAvailable;
                        
                       // r = nodes[newInd];   /* new ancestor */
                      //  r->index = nextAvailable;
                       // r->label = labelNodes;
                       // labelNodes = labelNodes+1;
                        
                       // r->indexCurrentClone = popI->index;
                       
                       // r->orderCurrentClone = popI->order;
                        //r->nodeClass = 4;
                        
                        //    p = *nodes + MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration]; // root of younger clone
                        //incommingPop = *((popI->immigrantsPopOrderedModelTime) + indexNextMigration );
                        
//                        if (!incommingPop){
//                            fprintf (stderr, "\nError. The incoming population to  poulation %d is empty \n", popI->index);
//                            exit (-1);
//                        }
                        
                        incomingPop = popI->immigrantsPopOrderedByModelTime[indexNextMigration].second;

                        //r->indexOldClone = incommingPop->index;
                        p = nodes[incomingPop->nodeIdAncestorMRCA];
                        
                        printf( "\n The incoming population %d to  population %d with node %d and time %lf", incomingPop->order, popI->order, p->index, p->timePUnits );
                        
                        //p = incomingPop->MRCA;
                        //p = *nodes + (incommingPop->nodeIdAncesterMRCA); // root of younger clone
                        indexNextMigration = indexNextMigration + 1;
                        p->indexCurrentClone = popI->index;
                        p->indexOldClone = incomingPop->index;
                        p->orderCurrentClone = popI->order;
                        // link the nodes
                       // r->left = p;

                       // r->right = NULL;
                        //choosePairIndividuals = NO;
                        
                        //ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals);
                        //q=*nodes + firstInd;
                        //r->right = q;//choose another random living individual of the population
                        
                        //p->anc1 = r;
                        //q->anc1 = r;
                        
                        //connectNodes(p, NULL, r);
                        //p->time = *currentTime;
                        // p->timePUnits = *currentTime * (popI->effectPopSize);
                        
                        //r->time = currentTime;// this is equal to the time of the migration
                        //r->timePUnits = currentTime * (popI->effectPopSize);
                       // nextAvailable=nextAvailable+1; /* 1 node more is available */
                        
                        k = p->indexCurrentClone;
                        incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                        // remove node from old clone in list of active gametes and add the new node of the current clone
                      // popI->idsActiveGametes[popI->numActiveGametes]=r->index;//adding the superfluos node
                      popI->idsActiveGametes[popI->numActiveGametes]=p->index;//adding the superfluos node
                        popI->numActiveGametes = popI->numActiveGametes + 1; /* now this clone has 1 more node */
                        
                        fprintf (stderr, "\n After inmigration. pop of order  %d, number of Active gametes %d ", popI->order, popI->numActiveGametes);
                        
                        for(int i=0; i < popI->idsActiveGametes.size();i++)
                            fprintf (stderr, "\n pop of order  %d Active gamete id %d", popI->order, popI->idsActiveGametes[i]);
                        //                if (noisy > 1)
                        //                    fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, popI->index, incommingPop->nodeIdAncesterMRCA, k);
                        //fprintf (stderr, "Migration, creating node %d (clone %d) derived from node %d (clone %d)", newInd, ThisCloneNumber, MatrixMigrationIDnodeMRCA[ThisCloneNumber][doAmigration], k);
                        if (programOptions.noisy > 1)
                            fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", p->timePUnits);
                        // fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                        /* memory for number of nodes */
                        if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                        {
                            /* ReallocNodes(&numNodes, activeGametes); */
                            if (programOptions.noisy == 4)
                                fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                            numNodes += INCREMENT_NODES;
                            /* realloc */
//                            *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
//                            if (!(*nodes))
//                            {
//                                fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
//                                exit (-1);
//                            }
                            //                    activeGametes = (int *) realloc (activeGametes, *numNodes * (long) sizeof(int));
                            //                    if (!activeGametes)
                            //                    {
                            //                        fprintf (stderr, "Could not reallocate activeGametes (%lu bytes)\n", *numNodes * (long) sizeof(int));
                            //                        exit (-1);
                            //                    }
                        }
                        
                    }
                    else {
                        //origin reached -- create the root
                        currentTime = timeNextMigration;
                        indexNextMigration = indexNextMigration + 1;
                        isCoalescence = NO;
                        isMigration = NO;
                        eventTime = currentTime;
                        if (programOptions.noisy > 1)
                            fprintf (stderr, "\n\n*** Clone origin ***\n");
                        if (programOptions.noisy == 4)
                            fprintf (stderr, "Clone origin %d at time (model units) = %lf\n", popI->index, currentTime);
                        if (popI->order < numClones - 1) // do not do it for the last clone
                        {
                            newInd = nextAvailable;
                            r = nodes[newInd];    /* new ancestor */
                            r->index = nextAvailable;
                            r->label = labelNodes;
                            labelNodes=labelNodes+1;
                            r->indexOldClone =r->indexCurrentClone = popI->index;
                            r->orderCurrentClone = popI->order;
                            r->effectPopSize=popI->effectPopSize;
                            popI->nodeIdAncestorMRCA=newInd;
                            
                            r->nodeClass = 4;
                            
                            firstInd = nextAvailable - 1;
                            //p = nodes + activeGametes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                            p = nodes[firstInd]; // descendant node (previously generated node, nextAvailable - 1)
                            // link the nodes
                            r->left = p;
                            r->right = NULL;
                            p->anc1 = r;
                            r->time = currentTime;
                            r->timePUnits = currentTime * popI->effectPopSize;
                            popI->MRCA = p;
                            
                            //connectNodes(p, NULL, r);
                            //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                            /* readjust active nodes */
                            nextAvailable=nextAvailable+1; /* 1 node more is available */
                            popI->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                            
                            //popI->idsGametes[popI->numGametes] = newInd; r is a superflous node and it will be removed so no need to add it
                            //popI->numGametes = popI->numGametes +1;
                            
                            if (programOptions.noisy > 1)
                                fprintf (stderr, "Creating origin node, it creates node %d derived from node %d", newInd, firstInd);
                            if (programOptions.noisy > 1)
                                fprintf (stderr, "\t|\tCurrentTime (input units) = %lf", r->timePUnits);
                            /* memory for number of nodes */
                            if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                            {
                                /* ReallocNodes(&numNodes); */
                                if (programOptions.noisy == 4)
                                    fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
                                numNodes += INCREMENT_NODES;
                                /* realloc */
//                                *nodes = (TreeNode *) realloc (*nodes, *numNodes  * (long) sizeof(TreeNode));
//                                if (!(*nodes))
//                                {
//                                    fprintf (stderr, "Could not reallocate nodes (%lu bytes)\n", *numNodes  * (long) sizeof(TreeNode));
//                                    exit (-1);
//                                }
//
                            }
                        }
                        else  {//origin of oldest pop reached
                            popI->nodeIdAncestorMRCA=nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                            r = nodes[nextAvailable-1];//popI->idsActiveGametes[0]
                            r->indexOldClone = r->indexCurrentClone = popI->index;
                            r->orderCurrentClone = popI->order;
                            popI->MRCA= r;
                            fprintf (stderr, "\n origin of the oldest population  %d", popI->nodeIdAncestorMRCA);
                        }
                    }
                }
                if (programOptions.noisy > 3)
                {
                    fprintf (stderr, "\nActive nodes (%d):",  popI->numActiveGametes);
                    for (i = 0; i < popI->numActiveGametes; i++)
                        fprintf (stderr, " %d", popI->idsActiveGametes[i]);
                    fprintf (stderr, "\t|\tNext node available = %d", nextAvailable);
                }
            }
            if (programOptions.noisy > 1)
                fprintf (stderr, "\n\nEvolutionary history of clone %d is completed \n", popI->index);
        }
/**************** RelabelNodes **************/
/*  After getting rid of superfluos node, we
 need to relabel those so they are consecutive
 Use the indexes as labels when there
 is recombination */
void RelabelNodes(TreeNode *p, int &intLabel)
{
    if (p != NULL)
    {
        RelabelNodes (p->left, intLabel);
        RelabelNodes (p->right, intLabel);
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
            p->label = intLabel;
            intLabel=intLabel+1;
        }
    }
}
        
        
/******************* toNewickString2 ****************/
/*  Build a string with the Newick representation using the left, right, anc1 pointers  */
char * toNewickString2 ( TreeNode *p, double mutationRate,     int doUseObservedCellNames)
{
    char buffer[1024];
    
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
            
          
            if (asprintf(&newickString,  "healthycell:%10.9lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate)<0)
                return NULL;
         
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
/************************************************************/
/********************* ProbabilityCloneiFromClonej2 ********************/
/* Obtain the probability that clone i is originated from clone j
 */
double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, vector<Population*> &populations, int numClones)
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
    h = Population::CalculateH(t, PopJ->timeOriginSTD, PopJ->delta);
    AboveTerm = ( PopJ->popSize) * h;
    //fprintf (stderr, "AboveTerm = %lf\n", AboveTerm);
    j=0;
    for (l = PopI->order + 1; l < numClones; l++)
    {    p = populations[l];
        //fprintf (stderr, "\ni = %d, j = %d, l = %d\n", i, j, l);
        t = (PopI->timeOriginSTD ) * (PopI->effectPopSize) / ( p->effectPopSize);
        h = Population::CalculateH(t, p->timeOriginSTD, p->delta);
        
        cum = cum + ( ( p->popSize) * h);
    }
    
    BelowTerm = cum;
    //fprintf (stderr, "BelowTerm = %lf\n", BelowTerm);
    
    ProbabilityIJ = AboveTerm / BelowTerm;
    //fprintf (stderr, "ProbabilityIJ = %lf\n", ProbabilityIJ);
    
    return ProbabilityIJ;
}

//void connectNodelets(TreeNode *node )
//{
//    if (node != NULL)
//    {
//        if (node->left == NULL && node->right== NULL)
//        {
//            char * temp;
//            node->isLeaf=YES;
//            node->nodeLeft= NULL;
//            node->nodeRight= NULL;
//            node->nodeBack->next = NULL;
//
//            node->nodeBack->node_index= node->index;
//            if (asprintf(&temp,  "%d_back",  node->label)<0)
//                return;
//            node->nodeBack->label=temp;
//        }
//        else
//        {
//            char * temp1;
//            char * temp2;
//            char *temp3;
//            node->isLeaf=NO;
//            node->nodeBack->next=node->nodeLeft;
//            node->nodeLeft->next=node->nodeRight;
//            node->nodeRight->next =node->nodeBack;
//
//            node->nodeLeft->node_index= node->index;
//
//            node->nodeRight->node_index= node->index;
//
//            node->nodeBack->node_index= node->index;
//            if (asprintf(&temp1,  "%d_back",  node->label)<0)
//                return;
//            node->nodeBack->label=temp1;
//            if (asprintf(&temp2,  "%d_left",  node->label)<0)
//                return;
//            node->nodeLeft->label= temp2;
//            if (asprintf(&temp3,  "%d_right",  node->label)<0)
//                return;
//            node->nodeRight->label= temp3;
//
//        }
//    }
//}

/********************* connectNodes **********************/
/* connectNodes*/
//void connectNodes(TreeNode *left, TreeNode *right, TreeNode *ancester  ){
//    if (left!=NULL && right!= NULL && ancester!=NULL )
//    {
//        connectNodelets(left);
//        connectNodelets(right);
//        connectNodelets(ancester);
//        //connect the child nodes
//        left->nodeBack->back =ancester->nodeLeft;
//        right->nodeBack->back =ancester->nodeRight;
//        
//        //connect the ancester node
//        ancester->nodeLeft->back =left->nodeBack;
//        ancester->nodeRight->back =right->nodeBack;
//        
//        //connect the edges
//        left->edgeBack->edge.utree.parent =ancester->nodeLeft;
//        right->edgeBack->edge.utree.parent=ancester->nodeRight;
//        
//        ancester->edgeLeft->edge.utree.child =left->nodeBack;
//        ancester->edgeRight->edge.utree.child=right->nodeBack;
//        
//        ancester->isLeaf=NO;
//    }
//    else if(left==NULL && right== NULL && ancester!=NULL )
//    { // the ancester node is a leaf
//        connectNodelets(ancester);
//        //connect the child nodes
//        ancester->nodeLeft =NULL;
//        ancester->nodeRight =NULL;
//        ancester->isLeaf=YES;
//        
//        ancester->edgeLeft=NULL;
//        ancester->edgeRight=NULL;
//    }
//    else if(left!=NULL && right== NULL && ancester!=NULL )
//    {
//        connectNodelets(left);
//        connectNodelets(ancester);
//        //connect the child nodes
//        left->nodeBack->back =ancester->nodeLeft;
//        //connect the ancester node
//        ancester->nodeLeft->back =left->nodeBack;
//        
//        //connect the edges
//        left->edgeBack->edge.utree.parent =ancester->nodeLeft;
//        
//        ancester->edgeLeft->edge.utree.child =left->nodeBack;
//        
//        ancester->isLeaf=NO;
//    }
//    else if(left==NULL && right!= NULL && ancester!=NULL )
//    {
//        connectNodelets(right);
//        connectNodelets(ancester);
//        
//        //connect the child nodes
//        right->nodeBack->back =ancester->nodeRight;
//        //connect the ancester node
//        ancester->nodeRight->back =right->nodeBack;
//        //connect the edges
//        right->edgeBack->edge.utree.parent=ancester->nodeRight;
//        
//        ancester->edgeRight->edge.utree.child=right->nodeBack;
//        
//        ancester->isLeaf=NO;
//    }
//}
///********************* setLength **********************/
///* setLength*/
//void setLength(TreeNode *node )
//{
//    double lengthEdge;
//    if(node->anc1!=NULL)
//    {
//        //node->nodeBack->length= node->lengthModelUnits;//this takes into account  the mutation rate
//        // lengthEdge = node->anc1->timePUnits - node->timePUnits; // this doesnt take into account the mutation rate
//        lengthEdge = node->anc1->timePUnits - node->timePUnits;//not in model time, already includes the effect pop size
//        node->nodeBack->length= lengthEdge;
//        
//        if (node->isLeaf==NO)
//        {
//            //node->nodeLeft->length=  node->timePUnits - node->left->timePUnits;
//            //node->nodeRight->length= node->timePUnits - node->right->timePUnits;
//            
//            node->nodeLeft->length=  node->timePUnits - node->left->timePUnits;
//            node->nodeRight->length= node->timePUnits - node->right->timePUnits;
//            
//            //node->edgeLeft->length  =node->timePUnits - node->left->timePUnits;
//            //node->edgeRight->length  =node->timePUnits - node->right->timePUnits;
//            node->edgeLeft->length  =node->timePUnits - node->left->timePUnits;
//            node->edgeRight->length  =node->timePUnits - node->right->timePUnits;
//        }
//        
//        if(node->edgeBack!=NULL)
//        {
//            node->edgeBack->length=lengthEdge;
//            
//        }
//        
//    }
//}

/***************************** Initialize*******************************/
/* Initialize*/
void Initialize( double (*Eij)[4], double (*Mij)[4], double *freq,  ProgramOptions &programOptions ) {
    programOptions.numDataSets = 10;            /* the number of samples to simulate */
    programOptions.numCells = 8;                /* number of cells in each data set */
    programOptions.ploidy = 2;                 /* we assume diploid genomes */
    programOptions.numSites = 10000;                /* number of sites (markers, loci) to simulate = length of the chromosomes */
    //  N = 1000;                    /* effective population size */
    programOptions.numPeriods = 0;                /* number of distinct demographic periods */
    programOptions.doDemographics = NO;        /* whether to implement demographics */
    programOptions.doExponential = YES;            /* whether to do exponential growth */
    programOptions.growthRate = 0;                /* rate for the exponential population growth */
    programOptions.mutationRate = 1.0e-7;        /* nucleotide mutation rate per site per generation */
    programOptions.rateVarAmongLineages = NO;    /* modify rate variation among branches  (to non-clock) */
    programOptions.alphabet = BINARY;          /* alphabet 0/1 or DNA") */
    programOptions.altModel = 0;                /* by default the alternative model will be ISM haploid */
    programOptions.propAltModelSites = 0;        /* proportion of sites that will mutate according to alternative model */
    programOptions.nonISMRelMutRate = 1.0;        /* relative rate alternative/default model for sites */
    programOptions.equalBaseFreq = YES;        /* DNA base frequencies */
    freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
    programOptions.titv = 0.5;                    /* transition/transversion rate ratio */
    programOptions.thereIsMij = NO;            /* mutation rate matrix*/
    Mij[0][0] = Mij[1][1] = Mij[2][2] = Mij[3][3] = 0;  /* mutation probabilities */
    Mij[0][1] = Mij[0][2] = Mij[0][3] = 1.0/3;
    Mij[1][0] = Mij[1][2] = Mij[1][3] = 1.0/3;
    Mij[2][0] = Mij[2][1] = Mij[2][3] = 1.0/3;
    Mij[3][0] = Mij[3][1] = Mij[3][2] = 1.0/3;
    programOptions.thereIsEij = NO;            /* error rate matrix*/
    Eij[0][0] = Eij[1][1] = Eij[2][2] = Eij[3][3] = 0;  /* sequencing error probabilities */
    Eij[0][1] = Eij[0][2] = Eij[0][3] = 1.0/3;
    Eij[1][0] = Eij[1][2] = Eij[1][3] = 1.0/3;
    Eij[2][0] = Eij[2][1] = Eij[2][3] = 1.0/3;
    Eij[3][0] = Eij[3][1] = Eij[3][2] = 1.0/3;
    
    programOptions.doJC = YES;
    programOptions.doHKY = NO;
    programOptions.doGTR = NO;
    programOptions.doGTnR = NO;
    programOptions.rateVarAmongSites = NO;         /* rate variation among different sites along the genome */
    programOptions.alphaSites = infinity;          /* alpha shape of the gamma distribution for rate variation among sites */
    programOptions.alphaBranches = infinity;       /* alpha shape of the gamma distribution for rate variation among lineages */
    programOptions.alphaCoverage = infinity;       /* alpha shape of the gamma distribution for coverage */
    programOptions.doSimulateFixedNumMutations = NO;    /* whether to simulate a fixed number of mutations */
    programOptions.doUserTree = NO;                /* whether to assume a user tree instead od making the coalescent */
    programOptions.doUserGenome = NO;                /* whether to use a user genome instead of a simulated one */
    programOptions.doPrintSNVgenotypes = NO;        /* whether to print SNVs */
    programOptions.doPrintSNVhaplotypes = NO;      /* whether to print haplotypes */
    programOptions.doPrintTrueHaplotypes = NO;      /* whether to print haplotypes without errors */
    programOptions.doPrintMLhaplotypes = NO;          /* whether to print ML haplotypes */
    programOptions.doPrintFullHaplotypes = NO;        /* whether to print sequences */
    programOptions.doPrintFullGenotypes = NO;        /* whether to print all genotypes (variable + invariable) */
    programOptions.doPrintTrees = NO;                /* whether to print the coalescent tree */
    programOptions.doPrintTimes = NO;                /* whether to print coalescent times */
    programOptions.doPrintAncestors = NO;          /* whether to print data for ancestral cells */
    programOptions.doSimulateReadCounts = NO;      /* do not produce reads by default */
    programOptions.doPrintCATG = NO;                /* whether to print read counts for SNVs in CATG format*/
    programOptions.doSimulateData = YES;            /* whether to simulate any data or do inference from real data */
    programOptions.doPrintSeparateReplicates = YES; /* whether to put every replica in its own file */
    
    programOptions.doPrintIUPAChaplotypes = YES;    /* whether to print IUPAC halotypes */
    programOptions.doGeneticSignatures = NO;        /* whether to use a genetic signature to model trinucleotide mutations */
    programOptions.numUserSignatures = 0;            /* by default we do not use a genetic signature */
    programOptions.healthyTipBranchLength = 0;     /* length of the branch leading to the healthy cell */
    programOptions.transformingBranchLength = 0;     /* length of the transforming branch leading to the healthy ancestral cell */
    programOptions.coverage = 0;                    /* NGS  depth for read counts */
    programOptions.rateVarCoverage = NO;            /* there is coverage dispersion */
    programOptions.ADOrate = 0;                    /* allelic dropout */
    programOptions.sequencingError = 0;            /* NGS error rate */
    programOptions.genotypingError = 0;            /* add errors directly in the genotypes */
    programOptions.SNPrate = 0.0;                    /* germline variation rate for rooth healthy genome */
    programOptions.meanAmplificationError = 0;     /* mean of beta distribution for WGA errors */
    programOptions.varAmplificationError = 0;      /* variance of beta distribution for WGA errors */
    programOptions.simulateOnlyTwoTemplates = NO;    /* whether simualate maximum of two templates after single-cell amplification, or there can be all four */
    programOptions.haploidCoverageReduction = 0.5; /* proportion of reads produced when a single allele is present */
    programOptions.allelicImbalance = 0.5;            /* proportion of maternal/ paternal reads */
    programOptions.doubletRate = 0.0;                /* no doublets by default */
    programOptions.numNodes = 3000;                /* initial number of nodes allocated to build the coalescent trees */
    programOptions.seed = time(NULL);                 /* seed for random numbers */
    programOptions.userSeed = 0;                    /* seed entered by the user */
    programOptions.noisy = 1;                        /* level of information to be printed in the screen (see below) */
    programOptions.doNGS = NO;
    programOptions.MutationAssignNum =1;
    
    
}

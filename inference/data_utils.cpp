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
#include <ctime>
#include <chrono>
#include <ctype.h>
#include <sys/timeb.h>
#include <sys/stat.h>
#include <float.h>


#include <boost/filesystem.hpp>
#include <zlib.h>
#include <unistd.h>

extern "C"
    {
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
#include <libpll/pll.h>
    
    }
#include "constants.hpp"
#include "mutationModel.h"
#include "output_functions.hpp"
#include "random.h"
#include "utils.hpp"
#include "treeLikelihood.hpp"
#include "genotype_error_model.hpp"
#include <algorithm>



/***************************** ReadParametersFromFile *******************************/
/* Reads parameter values from the parameter file */

void ReadParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths,
                            std::vector<int> &CloneNameBegin,
                            std::vector<int> &CloneSampleSizeBegin,
                            std::vector<int> &ClonePopSizeBegin,
                            std::vector<double> &CloneBirthRateBegin,
                            std::vector<double> &CloneDeathRateBegin,
                            std::vector<double> &CloneTimeOriginInput,
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
    long double argumentLong;
    int argumentInt;
    long int argumentLongInt;
    
    /* Used: N # X C U B Y D O I R P M T A E K J F ? V Z H S */
    
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
        //ch = toupper(ch);
        switch (ch)
        {
                
            case 'n':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                programOptions.numDataSets =argumentInt;
                break;
            case 'q':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt > 1  || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad option for simulate paramaters from priors (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                programOptions.doSimulateFromPriors =argumentInt;
                break;
            case '#':
                if (fscanf(stdin, "%lu bytes", &argumentLongInt) != 1 || argumentLongInt < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad seed (#) (%d)\n\n", (int)argumentLongInt);
                    Output::PrintUsage();
                }
                programOptions.userSeed =argumentLongInt;
                break;
                
                
            case 'x':
                if (fscanf(stdin, "%f", &argument) != 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad haplid/diploid chosen (x) (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                //*Nscaling = (int) argument;
                programOptions.Nscaling =(int) argument;
                if (programOptions.Nscaling < 1 || programOptions.Nscaling > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Haploid/diplod option (x) (1-2) (%d)\n\n", programOptions.Nscaling);
                    Output::PrintUsage();
                }
                break;
                
            case 'c':
                if (fscanf(stdin, "%f", &argument) != 1 || argument <= 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of clones (must be 1 or higher) (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                //*numClones = (int) argument;
                programOptions.numClones=(int) argument;
                
                
                
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
                                Output::PrintUsage();
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
                                Output::PrintUsage();
                            }
                        }
                        if (z == 3)
                        {
                            fscanf(stdin, "%f", &argument);
                            ClonePopSizeBegin.push_back((int)argument);
                            if (ClonePopSizeBegin[j] < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad population size for clone %d (should be higher than 0) (%d)\n\n", j, ClonePopSizeBegin[j] );
                                Output::PrintUsage();
                            }
                        }
                        if (z == 4)
                        {
                            fscanf(stdin, "%f", &argument);
                            CloneBirthRateBegin.push_back((double) argument);
                        }
                        if (z == 5)
                        {
                            fscanf(stdin, "%Lf", &argumentLong);
                            CloneDeathRateBegin.push_back((double) argumentLong);
                        }
                        if (z == 6)
                        {
                            fscanf(stdin, "%f", &argument);
                            CloneTimeOriginInput.push_back((double) argument);
                            if (CloneTimeOriginInput[j]  < 0)
                            {
                                fprintf (stderr, "PARAMETER ERROR: Bad time to origin for clone %d (should not be negative) (%lf)\n\n", j,  CloneTimeOriginInput[j] );
                                Output::PrintUsage();
                            }
                            if (CloneTimeOriginInput[j]==0)
                                programOptions.doEstimateTimesOriginClones = YES;
                        }
                        
                    }
                    
                }
                break;
                
                
            case 'u':
                
                if (fscanf(stdin, "%lf", &argumentDouble) != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mutation rate (%f) \n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.mutationRate=  (double) argumentDouble;
                break;
            case 'b':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0  || argument > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alphabet (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                //*alphabet = (int) argument;
                programOptions.alphabet =(int) argument;
                break;
                
            case 'y':
                if (fscanf(stdin, "%d",  &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", (int) argumentInt);
                    Output::PrintUsage();
                }
                programOptions.noisy =(int) argumentInt;
                break;
                
            case 'D':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 ||argumentDouble < 0 ||argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic dropout rate (%f)\n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.fixedADOrate =(double) argumentDouble;
                break;
            case 'o':
                if (fscanf(stdin, "%d", &argumentInt) < 0 )
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    Output::PrintUsage();
                }
                programOptions.outgroupSelection = argumentInt;
                if ( programOptions.outgroupSelection != 0 &&  programOptions.outgroupSelection != 1 &&  programOptions.outgroupSelection != 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n",  programOptions.outgroupSelection);
                    Output::PrintUsage();
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
                        Output::PrintUsage();
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
                        Output::PrintUsage();
                    }
                    programOptions.outgroupBranchLength_Root1Root2=argumentDouble;
                    if (fscanf(stdin, "%lf", &argumentDouble) < 0)
                    {
                        fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch length (Root2-Sample) value (%f)\n\n", argumentDouble);
                        Output::PrintUsage();
                    }
                    programOptions.outgroupBranchLength_RootSample=argumentDouble;
                }
                else
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad selection for outgroup %d (0: No outgroup, 1 outgroup with one branch, 2 outgroup with two branches\n\n", argumentInt);
                    Output::PrintUsage();
                }
                //outgroupBranchLength_RootSample = 0;
                break;
            case 'I':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad allelic imbalance (%f)\n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.allelicImbalance=argumentDouble;
                break;
            case 'r':
                if (programOptions.doHKY == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: Cannot specify a mutation matrix (GTR model) and a ti/tv (HKY model) at the same time\n\n");
                    Output::PrintUsage();
                }
                if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf",
                           &Mij[0][0], &Mij[0][1], &Mij[0][2], &Mij[0][3],
                           &Mij[1][0], &Mij[1][1], &Mij[1][2], &Mij[1][3],
                           &Mij[2][0], &Mij[2][1], &Mij[2][2], &Mij[2][3],
                           &Mij[3][0], &Mij[3][1], &Mij[3][2], &Mij[3][3])!=16)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix (-rx x x x x x x x x x x) (AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)\n\n");
                    Output::PrintUsage();
                }
                
                if (Mij[0][0] != 0  || Mij[1][1] != 0 || Mij[2][2] != 0 || Mij[3][3] != 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix: diagonals should be 0 \n\n");
                    Output::PrintUsage();
                }
                //*thereIsMij = YES;
                programOptions.thereIsMij=YES;
                if (Utils::CheckMatrixSymmetry (Mij) == YES)
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
            case 'R':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr,  "\n PARAMETER ERROR: Bad haploid coverage reduction (%f)\n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.haploidCoverageReduction = argumentDouble;
                
            case 'p':
                if (fscanf(stdin, "%lf", &argumentDouble) !=1 || argumentDouble < 0 || argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f)\n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.propAltModelSites=argumentDouble;
                if (programOptions.propAltModelSites < 0 || programOptions.propAltModelSites > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad proportion of alternative model sites (%f). It has to be between 0 and 1\n\n", programOptions.propAltModelSites);
                    Output::PrintUsage();
                }
                if (programOptions.propAltModelSites > 0 && programOptions.altModel != ISMhap && programOptions.doSimulateFixedNumMutations == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a proportion of non-ISM  sites bigger than zero if the number of mutations is fixed\n\n");
                    Output::PrintUsage();
                }
                if (programOptions.alphabet == DNA && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the alt model (%d) specified are incompatible", (int) programOptions.altModel);
                        Output::PrintUsage();
                    }
                }
                else if (programOptions.alphabet == BINARY && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the alt model (%d) specified are incompatible", (int) programOptions.altModel);
                        Output::PrintUsage();
                    }
                }
                break;
            case 'm':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 || argument > 2)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad alternative mutation model (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                programOptions.altModel = (int) argument;
                if (programOptions.alphabet == DNA && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == Mk)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The DNA alphabet and the model (%d) specified are incompatible", (int) argument);
                        Output::PrintUsage();
                    }
                }
                else if (programOptions.alphabet == BINARY && programOptions.propAltModelSites > 0)
                {
                    if (programOptions.altModel == finiteDNA)
                    {
                        fprintf (stderr, "PARAMETER ERROR: The binary alphabet and the model (%d) specified are incompatible", (int) argument);
                        Output::PrintUsage();
                    }
                }
                break;
            case 't':
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
                    Output::PrintUsage();
                }
                programOptions.meanAmplificationError= argumentDouble;
                programOptions.varAmplificationError= argumentDouble1;
                programOptions.simulateOnlyTwoTemplates=argumentInt;
                
                if ( programOptions.meanAmplificationError < 0 ||  programOptions.meanAmplificationError > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean amplification error (%f)\n\n",  programOptions.meanAmplificationError);
                    Output::PrintUsage();
                }
                if ( programOptions.varAmplificationError < 0 || ( programOptions.meanAmplificationError > 0 &&  programOptions.varAmplificationError >= ( programOptions.meanAmplificationError * (1.0 -  programOptions.meanAmplificationError))))
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n",  programOptions.meanAmplificationError);
                    Output::PrintUsage();
                }
                if ( programOptions.simulateOnlyTwoTemplates != 0 &&  programOptions.simulateOnlyTwoTemplates != 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad simulateOnlyTwoTemplates error (%d); it has to be 0 (assume 4 templates) or 1 (assume 2 templates)",  programOptions.simulateOnlyTwoTemplates);
                    Output::PrintUsage();
                }
                break;
            case 'E':
                if (fscanf(stdin, "%lf",  &argumentDouble) !=1 ||  argumentDouble < 0 ||  argumentDouble > 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing error (%f)\n\n",  argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.sequencingError= argumentDouble;
                programOptions.doNGS = NO;
                break;
            case 'G':
                if (fscanf(stdin, "%lf %lf", &argumentDouble, &argumentDouble1) != 2)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean/var genotyping error (%f ; %f)\n\n",argumentDouble, argumentDouble1);
                    Output::PrintUsage();
                }
                
                programOptions.meanGenotypingError= argumentDouble;
                programOptions.varGenotypingError= argumentDouble1;
                if ( programOptions.meanGenotypingError < 0 ||  programOptions.meanGenotypingError > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean genotyping error (%f)\n\n",  programOptions.meanGenotypingError);
                    Output::PrintUsage();
                }
                if ( programOptions.varGenotypingError < 0 ||  programOptions.varGenotypingError > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad var genotyping error (%f)\n\n",  programOptions.varGenotypingError);
                    Output::PrintUsage();
                }
                break;
            case 'X':
                if (fscanf(stdin, "%lf %lf", &argumentDouble, &argumentDouble1) != 2)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean/var ADO error (%f ; %f)\n\n",argumentDouble, argumentDouble1);
                    Output::PrintUsage();
                }
                programOptions.meanADOcell= argumentDouble;
                programOptions.varADOcell= argumentDouble1;
                
                if ( programOptions.meanADOcell < 0 ||  programOptions.meanADOcell > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean genotyping error (%f)\n\n",  programOptions.meanADOcell);
                    Output::PrintUsage();
                }
                if ( programOptions.varADOcell < 0 || ( programOptions.varADOcell > 0 &&  programOptions.varAmplificationError >= ( programOptions.meanADOcell * (1.0 -  programOptions.meanADOcell))))
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n",  programOptions.meanADOcell);
                    Output::PrintUsage();
                }
                
                break;
            case 'F':
                if (fscanf(stdin, "%lf %lf", &argumentDouble, &argumentDouble1) != 2)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean/var ADO error (%f ; %f)\n\n",argumentDouble, argumentDouble1);
                    Output::PrintUsage();
                }
                programOptions.meanADOsite= argumentDouble;
                programOptions.varADOsite= argumentDouble1;
                
                if ( programOptions.meanADOsite < 0 ||  programOptions.meanADOsite > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad mean genotyping error (%f)\n\n",  programOptions.meanADOsite);
                    Output::PrintUsage();
                }
                if ( programOptions.varADOsite < 0 || ( programOptions.varADOsite > 0 &&  programOptions.varADOsite >= ( programOptions.meanADOsite * (1.0 -  programOptions.meanADOsite))))
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad variance amplification error (%f); it has to be < mean*(1-mean)\n\n",  programOptions.meanADOsite);
                    Output::PrintUsage();
                }
                
                break;
            case 's':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 0 )
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of sites (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                programOptions.numSites= (int)argument;
                break;
            case 'k':
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
            case 'j':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad number of mutations (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                programOptions.numFixedMutations = (int) argument;
                //programOptions.doSimulateFixedNumMutations = YES;
                programOptions.doSimulateFixedNumMutations = YES;
                if (programOptions.propAltModelSites > 0 && programOptions.altModel != ISMhap)
                {
                    fprintf (stderr, "PARAMETER ERROR: You cannot specify a fixed number of mutations if there is any non-ISM  site. Set the proportion of non-ISM diploid sites to zero\n\n");
                    Output::PrintUsage();
                }
                break;
            case 'f':
                if (fscanf(stdin, "%lf %lf %lf %lf", &freq[0], &freq[1], &freq[2], &freq[3])!=4)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad Base Frequencies\n\n");
                    Output::PrintUsage();
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
                Output::PrintUsage();
                break;
                /*case 'H':
                 Output::PrintUsage();
                 break;*/
            case 'V':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble <= 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad coverage dispersion (%f)\n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.rateVarCoverage = YES;
                break;
            case 'z':
                if (fscanf(stdin, "%lf", &argumentDouble)!=1 || argumentDouble < 0)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad germline SNP rate (%f)\n\n", argumentDouble);
                    Output::PrintUsage();
                }
                programOptions.SNPrate = YES;
                break;
                
            case 'H':
                if (fscanf(stdin, "%f", &argument) !=1 || argument < 1)
                {
                    fprintf (stderr, "PARAMETER ERROR: Bad sequencing coverage (%d)\n\n", (int) argument);
                    Output::PrintUsage();
                }
                programOptions.coverage = (int) argument;
                if (programOptions.coverage > 0){
                    programOptions.doSimulateReadCounts=YES;
                    // *doSimulateReadCounts = YES;
                    
                }
                if (programOptions.genotypingError > 0 && programOptions.doSimulateReadCounts == YES)
                {
                    fprintf (stderr, "PARAMETER ERROR: Cannot specify a coverage larger than 0, which implies read count generation, and a genotyping error at the same time\n\n");
                    Output::PrintUsage();
                }
                break;
            default :
                fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
                Output::PrintUsage();
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
void ReadUntil(FILE *fv, char stopChar, std::string what)
{
    char ch;
    
    ch = fgetc(fv);
    while (!feof(fv) && ch != stopChar)
        ch = fgetc(fv);
    
    if (feof(fv) || ch != stopChar)
    {
        std::cerr << what << " missing" << std::endl;
        exit(0);
    }
}

/***************************** InitListClones*******************************/
/* InitListClones*/
void InitListClones(std::vector<Population *> &populations, int numClones, int verbose, const std::vector<int> &CloneNameBegin, const std::vector<int> &CloneSampleSizeBegin, const std::vector<double> &CloneBirthRateBegin,  const std::vector<double> &CloneDeathRateBegin, const std::vector<int> &ClonePopSizeBegin,
                    const std::vector<double> &CloneTimeOriginInput,  int &TotalNumSequences,int  doEstimateTimesOriginClones, ProgramOptions &programOptions, std::vector<gsl_rng *> rngGslvector )
{
    int z;
    int totalSampleSize=0;
    int totalPopulationSize=0;
    for (z = 0; z < numClones; z++){
        totalSampleSize+=CloneSampleSizeBegin[z];
        totalPopulationSize+=ClonePopSizeBegin[z];
    }
    TotalNumSequences=totalSampleSize;
    for (z = 0; z <= (numClones - 1); z++)
    {
        int ind = CloneNameBegin[z];
        int ord = 0;
        double timeOriginInput = CloneTimeOriginInput[z];
        int sampleSize = CloneSampleSizeBegin[z];
        int popSize = ClonePopSizeBegin[z];
        double birthRate = CloneBirthRateBegin[z];
        double deathRate = CloneDeathRateBegin[z];
        
        Population* pop;
        if (programOptions.doSimulateFromPriors==0)
            pop = new Population(ind, ord, timeOriginInput, sampleSize, totalSampleSize, popSize, totalPopulationSize, birthRate, deathRate, doEstimateTimesOriginClones);
        else{
            sampleSize= Random::randomUniformIntegerInterval(rngGslvector[1], programOptions.minSampleSize, programOptions.maxSampleSize);
            long double delta = Random::RandomExponential(1, NULL, true, rngGslvector[1], NULL);
            long double theta =  Random::RandomExponential(1, NULL, true, rngGslvector[1], NULL);
            pop = new Population(ind, ord,sampleSize, delta, theta,   programOptions);
            
        }
        
        
        populations.push_back(pop);
        
        if (verbose > 1)
            printf("\t%d\t\t", CloneNameBegin[z ]);
        if ((z + 1) != CloneNameBegin[z])
        {
            fprintf (stderr, "PARAMETER ERROR: Check order of clones. Clone (%d) in order is different to (%d). (d)\n\n", z, CloneNameBegin[z]);
            Output::PrintUsage();
        }
        if ( pop->delta <= 0)
        {
            fprintf (stderr, "PARAMETER ERROR: The growth rate cannot be lower than the death rate(Delta parameter negative, Delta=(%10.9Lf)) for population %d\n\n", pop->delta, z);
            Output::PrintUsage();
        }
        
        if (verbose > 1)
        {
            printf("\t\t\t%d\t\t\t", CloneSampleSizeBegin[z ]);
            printf("\t%d\t\t\t", ClonePopSizeBegin[z ]);
            printf("\t%Lf\t\t",  populations[z]->x);
            printf("\t\t\t%f\t", CloneBirthRateBegin[z ]);
            printf("\t%f\t", CloneDeathRateBegin[z ]);
            printf("\t%Lf\t\t",populations[z]->growthRate);
            //printf("\t%d\t\t", z);
            printf("\n");
        }
    }
}

void InitNumberNodes(double &TotalBirthRate, double &TotalDeathRate, int &TotalN,  std::vector<Population *> &populations, ProgramOptions &programOptions) {
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
        //TotalN = TotalN + popI->popSize;
        TotalBirthRate = TotalBirthRate + popI->birthRate;
        TotalDeathRate = TotalDeathRate + popI->deathRate;
    }
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    
    programOptions.numCells = programOptions.TotalNumSequences;
}
/***************************** ListClonesAccordingTimeToOrigin*******************************/
/* ListClonesAccordingTimeToOrigin*/
void ListClonesAccordingTimeToOrigin(std::vector<Population *> &populations, int numClones)
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
    strcpy(filePaths.likelihoodOuput, "loglik_file");
    strcpy(filePaths.tempInputStan, "stan_dump_file");
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
void InitFiles(Files &files){
    
    files.fplog = new FilePath();
    files.fpTrees = new FilePath();
    files.fpTrees2 = new FilePath();
    files.fpTimes = new FilePath();
    files.fpTimes2 = new FilePath();
    
    files.fpSNVgenotypes= new FilePath();
    files.fpSNVhaplotypes= new FilePath();
    files.fpTrueHaplotypes= new FilePath();
    files.fpFullGenotypes= new FilePath();
    files.fpFullHaplotypes= new FilePath();
    files.fpVCF= new FilePath();
    files.fpCATG= new FilePath();
    files.fpMLhaplotypes= new FilePath();
    files.fplog= new FilePath();
    files.fpTreeOutput= new FilePath();
    files.fpLikelihood = new FilePath();
    files.fpStanDump= new FilePath();
    
    Utils::init_to_empty_str(files.fplog->path);
    Utils::init_to_empty_str(files.fpTrees->path);
    Utils::init_to_empty_str(files.fpTrees2->path);
    Utils::init_to_empty_str(files.fpTimes->path);
    Utils::init_to_empty_str(files.fpTimes2->path);
    
    Utils::init_to_empty_str(files.fpSNVgenotypes->path);
    Utils::init_to_empty_str(files.fpSNVhaplotypes->path);
    Utils::init_to_empty_str(files.fpTrueHaplotypes->path);
    Utils::init_to_empty_str(files.fpFullGenotypes->path);
    Utils::init_to_empty_str(files.fpFullHaplotypes->path);
    
    Utils::init_to_empty_str(files.fpVCF->path);
    Utils::init_to_empty_str(files.fpMLhaplotypes->path);
    Utils::init_to_empty_str(files.fplog->path);
    Utils::init_to_empty_str(files.fpFullGenotypes->path);
    Utils::init_to_empty_str(files.fpTreeOutput->path);
    Utils::init_to_empty_str(files.fpLikelihood->path);
    Utils::init_to_empty_str(files.fpStanDump->path);
    
}
int SimulateData(ProgramOptions &programOptions, std::vector<int> &CloneNameBegin, std::vector<int> &CloneSampleSizeBegin,
                 std::vector<int> &ClonePopSizeBegin,
                 std::vector<Population *> &populations,
                 FilePaths &filePaths,
                 Files &files,
                 double freq[4],
                 double Mij[4][4],
                 double Eij[4][4],
                 std::vector<gsl_rng *> &rngGslvector,
                 std::vector<boost::random::mt19937 *> &rngBoostvector)

{
    int i,j,z;
    std::vector<TreeNode *> nodes;
    //char *newickString2;
    long double totalTreeLength;
    int    HEALTHY_ROOT, TUMOR_ROOT;
    int    numISMdeletions, numISMCNLOH;
    int cumNumMUperTree;
    int     numAltModelSites = 0, numDefaultModelSites = 0, numISMmutations = 0;
    //int altModel = 0;
    
    double cumfreq[4];
    long double cumMij[4][4];
    double cumEij[4][4];
    int  numSNVs = 0;
    long double numMaternalMU = 0;
    long double numPaternalMU = 0;
    
    //double    kappa = 0.0, beta = 0.0;
    long double freqR = 0.0, freqY = 0.0, freqAG = 0.0, freqCT = 0.0;
    //double    Rmat[6], NRmat[12];
    double Cijk[256], Root[4];
    long double *triNucFreq=NULL;
    double   cumNumSNVs = 0;
    double cumNumMU = 0;
    long double cumNumMaternalMU = 0;
    long double cumNumPaternalMU = 0;
    double cumNumDEL = 0;
    double cumNumCNLOH = 0;
    // double cumCountMLgenotypeErrors = 0;
    long double  cumCNLOHbranchLength= 0.0;
    double cumNumMUSq = 0;
    //double cumNumSNVsSq = 0;
    double cumNumDELSq = 0;
    double cumNumCNLOHSq = 0;
    //CellStr     *cell;
    
    
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
    
    
    
    std::vector<SiteStr> allSites(programOptions.numSites);
    
    for (i=0; i< programOptions.numSites; i++)
    {
        allSites[i].alternateAlleles = new int[4];
        if (!allSites[i].alternateAlleles)
        {
            fprintf (stderr, "Could not allocate the allSites[i].alternateAlleles[] structure\n");
            exit (-1);
        }
    }
    
    programOptions.altModelMutationRate = programOptions.mutationRate*programOptions.nonISMRelMutRate;
    
    std::vector<int> SNVsites(programOptions.numSites);
    std::vector<int> SFS(programOptions.numSites);
    std::vector<int> variantSites(programOptions.numSites);
    std::vector<int> DefaultModelSites(programOptions.numSites);
    std::vector<int> AltModelSites(programOptions.numSites);
    
    std::vector<TreeNode *> treeTips;
    std::vector<double> proportionsVector(programOptions.numClones);
    std::vector<int> varEvent(programOptions.numDataSets);
    std::vector<double> varTimeGMRCA(programOptions.numDataSets);
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    TreeNode *init_root = new TreeNode(programOptions.numSites);
    
    double cumNumCA = 0.0, meanNumCA = 0.0, cumNumMIG = 0.0, meanNumMIG = 0.0;
    double numEventsTot = 0.0, countTMRCA = 0.0,TMRCA = 0.0;
    int        numCA = 0,numCNLOH = 0, numMIG = 0;
    int dataSetNum;
    long int seedFirst =  programOptions.seed;
    int  numMU = 0, numDEL = 0, numProposedMU = 0, numFixedMutations = 0, numSNVmaternal = 0;
    programOptions.doUseObservedCellNames=NO;
    //int ***data;
    ValidateParameters(programOptions,CloneNameBegin , CloneSampleSizeBegin, ClonePopSizeBegin);
    
    
    std::vector<MSA*> msaList(programOptions.numDataSets* programOptions.MutationAssignNum);
    std::vector<RootedTree*> treesList(programOptions.numDataSets);
    std::vector<Partition*> partitionList(programOptions.numDataSets* programOptions.MutationAssignNum);
    std::vector<TreeLikelihood*> treeLikList(programOptions.numDataSets* programOptions.MutationAssignNum);
    
    PrepareLikelihoodOutputFile(filePaths, programOptions, files);
    writeHeaderLikelihoodFile(filePaths, programOptions,files, programOptions.numClones );
    
    double totalTimeMRCAPUnits =0.0;
    double totalTimeMRCAModelTime = 0.0;
    double totalRootModelTime = 0.0;
    double totalRootTimePUnits = 0.0;
    boost::random::mt19937 * rngBoost;
    for (dataSetNum = 0; dataSetNum < programOptions.numDataSets; dataSetNum++)// dataSetNum refers to a simulated tree number
    {
        
        rngBoost = rngBoostvector.at(dataSetNum);
        
        if(programOptions.doSimulateFromPriors==YES){
            
            SetPopulationParametersFromPriors( populations, programOptions.numClones,rngGslvector.at(dataSetNum), programOptions);
            
        }
        SetPopulationTimeOriginSTD(populations, programOptions.numClones, rngGslvector.at(dataSetNum), programOptions.doEstimateTimesOriginClones);
        InitNumberNodes( populations, programOptions);
        //        programOptions.meanADOsite = 0.1;
        //        programOptions.varADOsite=0.01;
        //        programOptions.meanADOcell = 0.1;
        //        programOptions.varADOcell=0.01;
        //        programOptions.meanGenotypingError= 0.01;
        //        programOptions.varGenotypingError=0.001;
        
        InitListPossibleMigrations(populations, programOptions.numClones);
        InitPopulationsCoalescentEvents( programOptions.numClones,  populations);
        
        ListClonesAccordingTimeToOrigin(populations, programOptions.numClones);
        
        if (programOptions.doPrintSeparateReplicates == YES)
            PrepareSeparateFiles(0, 1, dataSetNum,  filePaths, programOptions, files, populations);
        /* Reorganize seed per replicate */
        //seed = seedFirst+dataSetNum+10;
        numCA = numMU = numDEL = numProposedMU = TMRCA = numSNVmaternal = 0;
        
        programOptions.seed = seedFirst + dataSetNum * 1000;
        numSNVs=0;
        /* reset variables per replicate */
        if (programOptions.noisy > 0)
        {
            std::cout << "\rReplicate # " << (dataSetNum + 1)<< "/" << programOptions.numDataSets<< std::endl;
            fflush (stdout);
        }
        
        
        varEvent[dataSetNum] = 0;
        varTimeGMRCA[dataSetNum] = 0.0;
        
        numCA = numMIG = 0;
        countTMRCA = 0.0;
        /* coalescent tree */
        
        if (programOptions.noisy > 1)
            std::cout << "\n>> Start coalescent tree .. " << std::endl;
        
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
                                               init_root,
                                               programOptions.K,
                                               rngGslvector.at(dataSetNum),
                                               rngBoost
                                               ) ;
        if (programOptions.noisy > 1)
            std::cout <<  "\n>> Finishing coalescent tree ... DONE"<< std::endl;
        cumNumCA += numCA;
        cumNumMIG += numMIG;
        countTMRCA = root->left->timePUnits;
        
        totalTimeMRCAPUnits +=root->left->timePUnits * programOptions.mutationRate;
        totalTimeMRCAModelTime += root->left->time;
        totalRootModelTime+= root->time;
        totalRootTimePUnits+=root->timePUnits * programOptions.mutationRate;
        
        //fprintf ( stderr, "\n countTMRCA = %lf\n", countTMRCA);
        varTimeGMRCA[dataSetNum] = countTMRCA;
        varEvent[dataSetNum] = numCA + numMIG;
        // counterTime = counterTime + counterTimeInit;
        
        /*************** output files *************/
        //newickString2=NULL;
        // newickString2 = toNewickString2 ( root, programOptions.mutationRate,     programOptions.doUseObservedCellNames);
        // printf("\n newick = %s  \n", newickString2);
        
        if (programOptions.doPrintTrees == YES)
        {
            
            Output::PrintTrees(dataSetNum, root, files.fpTrees->f, programOptions.mutationRate, programOptions.doUseObservedCellNames);
            Output::PrintTrees2(dataSetNum, root, files.fpTrees2->f, programOptions.mutationRate, NULL, NO);
        }
        if (programOptions.doPrintTimes == YES)
        {
            Output::PrintTimes(dataSetNum, files.fpTimes->f, programOptions.mutationRate, nodes, programOptions.thereisOutgroup);
            Output::PrintTimes2(dataSetNum, files.fpTimes2->f, programOptions.mutationRate, nodes, programOptions.thereisOutgroup);
        }
        if (programOptions.doPrintTrees ==YES && programOptions.doPrintSeparateReplicates == YES)
        {
            
            fclose(files.fpTrees->f);
            fclose(files.fpTrees2->f);
        }
        if (programOptions.doPrintTimes ==YES && programOptions.doPrintSeparateReplicates == YES)
        {
            
            fclose(files.fpTimes->f);
            fclose(files.fpTimes2->f);
        }
        
        treesList[dataSetNum]= new RootedTree(files.fpTrees->path, true);
        
        long double logLikCoalTree=0;
        logLikCoalTree+= StructuredCoalescentTree::SumLogDensitiesTimeOriginSTDPopulations(populations, programOptions.numClones, programOptions.K);
        
        logLikCoalTree+= StructuredCoalescentTree::SumLogProbFatherPopulations(populations, programOptions.numClones, programOptions.K);
        
        logLikCoalTree+= StructuredCoalescentTree::LogDensityCoalescentTimesForPopulation(populations, programOptions.numClones, programOptions.K);
        
        if (programOptions.noisy >1)
            std::cout << "\n The coalescent tree likelihood is " << logLikCoalTree<< std::endl;
        
        if (programOptions.noisy > 1)
        {
            fprintf (stderr, "\nData set %d", dataSetNum + 1);
            fprintf (stderr, "\n\tNumber of coalescence events   =   %d", numCA);
            fprintf (stderr, "\n\tNumber of migration events     =   %d", numMIG);
        }
        
        totalTreeLength = SumBranches(root, programOptions.mutationRate, programOptions.healthyTipLabel);
        cumNumMUperTree=0;
        
        //free(newickString2);
        //newickString2=NULL;
        
        for (z = 0; z < programOptions.MutationAssignNum; z++)
        {
            numMU=0;
            numSNVs=0;
            numMaternalMU = 0;
            numPaternalMU = 0;
            
            
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
            
            InitializeGenomes (root, &(programOptions.seed), programOptions.alphabet, programOptions.doUserGenome,programOptions.numSites,  allSites, programOptions.doGeneticSignatures,cumfreq, triNucFreq, NULL,
                               rngGslvector.at(dataSetNum),
                               rngBoost);
            
            HEALTHY_ROOT=root->label;
            TUMOR_ROOT=root->label;
            
            //if (SNPrate > 0)
            //           AddGermlineVariation (treeRootInit[0], &(programOptions.seed),  programOptions.numSites, SNPrate, allSites, programOptions.alphabet,  data,   HEALTHY_ROOT, cumMij );
            
            EvolveSitesOnTree (root, MATERNAL, &(programOptions.seed), programOptions.rateVarAmongSites,  programOptions.numSites,  allSites, programOptions.doGeneticSignatures, programOptions.alphaSites, programOptions.propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , numISMmutations, programOptions.numFixedMutations, numSNVmaternal,  programOptions.doSimulateFixedNumMutations,  programOptions.alphabet,  numMU, cumMij,  programOptions.altModel, programOptions.altModelMutationRate, programOptions.doUserTree,  programOptions.doJC,  programOptions.doHKY,  programOptions.doGTR,
                               programOptions.doGTnR,  freqR,  freqY,
                               freqAG, freqCT, programOptions.titv, freq, Mij ,   Root,  Cijk, rngGslvector.at(dataSetNum),   rngBoost);
            EvolveSitesOnTree (root, PATERNAL, &(programOptions.seed), programOptions.rateVarAmongSites,  programOptions.numSites,  allSites, programOptions.doGeneticSignatures, programOptions.alphaSites, programOptions.propAltModelSites ,  numDefaultModelSites, numAltModelSites, DefaultModelSites, AltModelSites,  totalTreeLength , numISMmutations, numFixedMutations, numSNVmaternal,  programOptions.doSimulateFixedNumMutations,  programOptions.alphabet,  numMU, cumMij,  programOptions.altModel, programOptions.altModelMutationRate, programOptions.doUserTree,  programOptions.doJC,  programOptions.doHKY,  programOptions.doGTR,
                               programOptions.doGTnR,  freqR,  freqY,
                               freqAG, freqCT, programOptions.titv, freq, Mij,   Root,  Cijk, rngGslvector.at(dataSetNum),   rngBoost );
            cumNumMU += numMU;
            cumNumMUSq += pow(numMU,2);
            
            computeStatisticsNumberMutations(allSites, numMaternalMU, numPaternalMU);
            
            cumNumMaternalMU+=numMaternalMU;
            cumNumPaternalMU+=numPaternalMU;
            
            
            
            if (programOptions.doPrintSeparateReplicates == YES)
                PrepareSeparateFilesGenotypes(1, dataSetNum, z,
                                              filePaths, programOptions,files, populations, 1.0*numMU/programOptions.numSites );
            
            if (programOptions.CNLOHrate > 0)
            {
                /* evolve maternal CN_LOH */
                if (programOptions.noisy > 2){
                    std::cout << "\n>> Evolving maternal CN_LOH ... "<<std::endl;
                    EvolveCNLOHonTree (root, MATERNAL, numISMCNLOH, allSites, programOptions.numSites, &(programOptions.seed), programOptions.CNLOHrate,  programOptions.mutationRate,  totalTreeLength, cumCNLOHbranchLength,  rngGslvector.at(dataSetNum),  rngBoost );
                }
                
                /* evolve paternal CN_LOH  */
                if (programOptions.noisy > 2){
                    std::cout << "\n>> Evolving paternal CN_LOH ... " <<std::endl;
                    EvolveCNLOHonTree (root, PATERNAL, numISMCNLOH, allSites, programOptions.numSites, &(programOptions.seed),  programOptions.CNLOHrate,  programOptions.mutationRate,  totalTreeLength, cumCNLOHbranchLength,  rngGslvector.at(dataSetNum),  rngBoost );
                }
                
                cumNumCNLOH += numCNLOH;
                cumNumCNLOHSq += pow(numCNLOH,2);
            }
            
            if (programOptions.deletionRate > 0)
            {
                /* evolve maternal deletions */
                if (programOptions.noisy > 2)
                    std::cout <<"\n>> Evolving maternal deletions ... "<< std::endl;
                EvolveDeletionsOnTree (root, MATERNAL, allSites, numISMdeletions, programOptions.numSites,   totalTreeLength, programOptions.mutationRate, programOptions.deletionRate, numDEL,  &(programOptions.seed), rngGslvector.at(dataSetNum),  rngBoost);
                
                /* evolve paternal deletions  */
                if (programOptions.noisy > 2)
                    std::cout << "\n>> Evolving paternal deletions ... "<< std::endl;
                
                EvolveDeletionsOnTree (root, PATERNAL, allSites, numISMdeletions, programOptions.numSites,   totalTreeLength, programOptions.mutationRate, programOptions.deletionRate, numDEL,  &(programOptions.seed), rngGslvector.at(dataSetNum),  rngBoost);
                
                
                cumNumDEL += numDEL;
                cumNumDELSq += pow(numDEL,2);
            }
            
            
            numSNVs = countTrueVariants (nodes, programOptions.numSites, programOptions.numCells, root, allSites, variantSites, SNVsites );
            
            cumNumSNVs+=numSNVs;
            
            if (programOptions.doPrintTrueHaplotypes == YES)
            {
                if (programOptions.doPrintSeparateReplicates == NO)
                    fprintf (files.fpTrueHaplotypes->f, "[#%d]\n", z+1);
                //
                Output::PrintTrueFullHaplotypes (files.fpTrueHaplotypes->f,  nodes, root , programOptions.numNodes, programOptions.doPrintIUPAChaplotypes, programOptions.doPrintAncestors, programOptions.numSites,  programOptions.numCells, programOptions.alphabet, programOptions.doUserTree,    programOptions.doNGS,   NULL, NULL, HEALTHY_ROOT, TUMOR_ROOT, NULL, NO);
            }
            
            
            if (programOptions.doPrintTrueHaplotypes ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpTrueHaplotypes->f);
            }
            
            pll_phylip_t * phylip_file_true= pll_phylip_open(files.fpTrueHaplotypes->path,
                                                             pll_map_phylip);
            
            if (!phylip_file_true)
                fprintf (stderr, "\n ERROR: phylip file cannot be opened  to compute  Felsenstein likelihood!");
            int pos = (dataSetNum)*programOptions.MutationAssignNum+z;
            assert(pos>=0);
            msaList[pos] = new MSA( pll_phylip_parse_interleaved(phylip_file_true));
            
            
            /* compute true likelihoods: coalescent tree and sequence likelihoods*/
            
            double logLikGenotypes = 0.0;
            GenotypeErrorModel gtNoError("GT20", 0.0, 0.0, 16);
            
            if  (programOptions.computeLikelihoods){
            
            
            partitionList[pos] = new Partition(treesList[dataSetNum]->numTips(),// numberTips
                                               treesList[dataSetNum]->numInner(),//unsigned  int  clvBuffers
                                               16,// model->states,//numberStates
                                               (unsigned int)(msaList[pos]->getLength()),//unsigned  int  sites
                                               1,//unsigned  int numberRateMatrices
                                               treesList[dataSetNum]->numBranches(), // unsigned int probMatrices
                                               RATE_CATS,//RATE_CATS, // unsigned  int  numberRateCats
                                               treesList[dataSetNum]->numInner(), // unsigned  int numberScaleBuffers
                                               0, //int statesPadded
                                               false, false, false, false, false, false);
            
//            partitionList.emplace_back(
//                                       treesList[dataSetNum]->numTips(),// numberTips
//            treesList[dataSetNum]->numInner(),//unsigned  int  clvBuffers
//            16,// model->states,//numberStates
//            (unsigned int)(msaList[pos]->getLength()),//unsigned  int  sites
//            1,//unsigned  int numberRateMatrices
//            treesList[dataSetNum]->numBranches(), // unsigned int probMatrices
//            RATE_CATS,//RATE_CATS, // unsigned  int  numberRateCats
//            treesList[dataSetNum]->numInner(), // unsigned  int numberScaleBuffers
//            0, //int statesPadded
//            false, false, false, false, false, false);
//
           treeLikList[dataSetNum] = new TreeLikelihood(*(partitionList[pos]), *(treesList[dataSetNum]),  *(msaList[pos]), gtNoError)
            ;

            
             logLikGenotypes = treeLikList[pos]->computeRootLogLikelihood();
                
            // assert(!isnan(logLikGenotypes) && !isinf(logLikGenotypes));
                      
            if (programOptions.noisy >1)
                          std::cout << "\n The sequences likelihood without errors is " << logLikGenotypes<< std::endl;
            
            }
         
            
            //            if (msaList[pos]==NULL )
            //                fprintf (stderr, "\n ERROR: phylip file cannot be opened  to compute  Felsenstein likelihood!");
            
            
            if (programOptions.fixedADOrate > 0 || programOptions.doADOcell == YES || programOptions.doADOsite == YES)
                AllelicDropout (programOptions.TotalNumSequences, allSites,  programOptions.doADOcell,  programOptions.doADOsite,
                                programOptions.numSites, programOptions.fixedADOrate,  programOptions.meanADOcell, programOptions.varADOcell, programOptions.meanADOsite, programOptions.varADOsite, treeTips,   &(programOptions.seed), rngGslvector.at(dataSetNum),   rngBoost);
            
            
            /* introduce errors directly in the genotypes */
            if ( programOptions.meanGenotypingError > 0)
                GenotypeError (treeTips, allSites,  programOptions.alphabet,  programOptions.numSites,  programOptions.numCells, programOptions.meanGenotypingError , programOptions.varGenotypingError, Eij,  &(programOptions.seed),  rngGslvector.at(dataSetNum),  rngBoost);
            
            
            
            if (programOptions.doPrintFullGenotypes == YES)
            {
                if (programOptions.doPrintSeparateReplicates == NO)
                    fprintf (files.fpFullGenotypes->f, "[#%d]\n", z+1);
                Output::PrintFullGenotypes(files.fpFullGenotypes->f, nodes, root , programOptions.numNodes, programOptions.doPrintIUPAChaplotypes, programOptions.doPrintAncestors, programOptions.numSites,  programOptions.numCells, programOptions.alphabet, programOptions.doUserTree,    programOptions.doNGS,   NULL, NULL, HEALTHY_ROOT, TUMOR_ROOT, NULL, NO,  numSNVs, SNVsites);
            }
            if (programOptions.doPrintSNVgenotypes == YES && numSNVs > 0) /* we only print replicates with variation */
            {
                if (programOptions.doPrintSeparateReplicates == NO)
                    fprintf (files.fpSNVgenotypes->f, "[#%d]\n", z+1);
                Output::PrintSNVGenotypes(files.fpSNVgenotypes->f, nodes, root , programOptions.numNodes, programOptions.doPrintIUPAChaplotypes, programOptions.doPrintAncestors, programOptions.numSites,  programOptions.numCells, programOptions.alphabet, programOptions.doUserTree,    programOptions.doNGS,   NULL, NULL, HEALTHY_ROOT, TUMOR_ROOT, NULL, NO,  numSNVs, SNVsites);
            }
            
            if (programOptions.doPrintMLhaplotypes ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpMLhaplotypes->f);
                
            }
            if (programOptions.doPrintSNVgenotypes == YES && numSNVs > 0){
                
                fclose(files.fpSNVgenotypes->f);
                
            }
            if (programOptions.doPrintSNVhaplotypes ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpSNVhaplotypes->f);
                
            }
            if (programOptions.doPrintFullGenotypes ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpFullGenotypes->f);
                
            }
            
            if (programOptions.doPrintCATG ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpCATG->f);
            }
            if (programOptions.doSimulateReadCounts ==YES && programOptions.doPrintSeparateReplicates == YES)
            {
                
                fclose(files.fpVCF->f);
                
            }
            
            
            if  (programOptions.computeLikelihoods){
            pll_phylip_t * phylip_file_errors= pll_phylip_open(files.fpFullGenotypes->path,
                                                               pll_map_phylip);
            
            if (!phylip_file_errors)
                fprintf (stderr, "\n ERROR: phylip file cannot be opened  to compute  Felsenstein likelihood!");
            
            //            pll_msa_t *msa_errors =pll_phylip_parse_interleaved(phylip_file_errors);
            //
            //            if (!msa_errors)
            //              fprintf (stderr, "\n ERROR: phylip file cannot be opened  to compute  Felsenstein likelihood!");
            //
            
            
            GenotypeErrorModel gtErrorModel("GT20", programOptions.meanGenotypingError,  1.0 - sqrt (1.0 - programOptions.fixedADOrate), 16);
            
            
            
            treeLikList[pos]->changeGenotypeErrorModel(&gtErrorModel);
            
            double logLikGenotypeErrors = treeLikList[pos]->computeRootLogLikelihood();
            
            if (programOptions.noisy >1)
                std::cout << "\n The sequences likelihood with errors is " << logLikGenotypeErrors<< std::endl;
            
            writeLineLikelihoodFile( dataSetNum, filePaths, programOptions,files , populations,
                                    logLikCoalTree,  logLikGenotypes,  logLikGenotypeErrors, totalTreeLength);
            
 
            pll_phylip_close(phylip_file_true);
            pll_phylip_close(phylip_file_errors);
            
            }
          
            
        }/* end of mutation simulation process */
        
    }
    
    std::cout << "\nSimulation  statistics: \n "<< std::endl;
    std::cout << "\n The average time of MRCA in  time scaled by mutation rate is " << totalTimeMRCAPUnits/ programOptions.numDataSets << std::endl;
    
    std::cout << "\n The average time of MRCA in model time is "<< totalTimeMRCAModelTime/ programOptions.numDataSets<< std::endl;
    
    std::cout << "\n The average  time of tree root in time scaled by mutation rate is " << totalRootTimePUnits/ programOptions.numDataSets << std::endl;
    
    std::cout << "\n The average time of tree root in model time is "<< totalRootModelTime/ programOptions.numDataSets<< std::endl;
    
    std::cout << "\n The average number of mutations per dataset is  "<< cumNumMU/ programOptions.numDataSets<< std::endl;
    
    std::cout << "\n The average number of mutations per dataset per site is  "<< cumNumMU/ (programOptions.numDataSets *programOptions.numSites)<< std::endl;
    
    std::cout << "\n The average number of variable sites  is  "<< cumNumSNVs/ (programOptions.numDataSets)<< std::endl;
    
    std::cout << "\n The average number of maternal mutations per  branch per site is  "<< cumNumMaternalMU/ (programOptions.numDataSets*programOptions.numSites*(2*(programOptions.TotalNumSequences+1)-2))<< std::endl;
    
    std::cout << "\n The average number of paternal mutations per  branch per site is  "<< cumNumPaternalMU/ (programOptions.numDataSets*programOptions.numSites*(2*(programOptions.TotalNumSequences+1)-2))<< std::endl;
    
    fclose(files.fpLikelihood->f);
    
    
    /*release memory*/
    for (i=0; i< programOptions.numSites; i++)
    {
        free( allSites[i].alternateAlleles);
    }
    
    
    
    for (auto ptr : partitionList)
    {
        delete ptr;
    }
    partitionList.clear();
    
    for (auto ptr : treesList)
    {
        delete ptr;
    }
    treesList.clear();
    

    treeLikList.clear();
    
  
    
    return 0;
}
/***************************** ValidateParameters*******************************/
/* Validate parameters*/

void ValidateParameters(ProgramOptions &programOptions,
                        std::vector<int> CloneNameBegin , std::vector<int> CloneSampleSizeBegin, std::vector<int> ClonePopSizeBegin)
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
            Output::PrintUsage();
        }
        if ( CloneSampleSizeBegin[j] > ClonePopSizeBegin[j] )
        {
            fprintf (stderr, "PARAMETER ERROR: Clone (%d) cannot have sample size (%d) higher than population size (%d). (d)\n\n", j, CloneSampleSizeBegin[j] , ClonePopSizeBegin[j]);
            Output::PrintUsage();
        }
    }
}
/***************** bbinClones *****************/
/* binary search in the probabilities with clones */
int bbinClones (long double dat, long double *v, int n)
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
void  InitListPossibleMigrations(std::vector<Population *> &populations, int numClones)
{
    sort(populations.begin(), populations.end(), comparePopulationsByTimeOrigin);
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        p->InitListPossibleMigrations(i);
    }
}
void SetPopulationTimeOriginSTD(std::vector<Population *> &populations, int numClones,const gsl_rng* rngGsl, bool doEstimateTorigins){
    
    if (doEstimateTorigins){
        Population *p;
        for (size_t i = 0; i < numClones; ++i) {
            p = populations[i];
            p->setPopulationToriginConditionalDelta(rngGsl);
        }
    }
}
void InitNumberNodes( std::vector<Population *> &populations, ProgramOptions &programOptions) {
    programOptions.TotalNumSequences = 0;
    
    Population *popI;
    int j;
    for (j = 0; j < programOptions.numClones; j++)
    {
        popI = populations[j];
        popI->FatherPop = NULL;
        programOptions.TotalNumSequences = programOptions.TotalNumSequences + popI->sampleSize;
        //TotalN = TotalN + popI->popSize;
        
    }
    programOptions.numNodes = 2 * programOptions.TotalNumSequences + programOptions.numClones+ 10;
    
    programOptions.numCells = programOptions.TotalNumSequences;
}
void SetPopulationParametersFromPriors(std::vector<Population *> &populations, int numClones,const gsl_rng* rngGsl, ProgramOptions &programOptions){
    
    Population *p;
    programOptions.mutationRate=Random::RandomExponential(1, NULL, true, rngGsl, NULL);
    unsigned int  total_sample = 0;
    for (size_t i = 0; i < numClones; ++i) {
        p = populations[i];
        int maxSampleSize = ((programOptions.maxSampleSize /numClones) < programOptions.minSampleSize)? 2* programOptions.minSampleSize: (programOptions.maxSampleSize /numClones) ;
        p->sampleSize = Random::randomUniformIntegerInterval(rngGsl, programOptions.minSampleSize, maxSampleSize );
        p->numGametes = p->sampleSize;
        p->numActiveGametes = p->sampleSize;
        p->delta = Random::RandomExponential(0.01, NULL, true, rngGsl, NULL);
        
        total_sample+= p->sampleSize;
    }
    for (size_t i = 0; i < numClones; ++i) {
        p = populations[i];
        p->x  =  p->sampleSize / (1.0*total_sample);
        p->theta =  programOptions.mutationRate;
    }
}
/********************** resetMigrationsList ************************/
/* reset list of migrations*/
void resetMigrationsList(std::vector<Population*> &populations, int numClones){
    Population *p;
    int i;
    for (i = 0; i < numClones; i++) {
        p = populations[i];
        p->resetMigrationsList();
    }
}
/********************** ChooseFatherPopulation************************/
/* choose probabilistically  the father population of a  population  */
Population* ChooseFatherPopulation(std::vector<Population*> &populations, int numClones, Population  *PopChild,  long int *seed, int noisy, long double K, const gsl_rng* rngGsl) {
    
    Population  *p;
    double pij, ran;
    int  j, k;
    double cumProb[numClones  - (int)(PopChild->order)];
    cumProb[0] = 0.0;
    for (j = PopChild->order + 1; j < numClones; j++)
    {
        cumProb[j - PopChild->order] = 0.0;
        p = populations[j];
        pij = ProbabilityCloneiFromClonej2(PopChild, p, populations, numClones, K);
        
        if (noisy>1){
        std::cout<<  "\n the probability that  population with order  "<<p->order << " is the source of population with order "<< PopChild->order << " is " << pij<< std::endl;
        }
        cumProb[j - PopChild->order] = cumProb[j - 1 - PopChild->order] + pij;
        
        
    }
    // now selecting the ancestral clone
    ran = Random::randomUniformFromGsl2(rngGsl);
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
    if (noisy > 1)
        std::cout<<  "\nClone with order "<< PopChild->order <<" derived from clone with order "<< result->order << std::endl; // clone ThisCloneNumber (i) is originated from clone ThisOriginCloneNumber (j)
    /* Update list of migation times considering that clone i comes from clone j */
    if (noisy > 1)
        std::cout<< "\n*** Updating list of migration times (considering that clone "<< PopChild->index <<  " comes from clone " << result->index <<") .." << std::endl;
    return (result); //the father population has  order  (PopChild->order) + w
}
/********************** AssignCurrentSequencesToPopulation************************/
/* assign current sequences to  population  */
void AssignCurrentSequencesToPopulation(std::vector<Population *> &populations, std::vector<TreeNode*> &nodes,
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
    
    std::vector<int> CumSumNodes(numClones+1);
    
    
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
void ChooseRandomIndividual(int *firstInd,   int numClones, Population *popI,  int *secondInd, long *seed, int choosePairIndividuals, const gsl_rng* rngGsl)
{
    long double random;
    int k;
    long double *cumPopulPart = (long double *) malloc((popI->numActiveGametes + 1)* (long) sizeof(long double));
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
    random = Random::randomUniformFromGsl2(rngGsl);
    //fprintf (stderr, "\nran = %lf ", ran);
    *firstInd = bbinClones(random, cumPopulPart, popI->numActiveGametes)-1;
    
    
    if (*firstInd >= popI->numActiveGametes || *firstInd < 0 ) /* checking */
    {
        fprintf (stderr, "\n\nERROR: firstInd out of range!\n");
        exit (-1);
    }
    
    if (choosePairIndividuals== YES && popI->numActiveGametes > 1) {
        
        do//choose randomly another individual to coalesce
        {
            random = Random::randomUniformFromGsl2(rngGsl);
            *secondInd = bbinClones(random, cumPopulPart, popI->numActiveGametes)-1;
            
        } while (*firstInd == *secondInd  );
    }
    free (cumPopulPart);
    cumPopulPart=NULL;
}

/********************** MakeCoalescenceEvent************************/
/*  choose 2  active individuals  to make coalescent  */
void MakeCoalescenceEvent(std::vector<Population*> &populations, Population *popI, std::vector<TreeNode *> &nodes, int numClones, long int* seed, int noisy,   int &numActiveGametes, int &nextAvailable,
                          int &labelNodes, double &currentTime, int &numNodes, const gsl_rng * rngGsl)
{

    TreeNode  *p, *q, *r;
    int firstInd, secondInd=0, newInd=0;
    int choosePairIndividuals = YES;
    
    ChooseRandomIndividual(&firstInd, numClones, popI,  &secondInd, seed, choosePairIndividuals, rngGsl);
    
    newInd = nextAvailable;
    if (noisy > 1)
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
    r->nodeClass = 4;
    // link the nodes
    r->left = p;
    r->right = q;
    p->anc1 = r;
    q->anc1 = r;
    r->time = currentTime;
    
    assert(popI->x > 0);
    //if (popI->order == numClones-1)//oldest population or only one population and order=0
     //   r->timePUnits = currentTime;
    //else
        r->timePUnits = currentTime * (popI->x);

    if (noisy > 1)
        fprintf (stderr, "\t|\tCurrentTime (input units) = %Lf", r->timePUnits);
    /* readjust active nodes */
    
    popI->idsActiveGametes[firstInd] = newInd;
    popI->idsActiveGametes[secondInd] = popI->idsActiveGametes[popI->numActiveGametes - 1];;
    numActiveGametes = numActiveGametes - 1; /* less 1 active node */
    

    //update list ids nodes
    // popI->idsGametes[popI->numGametes] = newInd;
    popI->numGametes = popI->numGametes +1;
    
    nextAvailable=nextAvailable+1; /* 1 node more is available */
    
    popI->CoalescentEventTimes[ popI->numCompletedCoalescences]= r->timePUnits;//  r->time;
    popI->numActiveGametes = popI->numActiveGametes - 1; /* now this clone
                                                          has 1 less node */
    
    popI->numCompletedCoalescences= popI->numCompletedCoalescences+1;
    
    //fprintf (stderr, "\n pop of order  %d, number of Active gametes %d ", popI->order, popI->numActiveGametes);
    
    // for(int i=0; i < popI->idsActiveGametes.size();i++)
    //     fprintf (stderr, "\n pop of order  %d Active gamete id %d", popI->order, popI->idsActiveGametes[i]);
    /* memory for number of nodes */
    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
    {
        /* ReallocNodes(&numNodes, activeGametes); */
        if (noisy == 4)
            fprintf (stderr, "\n\n...Doing reallocation of nodes (Coalescence)\n");
        numNodes += INCREMENT_NODES;
        
    }
}
/********************** BuildTree************************/
/*  build tree */
TreeNode *BuildTree(std::vector<Population* > &populations,
                    Population *CurrentPop,
                    long int *seed,
                    ProgramOptions &programOptions,
                    std::vector<TreeNode *> &nodes,
                    std::vector<TreeNode *> &treeTips,
                    TreeNode *tumour_mrca,
                    int &nextAvailable,
                    int &newInd,
                    double &currentTime,// model time oldest population
                    int &labelNodes,
                    const gsl_rng *rngGsl)
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
                fprintf (stderr, "%d %d-- %d %d %d\n", Output::Index(p), j, Output::Index(p->left), Output::Index(p->right), Output::Index(p->anc1));
            }
            if (p->anc1 != NULL)
            {//update length field
                
                //   p->length = p->anc1->time- p->time;
                p->length = (p->anc1->timePUnits- p->timePUnits);// * programOptions.mutationRate;
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
            currentTime = (CurrentPop->timeOriginSTD -programOptions.outgroupBranchLength_Root1Root2 * tumour_mrca->time) / (1.0- programOptions.outgroupBranchLength_Root1Root2 )  ; // origin of the clone + time given by the user
            
        }
        else
        {
            fprintf (stderr, "\n\nError simulationg the outgroup. Check input settings\n");
            Output::PrintUsage();
        }
        //TreeNode* healthyRoot = *nodes + *nextAvailable;
        TreeNode* healthyRoot = nodes[nextAvailable];
        healthyRoot->index = nextAvailable;
        healthyRoot->label = labelNodes;
        labelNodes=labelNodes+1;
        healthyRoot->left = p;//coalTreeMRCA;
        //        coalTreeMRCA->anc = healthyRoot;
        p->anc1 = healthyRoot;
        
        // healthyRoot->timePUnits = p->timePUnits * healthyRoot->effectPopSize;
        healthyRoot->nodeClass = 5;
        //        coalTreeMRCA->length = transformingBranchLength/mutationRate;
        p->length = 0;
        
        //        coalTreeMRCA->branchLength = transformingBranchLength;
        p->lengthModelUnits = 0;
        
        //        healthyRoot->time = currentTime +  transformingBranchLength/mutationRate;
        healthyRoot->time = currentTime;
        
        //int transformingBranchLength=1.001;
        // healthyRoot->time = p->time * transformingBranchLength ;
        
        //healthyRoot->timePUnits = currentTime * healthyRoot->effectPopSize;
        
        //healthyRoot->timePUnits = currentTime;
        healthyRoot->timePUnits = currentTime * CurrentPop->x;
        
        //fprintf (stderr, "\n Time of the healthy root %Lf\n",  healthyRoot->timePUnits);
        
        p->length = (p->anc1->timePUnits- p->timePUnits) ;//* programOptions.mutationRate;
        //*mutationRate;
        p->lengthModelUnits = (p->anc1->time- p->time);
        //*mutationRate;
        
        healthyRoot->length = 0;
        //        healthyRoot->length = 0;
        
        
        nextAvailable++;
        //
        //        /* connect the healthy ancestral cell with the tip healthy cell*/
        
        TreeNode* healthyTip = nodes[nextAvailable];
        healthyTip->left = NULL;
        healthyTip->right = NULL;
        
        //connectNodes(NULL, NULL, healthyTip);
        
        healthyTip->anc1 = healthyRoot;
        healthyRoot->right = healthyTip;
        
        double  healthyTipBranchLengthRatio = Random::randomUniformFromGsl2(rngGsl);
        
        //we put the  time of healthy tip inside the tree root and 0
        healthyTip->time = healthyTipBranchLengthRatio * healthyRoot->time;
        healthyTip->timePUnits = healthyTipBranchLengthRatio * healthyRoot->timePUnits;
        healthyTip->length = 0;
        healthyTip->lengthModelUnits = 0;
        
        healthyTip->isOutgroup= YES;
        
        treeRootInit = healthyRoot;
        
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

void PrepareSeparateFiles(int ChainNumber, int paramSetNumber, int replicate,const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, std::vector<Population*> &populations)
{
    char File[1000];
    char dir[1000];
    /* contains the simulated tree in Newick format
     */
    //    if (boost::filesystem::exists( "Results"))
    //           boost::filesystem::remove_all("Results");
    //       boost::filesystem::create_directory("Results");
    //mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
    //mkdir("Results",0);
#ifdef MAC
    strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
    strcpy (dir, "Results/");
#endif
    //strcpy (resultsDir, dir);
   std::stringstream ss;
    Population *popI;
      for( size_t i = 0 ; i < populations.size(); i++)
      {
          popI = populations[i];
          ss << std::fixed << std::setprecision(2) <<"G" << i << "=" << popI->delta << "_";
          ss << "T" << i << "=" << popI->timeOriginSTD << "_";
          ss << "n" << i << "=" << popI->sampleSize << "_";
          ss << "m" << i << "=" << popI->theta << "_";
          ss << "x" << i << "=" << popI->x;
          ss << "_";
          
      }
    ss << "m=" << populations[0]->theta;
   
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
            sprintf(File,"%s/%s/%s_model_time_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile, paramSetNumber+1, replicate+1);
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (populations.size()==1){
            sprintf(File,"%s/%s/%s_model_time_Gamma=%.3Lf_T=%.3Lf_theta=%.3Lf_epsilon=%.3f_n=%d_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,programOptions.seqErrorRate,populations[0]->sampleSize, paramSetNumber+1, replicate+1);
                }
                else{
                    
                    std::string s = ss.str();
            //sprintf(File,"%s/%s/%s_%s_e=%.3f_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile,"trees",programOptions.seqErrorRate, paramSetNumber+1, replicate+1);
                    
                    sprintf(File,"%s/%s/%s_model_time_%s_e=%.3f_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile,s.c_str(),programOptions.seqErrorRate, paramSetNumber+1, replicate+1);
                    
                }
            }
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d/%s_model_time_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir,ChainNumber , filePaths.treeFile, paramSetNumber+1, replicate+1);
        }
        //if ((*fpTrees = fopen(File, "w")) == NULL)
        strcpy (files.fpTrees->path,File);
        if (openFile(&files.fpTrees->f, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        
        sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.treeDir);
        
        if (programOptions.doSimulateData ==YES)
        {
            sprintf(File,"%s/%s/%s_physical_time_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile, paramSetNumber+1, replicate+1);
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (populations.size()==1){ sprintf(File,"%s/%s/%s_physical_time_G=%.3Lf_T=%.3Lf_theta=%.3Lf_e=%.3f_n=%d_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,programOptions.seqErrorRate,populations[0]->sampleSize, paramSetNumber+1, replicate+1);
                }
                else{
                                   
                         std::string s = ss.str();
                           sprintf(File,"%s/%s/%s_physical_time_%s_e=%.3f_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, filePaths.treeFile,s.c_str(), programOptions.seqErrorRate, paramSetNumber+1, replicate+1);
                                   
                        }
               
            }
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d/%s_physical_time_%04d_%04d.tre", filePaths.resultsDir, filePaths.treeDir, ChainNumber, filePaths.treeFile, paramSetNumber+1, replicate+1);
            
        }
        //sprintf(File,"%s/%s/%s_2_%04d.tre", resultsDir, treeDir, treeFile, replicate+1);
        //if ((*fpTrees2 = fopen(File, "w")) == NULL)
        strcpy (files.fpTrees2->path,File);
        if (openFile(&files.fpTrees2->f, File) == -1)
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
            sprintf(File,"%s/%s/%s_model_time_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile, paramSetNumber+1, replicate+1);
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (populations.size()==1){ sprintf(File,"%s/%s/%s_model_time_Gamma=%.3Lf_T=%.3Lf_theta=%.3Lf_epsilon=%.3f_n=%d_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,programOptions.seqErrorRate,populations[0]->sampleSize, paramSetNumber+1, replicate+1);
                }
               else{
                                               
                    std::string s = ss.str();
                    sprintf(File,"%s/%s/%s_model_time_%s_e=%.3f_%04d_%04d.tre", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile,s.c_str(), programOptions.seqErrorRate, paramSetNumber+1, replicate+1);
                                               
                        }
            }
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d/%s_model_time_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, ChainNumber, filePaths.timesFile, paramSetNumber+1, replicate+1);
            
        }
        //if ((*fpTimes = fopen(File, "w")) == NULL)
        strcpy (files.fpTimes->path,File);
        if (openFile(&files.fpTimes->f, File) == -1)
        {
            fprintf (stderr, "Can't open \"%s\"\n", File);
            exit(-1);
        }
        sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.timesDir);
        if (programOptions.doSimulateData ==YES)
        {
            sprintf(File,"%s/%s/%s_physical_time_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile, paramSetNumber+1, replicate+1);
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (populations.size()==1)
                { sprintf(File,"%s/%s/%s_physical_time_Gamma=%.3Lf_T=%.3Lf_theta=%.3Lf_epsilon=%.3f_n=%d_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,programOptions.seqErrorRate,populations[0]->sampleSize, paramSetNumber+1, replicate+1);
                }
                else{
                                                            
                    std::string s = ss.str();
                    sprintf(File,"%s/%s/%s_physical_time_%s_e=%.3f_%04d_%04d.tre", filePaths.resultsDir, filePaths.timesDir, filePaths.timesFile,s.c_str(), programOptions.seqErrorRate, paramSetNumber+1, replicate+1);
                                                            
                    }
            }
        }
        else
        {
            sprintf(File,"%s/%s_Chain%d/%s_physical_time_%04d_%04d.txt", filePaths.resultsDir, filePaths.timesDir, ChainNumber,filePaths.timesFile, paramSetNumber+1, replicate+1);
            
        }
        // sprintf(File,"%s/%s/%s_2_%04d.txt", resultsDir, timesDir, timesFile, replicate+1);
        
        //if ((*fpTimes2 = fopen(File, "w")) == NULL)
        strcpy (files.fpTimes2->path,File);
        if (openFile(&files.fpTimes2->f, File) == -1)
        {
            fprintf(stderr, "Can't open %s.\n", File);
            exit(-1);
        }
    }
    
    
}
/********************* PrepareLikelihoodOutputFile **********************/
/* Open file for writing likelihood results */
void PrepareLikelihoodOutputFile(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files){
    
    char File[MAX_NAME];
    char dir[MAX_NAME];
    /* contains the simulated tree in Newick format
     */
    if (boost::filesystem::exists( "Results"))
        boost::filesystem::remove_all("Results");
    boost::filesystem::create_directory("Results");
    // mkdir("Results", S_IRWXU);
    
    /* Create "Results" folder (with type S_IRWXU (read, write and
     execute)) */
    //mkdir("Results",0);
#ifdef MAC
    strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
    strcpy (dir, "Results/");
#endif
    //strcpy (resultsDir, dir);
    
    
    sprintf(File,"%s/%s.txt", filePaths.resultsDir, filePaths.likelihoodOuput);
    
    strcpy (files.fpLikelihood->path,File);
    
    
    if (openFile(&files.fpLikelihood->f, File) == -1)
    {
        fprintf (stderr, "Can't open \"%s\"\n", File);
        exit(-1);
    }
    
}
/********************* PrepareTempFileInputStan **********************/
/* Open file for writing temp file for stan input  */
void PrepareTempFileInputStan(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int iter){
    
    char File[MAX_NAME];
    char dir[MAX_NAME];
    /* contains the simulated tree in Newick format
     */
    if (boost::filesystem::exists( "Results"))
        boost::filesystem::remove_all("Results");
    boost::filesystem::create_directory("Results");
    // mkdir("Results", S_IRWXU);
    
    /* Create "Results" folder (with type S_IRWXU (read, write and
     execute)) */
    //mkdir("Results",0);
#ifdef MAC
    strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
    strcpy (dir, "Results/");
#endif
    //strcpy (resultsDir, dir);
    
    
    sprintf(File,"%s/%s_%d.txt", filePaths.resultsDir, filePaths.tempInputStan, iter);
    
    strcpy (files.fpStanDump->path,File);
    
    
    if (openFile(&files.fpStanDump->f, File) == -1)
    {
        fprintf (stderr, "Can't open \"%s\"\n", File);
        exit(-1);
    }
    
}
void writeHeaderLikelihoodFile(const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, int numClones )
{
    std::string paramName;
    std::time_t timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    
    fprintf (files.fpLikelihood->f, "#Tumor_coal\t");
    fprintf (files.fpLikelihood->f, "#Generated %s\t\n", std::ctime(&timeNow));
    
    fprintf (files.fpLikelihood->f, "sim\t");
    fprintf (files.fpLikelihood->f, "nPop\t");
    fprintf (files.fpLikelihood->f, "log_lik_tree\t");
    fprintf (files.fpLikelihood->f, "log_lik_true_seq\t");
    fprintf (files.fpLikelihood->f, "log_lik_error_seq\t");
    fprintf (files.fpLikelihood->f, "fixed_ADO_rate\t");
    fprintf (files.fpLikelihood->f, "mean_cell_ADO_rate\t");
    fprintf (files.fpLikelihood->f, "var_cell_ADO_rate\t");
    fprintf (files.fpLikelihood->f, "mean_gen_error_rate\t");
    fprintf (files.fpLikelihood->f, "var_gen_error_rate\t");
    fprintf (files.fpLikelihood->f, "tot_tree_length\t");
    
    paramName =  "Theta";
    fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
    
    for (size_t  i = 0; i < numClones; ++i)
    {
        paramName =  "Delta_pop_" + std::to_string(i+1);
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        
        paramName =  "prop_effect_pop_size_pop_" + std::to_string(i+1) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        paramName =  "DeltaT_pop_" + std::to_string(i+1) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        paramName =  "Theta_pop_" + std::to_string(i+1) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        //paramName =  "mut_rate(pop" + std::to_string(i) + ")";
        //fprintf (files.fplog->f, "%s\t", paramName.c_str());
        paramName =  "sample_size_pop_" + std::to_string(i+1) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        paramName =  "T_input_pop_" + std::to_string(i) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        paramName =  "T_std_pop_" + std::to_string(i+1);
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        paramName =  "scaled_T_pop_" + std::to_string(i+1) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
        paramName =  "Delta*T_pop_" + std::to_string(i+1) ;
        fprintf (files.fpLikelihood->f, "%s\t", paramName.c_str());
    }
    fprintf (files.fpLikelihood->f, "\n");
}
void writeLineLikelihoodFile( int  simulationNumber, const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files ,  std::vector<Population *> &populations,
                             long double logLikCoalTree, long double logLikTrueSequences, long double logLikErrorSequences, double treeLength
                             )
{
    
    std::string paramName;
    Population *popI;
    
    fprintf (files.fpLikelihood->f, "%4d\t", simulationNumber);
    fprintf (files.fpLikelihood->f, "%2d\t", programOptions.numClones);
    fprintf (files.fpLikelihood->f, "%10.4Lf\t", logLikCoalTree);
    fprintf (files.fpLikelihood->f, "%10.4Lf\t", logLikTrueSequences);
    fprintf (files.fpLikelihood->f, "%10.4Lf\t", logLikErrorSequences);
    
    fprintf (files.fpLikelihood->f, "%10.4f\t", programOptions.fixedADOrate);
    fprintf (files.fpLikelihood->f, "%10.4f\t", programOptions.meanADOcell);
    fprintf (files.fpLikelihood->f, "%10.4f\t", programOptions.varADOcell);
    fprintf (files.fpLikelihood->f, "%10.4f\t", programOptions.meanGenotypingError);
    fprintf (files.fpLikelihood->f, "%10.4f\t", programOptions.varGenotypingError);
    fprintf (files.fpLikelihood->f, "%10.4f\t", treeLength);
    //
    fprintf (files.fpLikelihood->f, "%.4f\t", programOptions.mutationRate);
    
    for (unsigned int i = 0; i < populations.size(); ++i){
        popI =populations[i];
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->delta);
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->x);
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->deltaT);
        
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->theta);
        fprintf (files.fpLikelihood->f, "%3d\t", popI->sampleSize);
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->timeOriginInput);//multiplied by theta_i
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->timeOriginSTD);//T
        
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->scaledtimeOriginInput);//divided by theta_i
        fprintf (files.fpLikelihood->f, "%10.4Lf\t", popI->delta*popI->timeOriginSTD);//divided by theta_i
    }
    if (simulationNumber <= programOptions.numDataSets -1)
        fprintf (files.fpLikelihood->f, "\n");
}
/********************* PrepareSeparateFilesGenotypes **********************/
/* Open individual genotypes files to output results */
void PrepareSeparateFilesGenotypes(int paramSetNumber, int TreeNum,int MutationAssignNum,
                                   const FilePaths &filePaths, const ProgramOptions &programOptions,Files &files, std::vector<Population*> populations, double numMUperSite)
{
    char File[1000];
    char dir[1000];
    /* contains the simulated tree in Newick format
     */
    
    //mkdir("Results", S_IRWXU); /* Create "Results" folder (with type S_IRWXU (read, write and execute)) */
    //mkdir("Results",0);
#ifdef MAC
    strcpy (dir, ":Results:"); /* Copy the string in char variable dir = Results (char), is different mac vs windows */
#else
    strcpy (dir, "Results/");
#endif
    //strcpy (resultsDir, dir);
    std::stringstream ss;
       Population *popI;
         for( size_t i = 0 ; i < populations.size(); i++)
         {
             popI = populations[i];
             ss << std::fixed << std::setprecision(2) <<"G" << i << "=" << popI->delta << "_";
             ss << "T" << i << "=" << popI->timeOriginSTD << "_";
             ss << "n" << i << "=" << popI->sampleSize << "_";
             ss << "m" << i << "=" << popI->theta << "_";;
             ss << "x" << i << "=" << popI->x ;
             if (i!= populations.size()-1){
                 
                 ss << "_";
             }
         }
    std::string s = ss.str();
    if (programOptions.doSimulateData == YES)
    {
        /* contains SNV genotypes for every cell */
        if (programOptions.doPrintSNVgenotypes == YES)
        {
            sprintf(File,"%s/%s", filePaths.resultsDir, filePaths.SNVgenotypesDir);
            mkdir(File,S_IRWXU);
            sprintf(File,"%s/%s/%s_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVgenotypesDir, filePaths.SNVgenotypesFile, paramSetNumber+1, TreeNum+1, MutationAssignNum +1);
            
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (programOptions.fixedADOrate > 0 && programOptions.doADOcell == NO && programOptions.doADOsite == NO){
                    
                     if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVgenotypesDir, filePaths.SNVgenotypesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                     }
                     else{
                          
                         sprintf(File,"%s/%s/%s_%s_nMU=%.3f_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVgenotypesDir, filePaths.SNVgenotypesFile,s.c_str(),numMUperSite,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                     }
                    
                }
                /* variable ADO rate*/
                else{
                    
                    if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVgenotypesDir, filePaths.SNVgenotypesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    
                    }
                    else{
                        
                        sprintf(File,"%s/%s/%s_%s_nMU=%.3f_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVgenotypesDir, filePaths.SNVgenotypesFile,s.c_str(),numMUperSite,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                                           
                        
                    }
                    
                }
                
            }
            
            //sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVgenotypesDir, SNVgenotypesFile, replicate+1);
            //if ((*fpSNVgenotypes = fopen(File, "w")) == NULL)
            strcpy (files.fpSNVgenotypes->path,File);
            
            
            if (openFile(&files.fpSNVgenotypes->f, File) == -1)
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
            
            
            if (programOptions.doSimulateFromPriors==YES){
                
                
                if (programOptions.fixedADOrate > 0 && programOptions.doADOcell == NO && programOptions.doADOsite == NO){
                    
                    if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVhaplotypesDir, filePaths.SNVhaplotypesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    }
                    else{
                        
                        sprintf(File,"%s/%s/%s_%s_nMU=%.3f_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.SNVhaplotypesDir, filePaths.SNVhaplotypesFile,s.c_str(),numMUperSite,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                        
                    }
                    
                }
                /* variable ADO rate*/
                else{
                     if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir,filePaths.SNVhaplotypesDir, filePaths.SNVhaplotypesFile, populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                     }
                     else{
                         
                         sprintf(File,"%s/%s/%s_%s_nMU=%.3f_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir,filePaths.SNVhaplotypesDir, filePaths.SNVhaplotypesFile, s.c_str(),numMUperSite,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                         
                     }
                }
                
            }
            // sprintf(File,"%s/%s/%s.%04d", resultsDir, SNVhaplotypesDir, SNVhaplotypesFile, replicate+1);
            //if ((*fpSNVhaplotypes = fopen(File, "w")) == NULL)
            strcpy (files.fpSNVhaplotypes->path,File);
            if (openFile(&files.fpSNVhaplotypes->f, File) == -1)
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
            
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (populations.size()==1){
                
                sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.trueHaplotypesDir, filePaths.trueHaplotypesFile, populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                }
                else{
                    
                    sprintf(File,"%s/%s/%s_%s_nMU=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.trueHaplotypesDir, filePaths.trueHaplotypesFile, s.c_str(),numMUperSite, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    
                }
            }
            // sprintf(File,"%s/%s/%s.%04d", resultsDir, trueHaplotypesDir, trueHaplotypesFile, replicate+1);
            //if ((*fpTrueHaplotypes = fopen(File, "w")) == NULL)
            strcpy (files.fpTrueHaplotypes->path,File);
            if (openFile(&files.fpTrueHaplotypes->f, File) == -1)
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
            if (programOptions.doSimulateFromPriors==YES){
                if (populations.size()==1){
                sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.MLhaplotypesDir, filePaths.MLhaplotypesFile, populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                }
                else{
                    sprintf(File,"%s/%s/%s_%s_nMU=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.MLhaplotypesDir, filePaths.MLhaplotypesFile, s.c_str(),numMUperSite, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    
                }
            }
            
            strcpy (files.fpMLhaplotypes->path,File);
            if (openFile(&files.fpMLhaplotypes->f, File) == -1)
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
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (programOptions.fixedADOrate > 0 && programOptions.doADOcell == NO && programOptions.doADOsite == NO){
                    
                    if (populations.size()==1){
                    
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullGenotypesDir, filePaths.fullGenotypesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    }
                    else{
                        
                         sprintf(File,"%s/%s/%s_%s_nMU=%.3f_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullGenotypesDir, filePaths.fullGenotypesFile,s.c_str(),numMUperSite,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                        
                    }
                    
                }
                /* variable ADO rate*/
                else{
                    
                    if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir,filePaths.fullGenotypesDir, filePaths.fullGenotypesFile, populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    }
                    else{
                        sprintf(File,"%s/%s/%s_%s_nMU=%.3f_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir,filePaths.fullGenotypesDir, filePaths.fullGenotypesFile, s.c_str(),numMUperSite,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                        
                    }
                }
                
            }
            
            strcpy (files.fpFullGenotypes->path,File);
            if (openFile(&files.fpFullGenotypes->f, File) == -1)
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
            if (programOptions.doSimulateFromPriors==YES){
                
                
                if (programOptions.fixedADOrate > 0 && programOptions.doADOcell == NO && programOptions.doADOsite == NO){
                    
                    if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullHaplotypesDir, filePaths.fullHaplotypesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    }
                    else{
                        
                        sprintf(File,"%s/%s/%s_%s_nMU=%.3f_fADO=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullHaplotypesDir, filePaths.fullHaplotypesFile,s.c_str(),numMUperSite,programOptions.fixedADOrate,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    }
                    
                }
                /* variable ADO rate*/
                else{
                    
                    if (populations.size()==1){
                    sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullHaplotypesDir, filePaths.fullHaplotypesFile,populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    
                    }
                    else{
                        sprintf(File,"%s/%s/%s_%s_nMU=%.3f_mADOc=%.3f_vADOc=%.3f_mSeqE=%.3f_vSeqE=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.fullHaplotypesDir, filePaths.fullHaplotypesFile,s.c_str(),numMUperSite,programOptions.meanADOcell,programOptions.varADOcell,programOptions.meanGenotypingError,programOptions.varGenotypingError, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                                          
                        
                    }
                }
            }
            
            strcpy (files.fpFullHaplotypes->path,File);
            if (openFile(&files.fpFullHaplotypes->f, File) == -1)
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
            
            if (programOptions.doSimulateFromPriors==YES){
                
                if (populations.size()==1){
                
                sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.VCFdir, filePaths.VCFfile, populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                }
                else{
                    sprintf(File,"%s/%s/%s_%s_nMU=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.VCFdir, filePaths.VCFfile, s.c_str(),numMUperSite, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                    
                }
            }
            
            strcpy (files.fpVCF->path,File);
            if (openFile(&files.fpVCF->f, File) == -1)
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
            
            
            if (programOptions.doSimulateFromPriors==YES){
                
               if (populations.size()==1){
                sprintf(File,"%s/%s/%s_Delta=%.3Lf_T=%.3Lf_theta=%.3Lf_nMU=%.3f_n=%d_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.CATGdir, filePaths.CATGfile, populations[0]->delta,populations[0]->timeOriginSTD, populations[0]->theta,numMUperSite,populations[0]->sampleSize, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
               }
               else{
                   sprintf(File,"%s/%s/%s_%s_nMU=%.3f_%04d_%04d_%04d.txt", filePaths.resultsDir, filePaths.CATGdir, filePaths.CATGfile, s.c_str(),numMUperSite, paramSetNumber+1,  TreeNum+1, MutationAssignNum +1);
                   
               }
            }
            
            
            strcpy (files.fpCATG->path,File);
            if (openFile(&files.fpCATG->f, File) == -1)
            {
                fprintf (stderr, "Can't open \"%s\"\n", File);
                exit(-1);
            }
        }
    }
}
/***************************** InitPopulationsCoalescentEvents*******************************/
/* InitPopulationsCoalescentEvents*/
void InitPopulationsCoalescentEvents( int numClones,  std::vector<Population*> &populations) {
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
void InitializeGenomes (TreeNode *p, long int *seed,  int alphabet, int doUserGenome, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures, double cumfreq[4],long double *triNucFreq, char **cellNames,
                        const gsl_rng *randomGsl, boost::mt19937* rngBoost)
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
                        ran =  Random::randomUniformFromGsl2(randomGsl);
                        for (i=0; i<4; i++)
                        {
                            if (ran <= cumfreq[i])
                            {
                                //p->maternalSequence[site]=p->paternalSequence[site]=i;
                                p->maternalSequence.at(site)=i;
                                p->paternalSequence.at(site)=i;
                                
                                allSites[site].referenceAllele = p->maternalSequence.at(site);//   // then allSites[site].referenceAllele hosts the reference genome
                                break;
                            }
                        }
                    }
                }
                else /* initialize genome with trinucleotide frequencies */
                {
                    SimulateTriNucFreqGenome (cell, seed, p, alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, randomGsl,  rngBoost);
                }
            }
            else
            {
                anccell = p->anc1->label;
                for (site=0; site<numSites; site++)
                {
                    //                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                    //                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                    
                    p->maternalSequence.at(site)= p->anc1->maternalSequence.at(site);
                    p->paternalSequence.at(site)= p->anc1->paternalSequence.at(site);
                }
            }
        }
        else{
            //            for (site=0; site<numSites; site++){
            //                p->maternalSequence.at(site)=0;
            //                p->paternalSequence.at(site)=0;
            //                p->numbersMutationsUnderSubtreePerSite.at(site)=0;
            //                p->numbersMaternalMutationsPerSite.at(site)=0;
            //                p->numbersPaternalMutationsPerSite.at(site)=0;
            //
            //            }
        }
        InitializeGenomes (p->left, seed,   alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, cellNames, randomGsl, rngBoost);
        InitializeGenomes (p->right, seed,  alphabet,  doUserGenome,  numSites, allSites,  doGeneticSignatures,  cumfreq,triNucFreq, cellNames, randomGsl, rngBoost);
    }
}


/* Returns the sum of the branch lengths for a given tree */
double SumBranches (TreeNode *p, double mutationRate, std::string &healthyTipLabel)
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
        if (strlen(p->left->cellName) > 0 && std::strcmp(p->left->cellName, healthyTipLabel.c_str())==0){
            SumBranches (p->left,  mutationRate, healthyTipLabel);
        }
        if (strlen(p->left->cellName) > 0 && std::strcmp(p->right->cellName, healthyTipLabel.c_str())==0){
            SumBranches (p->right,   mutationRate, healthyTipLabel);
        }
        
    }
    
    return sum;
}
/************************* MakeCoalescenceTree2 ************************/
/* Builds a genealogy under the structured  ' */ /* this function go by events CA, MIG, CONV */

TreeNode *MakeCoalescenceTree2 (long int *seed, std::vector<Population *> &populations,
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
                                std::vector<TreeNode *> &nodes, // empty vector
                                std::vector<TreeNode *> &treeTips,
                                TreeNode    *treeRootInit,
                                long double K,
                                const gsl_rng* randomGsl,
                                boost::mt19937* rngBoost
                                ) {
    
    int      c,  i,  m, cumIndivid, isCoalescence, whichInd,
    newInd, eventNum, numActiveGametes,
    isMigration, whichClone, currentNumberAliveClones;
    int     labelNodes;
    double    currentTime, eventTime;
    TreeNode  *p;
    int     *numParcialActiveGametes;
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
    
    nodes.clear();
    for (i = 0; i < numNodes; i++)
    {
        p = new TreeNode(programOptions.numSites);
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
                           eventNum,
                           K,
                           randomGsl,
                           rngBoost);
        if (i< numClones-1)   //if it is not the last one
        {
            //choose the father population from which the population i came
            
            fatherPop= ChooseFatherPopulation(populations, numClones, currentPop, seed,  programOptions.noisy, programOptions.K, randomGsl);
            currentPop->FatherPop = fatherPop;
            //update list of migrant times
            Population::UpdateListMigrants(numClones, currentPop, fatherPop);
        }
        i = i + 1;
    }
    //BuildTree also updates the length field of nodes
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
                               labelNodes,
                               randomGsl
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
void SimulatePopulation( Population *popI, std::vector<Population*> &populations,
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
                        std::vector<TreeNode *> &nodes,
                        int &nextAvailable,
                        int &numActiveGametes,
                        int &labelNodes,
                        double &currentTime,
                        int &eventNum,
                        long double K,
                        const gsl_rng* randomGsl,
                        boost::mt19937* rngBoost
                        )
{
    int   i, k, isCoalescence,
    firstInd, newInd,
    isMigration, whichClone;
    double     eventTime;
    TreeNode  *p,  *r;
    
    eventNum = 0;
    
    double ThisRateCA = 0.0;
    double ThisTimeCA_W = 0.0;
    double ThisTimeCA_V1 = 0.0;
    double ThisTimeCA_V2 = 0.0;
    int numParcialActiveGametes;
    
    numParcialActiveGametes = popI->numActiveGametes;
    
    int numMigrations = (popI->numIncomingMigrations); //taking into account also the    time of origin as a migration
    
    double timeNextMigration;
    int indexNextMigration = 0;
    Population *incomingPop;
    
    if (programOptions.noisy > 1)
        std::cout << "\n\n>> Simulating evolutionary history of clone " << popI->index << " (number active gametes " << popI->numActiveGametes << " , original time to origin "<< popI->timeOriginInput << std::endl;
    if (programOptions.noisy > 1)
        std::cout << "\n\n> Simulating evolutionary history of clone  or order  "<< popI->order<< std::endl;
    currentTime=0;
    while (indexNextMigration < numMigrations) {
        timeNextMigration = popI->immigrantsPopOrderedByModelTime[indexNextMigration].first;
        if ( popI->numActiveGametes >= 2) {
            ThisRateCA = (double)  popI->numActiveGametes * ((double)  popI->numActiveGametes - 1) / 2.0;
            ThisTimeCA_W = Random::RandomExponential (ThisRateCA, seed, true, randomGsl,rngBoost) ;//this is Kingman coal time
            ThisTimeCA_V1 = Population::FmodelTstandard (currentTime , popI->timeOriginSTD, popI->delta,   K);//this is current time in Kingman coal time
            
            //ThisTimeCA_V2 = Population::GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD, popI->delta,  K);//this is current time in Kingman coal time
            
            //assert(ThisTimeCA_V1==ThisTimeCA_V2);
            
            ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;//this is new  time in Kingman coal time
            
            ThisTimeCA_V2 = Population::GstandardTmodel(ThisTimeCA_V1, popI->timeOriginSTD , popI->delta,  K);//this is new  time in Kingman coal time
        }
        else
        {
            ThisTimeCA_V2 = timeNextMigration + 1.0; // it reached a "provisional" MRCA
            // fprintf (stderr, "\n Only 1 active gamete and %d true migrations out of %d true migrations \n", indexNextMigration, numMigrations - 1);
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
                std::cout << "\n\n*** Event " << eventNum << "  currentTime (model time) = " << ThisTimeCA_V2 <<", currentTime (standard time) = " << ThisTimeCA_V1 << std::endl;
                std::cout << "\n\n*** Event " << eventNum <<  "  currentTime (input units) = " << ThisTimeCA_V2<< std::endl;
            }
            if (programOptions.noisy == 4)
                std::cout <<"* Coalescence *\n"<< std::endl;
            MakeCoalescenceEvent(populations, popI, nodes, numClones, seed, programOptions.noisy, numActiveGametes, nextAvailable, labelNodes, currentTime,  programOptions.numNodes, randomGsl);
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
                    std::cout << "\n\n*** Event "<< eventNum << "  currentTime (model units) = " <<  ThisTimeCA_V2<< std::endl;
                }
                if (programOptions.noisy == 4)
                {
                    std::cout << "* Migration *"<< std::endl;
                    
                }
                
                incomingPop = popI->immigrantsPopOrderedByModelTime[indexNextMigration].second;
                
                //r->indexOldClone = incommingPop->index;
                p = nodes[incomingPop->nodeIdAncestorMRCA];
                
                if (programOptions.noisy >1)
                    std::cout << "\n The incoming population "<< incomingPop->order << " to  population "<< popI->order << " with node "<< p->index<< " and time " <<  p->timePUnits << std::endl;
                
                //p = incomingPop->MRCA;
                //p = *nodes + (incommingPop->nodeIdAncesterMRCA); // root of younger clone
                indexNextMigration = indexNextMigration + 1;
                p->indexCurrentClone = popI->index;
                p->indexOldClone = incomingPop->index;
                p->orderCurrentClone = popI->order;
                
                
                k = p->indexCurrentClone;
                incomingPop->numActiveGametes = incomingPop->numActiveGametes - 1; /* now the other clone has 1 less node */
                
                popI->idsActiveGametes[popI->numActiveGametes]=p->index;//adding the superfluos node
                popI->numActiveGametes = popI->numActiveGametes + 1; /* now this clone has 1 more node */
                
                if (programOptions.noisy > 1)
                    std::cout << "\t|\tCurrentTime (input units) = "<< p->timePUnits<< std::endl;
                
                /* memory for number of nodes */
                if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                {
                    /* ReallocNodes(&numNodes, activeGametes); */
                    if (programOptions.noisy == 4)
                        std::cout << "\n\n...Doing reallocation of nodes (Coalescence)"<< std::endl;
                    numNodes += INCREMENT_NODES;
                    /* realloc */
                    
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
                    std::cout <<"\n\n*** Clone origin ***\n"<< std::endl;
                if (programOptions.noisy == 4)
                    std::cout << "Clone origin "<< popI->index <<" at time (model units) = " << currentTime<< std::endl;
                if (popI->order < numClones - 1) // do not do it for the last clone
                {
                    newInd = nextAvailable;
                    r = nodes[newInd];    /* new ancestor */
                    r->index = nextAvailable;
                    r->label = labelNodes;
                    labelNodes=labelNodes+1;
                    r->indexOldClone =r->indexCurrentClone = popI->index;
                    r->orderCurrentClone = popI->order;
                    
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
                    
                   // if (popI->order == numClones-1)//oldest population
                    //{
                   //     r->timePUnits = currentTime ;
                   // }
                    
                   // else{
                        r->timePUnits = currentTime * popI->x;
                        
                    //}
                    popI->MRCA = p;
                    
                    //connectNodes(p, NULL, r);
                    //fprintf (stderr, "\n r->index = %d, r->time = %lf\n", r->index, r->time);
                    /* readjust active nodes */
                    nextAvailable=nextAvailable+1; /* 1 node more is available */
                    popI->idsActiveGametes[0] = newInd;//always will be in the 0  position because there is only one left
                    
                    //popI->idsGametes[popI->numGametes] = newInd; r is a superflous node and it will be removed so no need to add it
                    //popI->numGametes = popI->numGametes +1;
                    
                    if (programOptions.noisy > 1)
                        std::cout << "Creating origin node, it creates node "<< newInd<< " derived from node " << firstInd<< std::endl;
                    if (programOptions.noisy > 1)
                        std::cout << "\t|\tCurrentTime (input units) = "<< r->timePUnits<< std::endl;
                    /* memory for number of nodes */
                    if (nextAvailable >= numNodes)  /* if there aren't enough nodes it go into and it addition more */
                    {
                        /* ReallocNodes(&numNodes); */
                        if (programOptions.noisy == 4)
                            std::cout << "\n\n...Doing reallocation of nodes (Coalescence)\n"<< std::endl;
                        numNodes += INCREMENT_NODES;
                        /* realloc */
                        
                    }
                }
                else  {//origin of oldest pop reached
                    popI->nodeIdAncestorMRCA=nextAvailable-1;//for the last population, nodeIdAncesterMRCA is the MRCA instead of ancester of MRCA
                    r = nodes[nextAvailable-1];//popI->idsActiveGametes[0]
                    r->indexOldClone = r->indexCurrentClone = popI->index;
                    r->orderCurrentClone = popI->order;
                    popI->MRCA= r;
                    //fprintf (stderr, "\n origin of the oldest population  %d", popI->nodeIdAncestorMRCA);
                }
            }
        }
        if (programOptions.noisy > 3)
        {
            std::cout << "\nActive nodes: " <<  popI->numActiveGametes << std::endl;
            for (i = 0; i < popI->numActiveGametes; i++)
                std::cout << " %d" << popI->idsActiveGametes[i] << std::endl;
            std::cout << "\t|\tNext node available = "<< nextAvailable<< std::endl;
        }
    }
    if (programOptions.noisy > 1)
        std::cout << "\n\nEvolutionary history of clone %d is completed " <<  popI->index << std::endl;
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
    if (p != NULL)
    {
        if (p->isOutgroup == YES)     /* Outgroup */
        {
            strcpy( p->cellName,"healthycell");
            
            
            if (asprintf(&newickString,  "healthycell:%10.9Lf",  (p->anc1->timePUnits - p->timePUnits) * mutationRate)<0)
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
                if (asprintf(&newickString,   "%s:%10.9Lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate)<0)
                    return NULL;
                //snprintf(newickString,  size,  "%s:%10.9lf",  p->observedCellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate);
                return newickString;
            }
            else{
                if (asprintf(&newickString,   "%s:%10.9Lf",  p->cellName, (p->anc1->timePUnits - p->timePUnits)*mutationRate)<0)
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
                
                if (asprintf(&newickString, "(%s,%s):%10.9Lf", left, right,  (p->anc1->timePUnits - p->timePUnits)*mutationRate )<0)
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
double ProbabilityCloneiFromClonej2 (Population *PopI, Population* PopJ, std::vector<Population*> &populations, int numClones, long double K)
{
    double  ProbabilityIJ, AboveTerm, BelowTerm;
    int     l, j;
    double  logh, h, a, b, c, d, e, t, cum;
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
    
    t = (PopI->timeOriginSTD ) * (PopI->x) / ( PopJ->x);
    logh = Population::LogCalculateH(t, PopJ->timeOriginSTD, PopJ->delta, K );
    
    //AboveTerm = ( PopJ->popSize) * h;
    h = exp(logh);
    AboveTerm = (PopJ->x) * h;
    j=0;
    for (l = PopI->order + 1; l < numClones; l++)
    {    p = populations[l];
        t = (PopI->timeOriginSTD ) * (PopI->x) / ( p->x);
        logh = Population::LogCalculateH(t, p->timeOriginSTD, p->delta, K);
        //cum = cum + ( ( p->popSize) * h);
        cum = cum +  (p->x)*exp(logh);
    }
    
    BelowTerm = cum;
    assert(BelowTerm>0);
    ProbabilityIJ = AboveTerm / BelowTerm;
    return ProbabilityIJ;
}


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
    programOptions.alphaSites = infinity2;          /* alpha shape of the gamma distribution for rate variation among sites */
    programOptions.alphaBranches = infinity2;       /* alpha shape of the gamma distribution for rate variation among lineages */
    programOptions.alphaCoverage = infinity2;       /* alpha shape of the gamma distribution for coverage */
    programOptions.doSimulateFixedNumMutations = NO;    /* whether to simulate a fixed number of mutations */
    programOptions.doUserTree = NO;                /* whether to assume a user tree instead od making the coalescent */
    programOptions.doUserGenome = NO;                /* whether to use a user genome instead of a simulated one */
    programOptions.doPrintSNVgenotypes = NO;        /* whether to print SNVs */
    programOptions.doPrintSNVhaplotypes = NO;      /* whether to print haplotypes */
    programOptions.doPrintTrueHaplotypes = NO;      /* whether to print haplotypes without errors */
    programOptions.doPrintMLhaplotypes = NO;          /* whether to print ML haplotypes */
    programOptions.doPrintFullHaplotypes = NO;        /* whether to print sequences */
    programOptions.doPrintFullGenotypes = YES;        /* whether to print all genotypes (variable + invariable) */
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
    
    programOptions.K=0.8;
    programOptions.minSampleSize = 10;
    programOptions.maxSampleSize = 100;
}
/***************************** ReadMCMCParametersFromFile *******************************/
/* Reads parameter values from the parameter file */

void ReadMCMCParametersFromFile(ProgramOptions &programOptions, FilePaths &filePaths, MCMCOptions &mcmcOptions)

{
    int   j;
    char  ch;
    float   argument;
    
    
    int argumentInt;
    
    std::string argumentString;
    
    /* Used:  */
    
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
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad number of chains (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                mcmcOptions.numChains =argumentInt;
                break;
            case 'I':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad number of iterations (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                mcmcOptions.Niterations =argumentInt;
                break;
            case 'S':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad number of thinni g iterations (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                mcmcOptions.thinning =argumentInt;
                break;
            case 'R':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad max number of proposal before rejection  (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                mcmcOptions.maxNumberProposalAttempts =argumentInt;
                break;
            case 'F':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0 || argumentInt > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad option user fixed tree  (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                programOptions.doUseFixedTree =argumentInt;
                break;
            case 'Z':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad  number of clones  (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                programOptions.numClones =argumentInt;
                break;
            case 'C':
                if (fscanf(stdin, "%d", &argumentInt) != 1 || argumentInt < 0 || argumentInt > 1)
                {
                    fprintf(stderr, "PARAMETER ERROR: Bad option number of clones known  (%d)\n\n", (int)argumentInt);
                    Output::PrintUsage();
                }
                programOptions.numberClonesKnown =argumentInt;
                break;
            case 'T':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths.inputTreeFile, "trees.tre");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths.inputTreeFile[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths.inputTreeFile[j] = '\0';
                }
                break;
            case 'G':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths.inputGenotypeFileFasta, "geno.fasta");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths.inputGenotypeFileFasta[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths.inputGenotypeFileFasta[j] = '\0';
                }
                break;
            case 'P':
                ch = fgetc(stdin);
                if (isspace(ch))
                {
                    strcpy(filePaths.inputGenotypeFilePhylip, "geno.phylip");
                }
                else
                {
                    j = 0;
                    do
                    {
                        filePaths.inputGenotypeFilePhylip[j] = ch;
                        j++;
                        ch = fgetc(stdin);
                    }
                    while (!isspace(ch));
                    filePaths.inputGenotypeFilePhylip[j] = '\0';
                }
                break;
            default :
                fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
                Output::PrintUsage();
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
/************************ computeUnfoldedISMSFS ***********************/
/* compute unfolded version of SFS */
void computeUnfoldedISMSFS(int numSites,std::vector<SiteStr> &allSites,int numSNVs, std::vector<int> &SNVsites, std::vector<int> &SFS, std::vector<int> &numberDifferences){
    int numberDiff;
    int site;
    int posSNV;
    for (site=0; site<numSites; site++){
        SFS[site]=0;
    }
    for (site=0; site<numSNVs; site++){
        // if (allSites[site].isVariant == YES){
        posSNV =SNVsites[site];
        numberDiff=allSites[posSNV].numberDiffReference;
        SFS[numberDiff]=SFS[numberDiff]+1;
        numberDifferences[site]=numberDiff;
        //}
    }
}
int countTrueVariants (std::vector<TreeNode *> &nodes,  int numSites, int numCells, TreeNode *HEALTHY_ROOT, std::vector<SiteStr> &allSites, std::vector<int> &variantSites, std::vector<int> &SNVsites )
{
    int        cell, site;
    int        nVariants = 0;
    int numberDiff=0;
    int isVariant=NO;
    TreeNode* p;
    for (site=0; site<numSites; site++)
    {   isVariant=NO;
        numberDiff=0;
        for (cell=0; cell<numCells; cell++)
        {    p= nodes.at(cell);
            if (p->maternalSequence[site] != Definitions::DELETION && p->maternalSequence[site] != HEALTHY_ROOT->maternalSequence[site])
            {
                isVariant=YES;
                numberDiff++;
            }
            if (p->paternalSequence[site] != Definitions::DELETION && p->paternalSequence[site] != HEALTHY_ROOT->paternalSequence[site])
            {
                isVariant=YES;
                numberDiff++;
            }
        }
        if (isVariant){
            allSites[site].isVariant = YES;
            SNVsites[ nVariants] = site;
            nVariants++;
        }
        allSites[site].numberDiffReference=numberDiff;
    }
    return nVariants;
}
void computeStatisticsNumberMutations(std::vector<SiteStr> allSites, long double &maternalMutationPerSite, long double &paternalMutationPerSite){
    
    
    long double maternalTotal=0.0;
    long double paternalTotal=0.0;
    
    for(size_t i =0; i< allSites.size();++i){
        
        maternalTotal+=allSites[i].numMutationsMaternal;
        paternalTotal+=allSites[i].numMutationsPaternal;
        
    }
    maternalMutationPerSite = maternalTotal;
    paternalMutationPerSite = paternalTotal;
}
long double computeParamPowerDistribQuantileUntil(long double areaUntilb, long double b, long double from)
{
    
    if (areaUntilb >=0 && areaUntilb <=1)
    {
        long double result=areaUntilb/ ((1.0 -areaUntilb)*(b-from) );
        return result;
    }
    else
        return 0;
}
void setDefaultOptions(ProgramOptions &programOptions, MCMCOptions &mcmcOptions ){
    //programOptions
    //programOptions.K=0.8;
    programOptions.K=0.0;
    programOptions.K_inference =2.0;
    programOptions.numberClonesKnown=YES;
    programOptions.populationSampleSizesKnown = YES;
    programOptions.doPrintTimes = YES;
    
    programOptions.doUseObservedCellNames =1;
    programOptions.doUseGenotypes = YES;
    programOptions.doUseFixedTree =NO;
    programOptions.seed = 1248697;
    programOptions.thereisOutgroup = YES;
    programOptions.outgroupSelection =1;
    
    programOptions.mutationRate= 0.01;
    programOptions.doUsefixedMutationRate=0;
    programOptions.doSimulateData=NO;
    
    programOptions.doUseFixedTree = NO;
    
    mcmcOptions.slidingWindowSizeTotalEffectPopSize = 20000;
    mcmcOptions.slidingWindowSizeGrowtRate=0.01;
    mcmcOptions.tuningParameter = 1;
    
    mcmcOptions.totalEffectPopSizefrom = 7;
    mcmcOptions.totalEffectPopSizeto = 11;
    //mcmcOptions.MutRatefrom = -20;
    //mcmcOptions.MutRateto = -10;
    mcmcOptions.MutRatefrom = -4;
    mcmcOptions.MutRateto = 1;
    mcmcOptions.Deltafrom = -4;
    mcmcOptions.Deltato = 8;
    mcmcOptions.fixedLambda=1;
    mcmcOptions.priorsType =1; //0: log uniform, 1: exponential priors, 2: power law
    mcmcOptions.kernelType = 0;//0: multiplier move, 1: normal
    
    
    mcmcOptions.sigmaNormalKernelTotalEffectivePopulationSize = 300;
    mcmcOptions.sigmaNormalKernelMutationRate =  0.0000000001 ;
    mcmcOptions.sigmaNormalKernelGrowthRate=0.00001;
    //mcmcOptions.lambdaExponentialPriorMutationRate=10000000;
    
    mcmcOptions.lambdaExponentialMutationRateSimulation=100;
    
    mcmcOptions.lambdaExponentialPriorMutationRate=1;
    mcmcOptions.lambdaExponentialPriorSeqError = 10;
    mcmcOptions.lambdaExponentialPriorDropoutError = 10;
    //mcmcOptions.lambdaExponentialPriorMutationRate=0.0000001;
    //mcmcOptions.lambdaExponentialPriorTotalEffectivePopSize=0.0001;
    
    mcmcOptions.lambdaExponentialPriorGrowthRate = 0.1;
    mcmcOptions.lambdaExponentialGrowthRateSimulation = 0.01;
    //mcmcOptions.lambdaExponentialPriorGrowthRate= 1;
    // mcmcOptions.lambdaExponentialPriorGrowthRate= 0.00001;
    mcmcOptions.GrowthRatefrom = -1;
    mcmcOptions.GrowthRateto = 6;
    mcmcOptions.parameterPowerLawDistributionTotalEffectPopSize=0.9;
    mcmcOptions.parameterPowerLawDistributionMutationRate= 5000;
    mcmcOptions.parameterPowerLawDistributionGrowthRate=100;
    mcmcOptions.parameterPowerLawDistributionTimeOriginInputOldestPop=1000;
    mcmcOptions.sigmaNormalKernelTimeofOrigin=0.001;
    mcmcOptions.lengthIntervalMultiplier = 0.9;
    //mcmcOptions.lengthIntervalMultiplierTheta = 1.6;
    //mcmcOptions.lengthIntervalMultiplierDeltaT = 1.9;
    mcmcOptions.lengthIntervalMultiplierTheta = 1.2;
    mcmcOptions.lengthIntervalMultiplierDeltaT = 1.1;
    
    programOptions.seqErrorRate=programOptions.sequencingError=0;
    programOptions.dropoutRate=programOptions.ADOrate=0;
    
    mcmcOptions.paramMultiplierMoveTheta = 2;
    mcmcOptions.paramMultiplierEffectPopSize = 2;
    //mcmcOptions.Niterations = 10000000;
    //mcmcOptions.numberWarmUpIterations = mcmcOptions.Niterations / 2.0;
    mcmcOptions.numberWarmUpIterations = 0;
    mcmcOptions.useSequencesLikelihood =0;
    
    mcmcOptions.verbose = 0;
    
    mcmcOptions.paramMultiplierGrowthRate = Utils::parameterMultiplierMCMCmove(mcmcOptions.lengthIntervalMultiplierDeltaT);
    mcmcOptions.paramMultiplierTheta =Utils::parameterMultiplierMCMCmove (mcmcOptions.lengthIntervalMultiplierTheta);
    
    // mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop = 0.008;
    //mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop = 1;
    mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop = 0.9;
    //mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop = 0.5;
    mcmcOptions.paramMultiplierTimeOriginOldestPop =Utils::parameterMultiplierMCMCmove (mcmcOptions.lengthIntervalMultiplierTimeOriginOldestPop);
    
    mcmcOptions.upperBoundTimeOriginInputOldestPop = 7;
    mcmcOptions.percentIterationsToComputeThinnig=0.1;
    mcmcOptions.doThinning=true;
    // mcmcOptions.splitThetaDeltaTmoves= true;
    mcmcOptions.thresholdAccceptanteRate=0.44;
    mcmcOptions.updateLengthMultiplierMCMCMove=0.2;
    
    mcmcOptions.thresholdAutoCorrelation=0.05;
    mcmcOptions.doThinning=true;
    mcmcOptions.numberChainsPerTree = 4;
}
void printProgramHeader()
{
    //auto end = std::chrono::system_clock::now();
    std::time_t timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "+                                                                                      " << std::endl;
    std::cout << "\n\n+Bayesian inference of tumour single cell growth rates using the structured coalescent " << std::endl;
    std::cout << "+                                                                                      " << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "+                         Version 1.0.0    [August-1-2020]                             " << std::endl;
    std::cout << "+                       Program started at  "   <<  std::ctime(&timeNow)                 << std::endl;
    
    
}
void printSimulatorProgramHeader()
{
    //auto end = std::chrono::system_clock::now();
    std::time_t timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "+                                                                                      " << std::endl;
    std::cout << "+Simulator of single cell DNA(sc-DNA) genotypes under the Birth-Death structured coalescent " << std::endl;
    std::cout << "+                                                                                      " << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "+                         Version 1.0.0    [August-1-2020]                             " << std::endl;
    std::cout << "+                       Program started at  "   <<  std::ctime(&timeNow)                 << std::endl;
    
    
}
void simulateTrees(int numberTrees,std::vector<StructuredCoalescentTree *> &structuredCoalTrees,  std::vector<pll_rtree_t *> &trees,          std::vector<long double> &realThetas,
                   std::vector<std::vector<long double>> &realDeltaTs,
                   std::vector<std::vector<long double>> &realTs,
                   std::vector<int> & sampleSizes, ProgramOptions &programOptions, MCMCOptions & mcmcOptions, std::vector<gsl_rng * > rngGsl, std::vector<boost::mt19937* > rngBoost, std::vector<std::vector<int> > &ObservedData,char* ObservedCellNames[], pll_msa_t *msa, std::string& healthyTipLabel, long double seqErrorRate,
                   long double dropoutRate  )
{
    
    int numNodes;
    std::vector<pll_rnode_t*> rnodes;
    std::vector<pll_rnode_t*> rtreeTips;
    std::vector<pll_tree_edge_t *> edges;
    std::vector<std::pair<double, pll_tree_edge_t *> > edgeLengths;
    if (programOptions.numberClonesKnown)
    {
        numNodes = programOptions.numNodes;
    }
    for (unsigned int i=0; i < numberTrees ; i++)
    {
        auto theta =sampleMutationRateSimulation( mcmcOptions, programOptions, rngGsl.at(i), rngBoost.at(i));
        
        realThetas.push_back(theta);
        
        auto structuredCoalTree = new StructuredCoalescentTree(programOptions.numClones, sampleSizes, theta, mcmcOptions, programOptions,   rngGsl.at(i),rngBoost.at(i), ObservedData, ObservedCellNames, msa, healthyTipLabel,  seqErrorRate,
                                                               dropoutRate );
        
        realDeltaTs.push_back(structuredCoalTree->getDeltaTs());
        realTs.push_back(structuredCoalTree->getTs());
        trees.push_back(structuredCoalTree->getTree());
        structuredCoalTrees.push_back(structuredCoalTree);
    }
}
long double  initMutationRate( MCMCOptions &mcmcOptions, ProgramOptions &programOptions, const gsl_rng * randomGsl, boost::mt19937* rngBoost) {
    long double theta=0.0;
    if (programOptions.doUsefixedMutationRate)
        theta = programOptions.mutationRate;
    else
    {
        if(mcmcOptions.priorsType==0)
        {
            theta = Random::RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, mcmcOptions.useGSLRandomGenerator, randomGsl, rngBoost);
        }
        else if(mcmcOptions.priorsType==1)
        {
            if (mcmcOptions.fixedValuesForSimulation)
                theta= 1.0 / mcmcOptions.lambdaExponentialPriorMutationRate;
            else
                theta = Random::RandomExponentialStartingFrom( mcmcOptions.lambdaExponentialPriorMutationRate, 0, mcmcOptions.useGSLRandomGenerator, randomGsl, rngBoost);
        }
        
        else{
            
            theta = Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionMutationRate, 0, mcmcOptions.useGSLRandomGenerator, randomGsl,  rngBoost);
        }
        std::cout << "True scaled mutation rate outside: " << theta << std::endl;
    }
    return theta;
}
long double  sampleMutationRateSimulation( MCMCOptions &mcmcOptions, ProgramOptions &programOptions,const gsl_rng *randomGsl, boost::mt19937* rngBoost) {
    long double theta=0.0;
    if (programOptions.doUsefixedMutationRate)
        theta = programOptions.mutationRate;
    else
    {
        if(mcmcOptions.priorsType==0)
        {
            theta = Random::RandomLogUniform(mcmcOptions.MutRatefrom, mcmcOptions.MutRateto, mcmcOptions.useGSLRandomGenerator, randomGsl, rngBoost);
            
        }
        else if(mcmcOptions.priorsType==1)
        {
            if (mcmcOptions.fixedValuesForSimulation)
                theta= 1.0 / mcmcOptions.lambdaExponentialMutationRateSimulation;
            else
                theta = Random::RandomExponentialStartingFrom( mcmcOptions.lambdaExponentialMutationRateSimulation, 0, mcmcOptions.useGSLRandomGenerator, randomGsl, rngBoost);
        }
        else{
            
            theta = Random::RandomPowerLawDistribution(mcmcOptions.parameterPowerLawDistributionMutationRate, 0, mcmcOptions.useGSLRandomGenerator, randomGsl,  rngBoost);
        }
        std::cout << "True scaled mutation rate outside: " << theta  << std::endl;
    }
    return theta;
}


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

//  mutationModel.cpp
//  simul
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//


#include "mutationModel.h"

#include "eigen.hpp"

#include "random.h"

#include "utils.hpp"

using namespace Definitions;

void SimulateISM (TreeNode *treeRoot, int genome, int doISMhaploid, long int *seed,  std::vector<int> &DefaultModelSites, int numDefaultModelSites, std::vector<int> &AltModelSites, int numAltModelSites, long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4], long double  mutationRate, const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
    int   i, trials, numMutations, mutationsSoFar;
    long double   totalBranchSum;
    long double  cumBranchLength =0;
    long double  uniform =0;
    int mutationAdded;
    long double  ran=0;
    int numModelSites = numDefaultModelSites;
    
    totalBranchSum = totalTreeLength * numModelSites * mutationRate;
    
    int ploidyFactor;
    
    if (doISMhaploid==NO){
        ploidyFactor=1;
    }
    else{
        ploidyFactor=2;
    }
    
    do{
        numMutations = Random::RandomPoisson (totalBranchSum, seed, randomGsl,  rngBoost );
    }while( numISMmutations + numMutations > (ploidyFactor * numModelSites));
    
    
    numISMmutations = numISMmutations + numMutations;
    trials = 0;
    mutationsSoFar = 0;
    mutationAdded=NO;
    while (mutationsSoFar < numMutations)
    {
        i = Random::RandomUniformTo(numModelSites, seed, true, randomGsl,  rngBoost);
        
        while (((doISMhaploid == NO)  && (allSites[i].numMutations != 0))  ||
               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[i].numMutationsMaternal != 0)) ||
               ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[i].numMutationsPaternal != 0)))
        {
            //            fprintf (stderr, "\n\n site %d again!", i);
            i = Random::RandomUniformTo(numModelSites, seed, true, randomGsl,  rngBoost);
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
            SimulateISMDNAforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij, mutationRate, uniform, cumBranchLength,  ran, randomGsl,  rngBoost);
        else{
            mutationAdded=NO;
            SimulateISMforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij,mutationRate ,  cumBranchLength,  uniform,  mutationAdded, randomGsl,  rngBoost);
        }
        mutationsSoFar++;

        std::cout<< "\nmutations = " << numMutations << " mutations so far = " << mutationsSoFar << std::endl;
        std::cout<<"\n position = "<< i<<std::endl;
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
    }
}
void SimulateGenotype (TreeNode *treeRoot, std::vector<int> numberOfSitesWithKMutations, int numberVariableSites, int doISMhaploid, long int *seed,  std::vector<int> &DefaultModelSites, int numDefaultModelSites, std::vector<int> &AltModelSites, int numAltModelSites, long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4], long double  mutationRate, const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
   long double  cumBranchLength = 0;

   long double  uniform =0;
   long double  ran=0;
   int cumulativeNumberSites = 0;
    
   for (unsigned int i=1; i< numberOfSitesWithKMutations.size(); i++)
    {
        for (unsigned int j=cumulativeNumberSites; j< (cumulativeNumberSites+numberOfSitesWithKMutations[i]); j++)
        {
            for (unsigned int k=1; k<= i; k++)
              SimulateISMGenotypeforSite(treeRoot, 1, j, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij, mutationRate, uniform, cumBranchLength,  ran, randomGsl,  rngBoost);
        }
        cumulativeNumberSites+=numberOfSitesWithKMutations[i];
    }
    
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

void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4],long double  mutationRate, long double  &uniform, long double  &cumBranchLength, long double  &ran, const gsl_rng *randomGsl, boost::mt19937* rngBoost )
{
    
    //    static long double     cumBranchLength, uniform, ran;
    
    int             j, cell, anccell, ancstate;
    
    
    
    if (p != NULL)
        
    {
        
        cell = p->label;
        
        
        
        if ( p->anc1 == NULL)
            
        {
            
            cumBranchLength = 0;
            
            uniform = Random::randomUniformFromGsl2(randomGsl) * totalTreeLength;
            
        }
        
        else
            
        {
            
            anccell = p->anc1->label;
            
            if(genome == MATERNAL )
                
                ancstate = p->anc1->maternalSequence[site];
            
            else
                
                ancstate = p->anc1->paternalSequence[site];
            
            //            cumBranchLength += p->length;// ->branchLength;
            
            cumBranchLength = cumBranchLength+ p->length;
            
            
            
            if ((cumBranchLength < uniform) || /* => there will be no change */
                
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
                
                ran = Random::randomUniformFromGsl2(randomGsl)  * cumMij[ancstate][3];
                
                for (j=0; j<4; j++)
                    
                {
                    
                    if (ran <= cumMij[ancstate][j])
                        
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
                
                numMU++;
                
            }
            
        }
        
        SimulateISMDNAforSite (p->left, genome, site, doISMhaploid, seed,  totalTreeLength, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran,  randomGsl,  rngBoost);
        
        SimulateISMDNAforSite (p->right, genome, site, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran,  randomGsl,  rngBoost);
        
    }
    
}
/********************************** SimulateISMGenotypeforSite ***********************************/
void SimulateISMGenotypeforSite (TreeNode *p, int numberMutations, int site, int doISMhaploid, long int *seed, long double  scaledTotalTreeLength, std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4],long double  mutationRate, long double  &uniform, long double  &cumBranchLength, long double  &ran, const gsl_rng *randomGsl, boost::mt19937* rngBoost )
{
    
    //    static long double     cumBranchLength, uniform, ran;
    
    int             cell, anccell, maternal_ancstate, paternal_ancstate, ancestral_genotype, newstate, maternal_newstate, paternal_newstate;
    
    
    
    if (p != NULL)
        
    {
        
        cell = p->label;
        
        
        
        if ( p->anc1 == NULL)
            
        {
            
            cumBranchLength = 0;
            
              uniform = Random::randomUniformFromGsl2(randomGsl) * scaledTotalTreeLength;
            
        }
        
        else
            
        {
            
            anccell = p->anc1->label;
            
             maternal_ancstate = p->anc1->maternalSequence[site];
             paternal_ancstate = p->anc1->paternalSequence[site];
            
            maternal_newstate = maternal_ancstate;
            paternal_newstate = paternal_ancstate;
            
            ancestral_genotype = Utils::WhichGenotypeIndex(maternal_ancstate, paternal_ancstate);
            
            //            cumBranchLength += p->length;// ->branchLength;
            
            cumBranchLength = cumBranchLength+ p->length*mutationRate;
            
            
            
            if ((cumBranchLength < uniform) || /* => there will be no change */
                
                (allSites[site].numMutations >= numberMutations))
                
            {
                
                
                    
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                
                
                    
                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                
                
                
                p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                
                p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                
                p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site];
                
            }
            
            else /* => there will be change */
                
            {
                long double  row_Qij[16] = { 0.0 };
                
                
                FillRowQGTJCSubstitutionMatrix ( row_Qij, ancestral_genotype);
                //FillGenotypeGTJCSubstitutionMatrix(Qij, p->length, mutationRate, freq, maxEntry, maxPerRow);
                
                //ran = Random::randomUniformFromGsl2(randomGsl)  * cumMij[ancestral_genotype][3];
                newstate = Random::ChooseUniformState(row_Qij, seed, true, randomGsl, rngBoost);
                
                Utils::WhichMaternalPaternalIndex (newstate, maternal_newstate, paternal_newstate);
                

                
                if (maternal_ancstate != maternal_newstate){
                    
                    allSites[site].numMutationsMaternal++;
                    
                    p->maternalSequence[site]= maternal_newstate;
                                   
                    p->paternalSequence[site]=p->anc1->paternalSequence[site];
                    
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site]+1;
                    
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site];
                    
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                        numMU++;
                }
                
                else if (paternal_ancstate != paternal_newstate){
                    
                    allSites[site].numMutationsPaternal++;
                    
                    p->paternalSequence[site]= paternal_newstate;
                                                     
                    p->maternalSequence[site]=p->anc1->maternalSequence[site];
                    
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                    
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site]+1;
                    
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                      numMU++;
                    
                }
                
                allSites[site].numMutations++;
              
                if (allSites[site].numMutations ==numberMutations)
                    return;
                  
                
            }
            
        }
        
        SimulateISMGenotypeforSite (p->left, numberMutations, site, doISMhaploid, seed,  scaledTotalTreeLength, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran,  randomGsl,  rngBoost);
        
        SimulateISMGenotypeforSite (p->right, numberMutations, site, doISMhaploid, seed, scaledTotalTreeLength, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran,  randomGsl,  rngBoost);
        
    }
    
}

/********************************** SimulateISMForSite ***********************************/
/*    Simulates a 0/1 mutation under an infinite sites model (ISM) for a given site. The branch
 where this mutation is placed is chosen according to its length.
 0 is the reference (healthy) allele
 */
void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, long double  totalTreeLength, std::vector<SiteStr> &allSites, int  &numMU, long double  cumMij[4][4], long double  mutationRate, long double  &cumBranchLength, long double  &uniform, int &mutationAdded, const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
 
    int             cell, anccell;
    //    if (*mutationAdded==YES)
    //        return;
    if (p != NULL)
    {
        cell = p->label;
        
        //        if (p->isOutgroup  == YES || p->anc1 == NULL)
        if ( p->anc1 == NULL)
        {
            cumBranchLength = 0;
            //            long double  rUniform=RandomUniform(seed) * totalTreeLength;
            uniform = Random::RandomUniform(seed) * totalTreeLength;
            
        }
        else
        {
            anccell = p->anc1->label;
           
            cumBranchLength = cumBranchLength+ p->length;
     
            if ((cumBranchLength < uniform) ||// (*mutationAdded==YES)/* => there will be no change */
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
                    numMU=numMU+1;
                    mutationAdded=YES;
                }
                else if(genome == PATERNAL && p->anc1->paternalSequence[site]==0 )
                {
                    p->paternalSequence[site]=1;
                    p->numbersPaternalMutationsPerSite[site]=p->anc1->numbersPaternalMutationsPerSite[site]+1;
                    p->numbersMaternalMutationsPerSite[site]=p->anc1->numbersMaternalMutationsPerSite[site];
                    p->numbersMutationsUnderSubtreePerSite[site]=p->anc1->numbersMutationsUnderSubtreePerSite[site]+1;
                    numMU=numMU+1;
                    mutationAdded=YES;
                }
            
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
      
        SimulateISMforSite (p->left, genome, site, doISMhaploid, seed, totalTreeLength, allSites, numMU, cumMij, mutationRate, cumBranchLength, uniform, mutationAdded, randomGsl, rngBoost);
        SimulateISMforSite (p->right, genome, site, doISMhaploid, seed,totalTreeLength, allSites, numMU,cumMij, mutationRate, cumBranchLength,  uniform, mutationAdded, randomGsl,  rngBoost);
        
    
    }
}
/************************************* SimulateFiniteDNA **********************************************/
/* Simulates the nucleotide substitution process under a 4-state Markov model including JC, HKY, GTR and GTRnr */
/* Note that beta is set such that mean substitution rate will be 1.
 E.g., for J-C model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, long double &freqR, long double   &freqY,  long double  &freqAG,  long double  &freqCT,  double  titv,  double  freq[4],  double  Mij[4][4], int numAltModelSites, std::vector<int> &AltModelSites, std::vector<SiteStr> &allSites,  int rateVarAmongSites, long double  altModelMutationRate, int &numMU,  double  Root[],  double  Cijk[], const gsl_rng *rngGsl, boost::mt19937* rngBoost)
{
    int     i, j;
 
    long double  beta, kappa;
    double  Qij[16];
     double  mr;
    
    if (doJC == YES)
    {
        beta = 4./3;
    }
    else if (doHKY == YES)
    {
       
        freqR = freq[Definitions::A] + freq[Definitions::G];
        freqY = freq[Definitions::C] + freq[Definitions::T];
        freqAG = freq[Definitions::A] * freq[Definitions::G];
        freqCT = freq[Definitions::C] * freq[Definitions::T];
        kappa = (titv* freqR*freqY)/(freqAG+freqCT);
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
        linalgebra::EigenREV(mr, Qij, Root, Cijk);
    }
  
    for (i=0; i<numAltModelSites; i++)
        SimulateFiniteDNAforSite (p,  genome, AltModelSites[i], allSites,  seed,  rateVarAmongSites,  altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,    beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk, rngGsl, rngBoost);
    
}
/************************************* SimulateFiniteDNAGenotype **********************************************/


void SimulateFiniteDNAGenotype (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, long double &freqR, long double   &freqY,  long double  &freqAG,  long double  &freqCT,  double  titv,  double  freq[4],  double  Mij[4][4], int numAltModelSites, std::vector<int> &AltModelSites, std::vector<SiteStr> &allSites,  int rateVarAmongSites, long double  altModelMutationRate, int &numMU,  double  Root[],  double  Cijk[], const gsl_rng *rngGsl, boost::mt19937* rngBoost)
{
    int     i, j;
 
    long double  beta = 0, kappa = 0;
    double  Qij[16];
    double  mr;
    
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
        linalgebra::EigenREV(mr, Qij, Root, Cijk);
   
  
    for (i=0; i<numAltModelSites; i++)
        SimulateFiniteDNAforSite (p,  genome, AltModelSites[i], allSites,  seed,  rateVarAmongSites,  altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,    beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk, rngGsl, rngBoost);
    
}
/*********************************** JC **************************************/
/*    JC performs Jukes-Cantor 69 correction */
/* Note that beta was set such that mean substitution rate will be 1.
 for the JC model, beta=4/3, where 12(1/4)(1/4)(4/3) = 1.      */
void JCmodel (long double  Pij[4][4], long double  branchLength, long double  beta )
{
    int i, j;
    
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if (i == j)
                Pij[i][j] = 0.25 + 0.75*exp(beta*(-branchLength));
            else
                Pij[i][j] = 0.25 - 0.25*exp(beta*(-branchLength));
        }
    }
}
/*********************************** GTJC16 **************************************/
void GTJCmodel (long double  Pij[16][16], long double  branchLength, long double  theta, double maxEntry, double  maxPerRow[16]  )
{
    int i, j;
    
    maxEntry = 0.0;
    double exp_4_3 = exp(-4.0*branchLength*theta/3.0);
    double exp_2_3 = exp(-2.0*branchLength*theta/3.0);
    double b = (1.0/16.0)*exp_4_3*(exp_2_3 - 1)*(exp_2_3 -1);
    
    for (i=0; i<16; i++)
    {
        for (j=0; j<16; j++)
        {
  
                Pij[i][j] = b;
        }
    }
    double     c = (1.0/16.0)*(1-3*exp_4_3+2*exp_2_3);
    double     a = (1.0/16.0)*exp_4_3*(exp_2_3 +3)*(exp_2_3 +3);
    maxEntry = std::max(b,c);
    
    for(j=0; j<16; j++){
          
             i=j+1;
             if (i== 5 || i== 6 || i== 7 || i== 11 || i== 12 || i== 13){
               Pij[j][0]= c;
              }
          
             if (i== 5 || i== 8 || i== 9 || i== 11 || i== 14 || i== 15){
                Pij[j][1]= c;

              }
          
             if (i== 6 || i== 8 || i== 10 || i== 12 || i== 14 || i== 16){
                Pij[j][2]= c;
              }
          
             if (i== 7 || i== 9 || i== 10 || i== 13 || i== 15 || i== 16){
                Pij[j][3]= c;
              }
        
             if (i== 1 || i== 2 || i== 6 || i== 7 || i== 14 || i== 15){
                Pij[j][4]= c;
              }
           
             if (i== 1 || i== 3 || i== 5 || i== 7 || i== 8 || i== 16){
                Pij[j][5]= c;
              
              }
            
             if (i== 1 || i== 4 || i== 5 || i== 6 || i== 9 || i== 10){
                Pij[j][6]= c;
               
              }
       
             if (i== 2 || i== 3 || i== 6 || i== 9 || i== 11 || i== 16){
                Pij[j][7]= c;
              }
          
             if (i== 2 || i== 4 || i== 7 || i== 8 || i== 10 || i== 11){
                Pij[j][8]= c;
     
              }
          
             if (i== 3 || i== 4 || i== 7 || i== 9 || i== 12 || i== 14){
                Pij[j][9]= c;
              }
            
             if (i== 1 || i== 2 || i== 8 || i== 9 || i== 12 || i== 13){
                Pij[j][10]= c;
               
              }
         
             if (i== 1 || i== 3 || i== 10 || i== 11 || i== 13 || i== 14){
               Pij[j][11]= c;
           
              }
           
             if (i== 1 || i== 4 || i== 11 || i== 12 || i== 15 || i== 16){
              Pij[j][12]= c;
              
              }
          
             if (i== 2 || i== 3 || i== 5 || i== 10 || i== 12 || i== 15){
                Pij[j][13]= c;
              
              }
          
             if (i== 2 || i== 4 || i== 5 || i== 13 || i== 14 || i== 16){
                Pij[j][14]= c;
               
              }
        
             if (i== 3 || i== 4 || i== 6 || i== 8 || i== 13 || i== 15){
                Pij[j][15]= c;
              }
           Pij[j][j]= a;
         }
    
}
/*********************************** HKY **************************************/
/*    HKY performs Hasegawa-Kishino-Yano 85 correction */

void HKYmodel (long double  Pij[4][4], long double  branchLength, long double  kappa, long double  freqR, long double  freqY, long double  beta,  double  freq[4])
{
    int            i, j;
    long double         AA1, t, PIj;
    
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
void GTRmodel (long double  Pij[4][4], long double  branchLength,  double  Root[],  double  Cijk[])
{
    int     i, j, k;
    long double     t, expt[4];
    
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
void FillSubstitutionMatrix (long double  ch_prob[4][4], long double  branchLength, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta, long double  kappa, long double  freqR, long double  freqY,  double  freq[4],  double  Root[],  double  Cijk[])
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
/********************* FillGenotypeGTJCSubstitutionMatrix **********************/

void FillGenotypeGTJCSubstitutionMatrix (long double  ch_prob[16][16], long double  branchLength, long double  theta,  double  freq[16], double maxEntry, double  maxPerRow[16] )
{
    int i, j;
 
    if (branchLength<1e-6)
    {
        for (i=0; i<16; i++)
        {
            for (j=0; j<16; j++)
            {
                if (i == j)
                    ch_prob[i][j] = 1.0;
                else
                    ch_prob[i][j] = 0.0;
            }
            maxPerRow[i]  = 0.0;
        }
        maxEntry = 0.0;
    }
    else
        GTJCmodel (ch_prob,  branchLength,  theta,  maxEntry,   maxPerRow );
    
  
}
/********************* FillGenotypeSubstitutionMatrix **********************/

void FillRowQGTJCSubstitutionMatrix (long double  row[16], int indexCurrentGenotype)
{
    
    long double one_over_6 = 1.0 / 6.0;
  
    switch(indexCurrentGenotype){
        case 0:
          row[4] = one_over_6;
          row[5] = one_over_6;
          row[6] = one_over_6;
          row[10] = one_over_6;
          row[11] = one_over_6;
          row[12] = one_over_6;
          break;
        case 1:
          row[4] = one_over_6;
          row[7] = one_over_6;
          row[8] = one_over_6;
          row[10] = one_over_6;
          row[13] = one_over_6;
          row[14] = one_over_6;
          break;
        case 2:
          row[5] = one_over_6;
          row[7] = one_over_6;
          row[9] = one_over_6;
          row[11] = one_over_6;
          row[13] = one_over_6;
          row[15] = one_over_6;
          break;
        case 3:
          row[6] = one_over_6;
          row[8] = one_over_6;
          row[9] = one_over_6;
          row[12] = one_over_6;
          row[14] = one_over_6;
          row[15] = one_over_6;
          break;
        case 4:
          row[0] = one_over_6;
          row[1] = one_over_6;
          row[5] = one_over_6;
          row[6] = one_over_6;
          row[13] = one_over_6;
          row[14] = one_over_6;
          break;
        case 5:
          row[0] = one_over_6;
          row[2] = one_over_6;
          row[4] = one_over_6;
          row[6] = one_over_6;
          row[7] = one_over_6;
          row[15] = one_over_6;
          break;
        case 6:
          row[0] = one_over_6;
          row[3] = one_over_6;
          row[4] = one_over_6;
          row[5] = one_over_6;
          row[8] = one_over_6;
          row[9] = one_over_6;
          break;
        case 7:
          row[1] = one_over_6;
          row[2] = one_over_6;
          row[5] = one_over_6;
          row[8] = one_over_6;
          row[10] = one_over_6;
          row[15] = one_over_6;
          break;
        case 8:
          row[1] = one_over_6;
          row[3] = one_over_6;
          row[6] = one_over_6;
          row[7] = one_over_6;
          row[9] = one_over_6;
          row[10] = one_over_6;
          break;
        case 9:
          row[2] = one_over_6;
          row[3] = one_over_6;
          row[6] = one_over_6;
          row[8] = one_over_6;
          row[11] = one_over_6;
          row[13] = one_over_6;
          break;
        case 10:
          row[0] = one_over_6;
          row[1] = one_over_6;
          row[7] = one_over_6;
          row[8] = one_over_6;
          row[11] = one_over_6;
          row[12] = one_over_6;
          break;
        case 11:
          row[0] = one_over_6;
          row[2] = one_over_6;
          row[9] = one_over_6;
          row[10] = one_over_6;
          row[12] = one_over_6;
          row[13] = one_over_6;
          break;
        case 12:
          row[0] = one_over_6;
          row[3] = one_over_6;
          row[10] = one_over_6;
          row[11] = one_over_6;
          row[14] = one_over_6;
          row[15] = one_over_6;
          break;
        case 13:
          row[1] = one_over_6;
          row[2] = one_over_6;
          row[4] = one_over_6;
          row[9] = one_over_6;
          row[11] = one_over_6;
          row[14] = one_over_6;
          break;
        case 14:
          row[1] = one_over_6;
          row[3] = one_over_6;
          row[4] = one_over_6;
          row[12] = one_over_6;
          row[13] = one_over_6;
          row[15] = one_over_6;
          break;
        case 15:
          row[2] = one_over_6;
          row[3] = one_over_6;
          row[5] = one_over_6;
          row[7] = one_over_6;
          row[12] = one_over_6;
          row[14] = one_over_6;
          break;
            
    }

}
/************************************* SimulateFiniteDNAforSite **********************************************/
/* Simulates JC, HKY, GTR or GTRnr for a given site */
void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,std::vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, long double  altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, long double  beta,  long double  kappa,  long double&  freqR,  long double&  freqY,  double  freq[4],  double  Root[],  double  Cijk[], const gsl_rng *randomGsl, boost::mt19937* rngBoost)
{
    long double     branchLength;
    long double Pij[4][4];
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
                    branchLength =  altModelMutationRate* p->length * allSites[site].rateMultiplier;
                }
                else{
                    branchLength =  altModelMutationRate * p->length ;
                    // branchLength = altModelMutationRate * p->length;
                    
                }
                
                
                FillSubstitutionMatrix (Pij, branchLength,
                                        doJC, doHKY,  doGTR,  doGTnR, beta,  kappa,  freqR,  freqY,  freq,  Root,  Cijk);
                
                
                if (genome ==MATERNAL )
                    newstate =p->maternalSequence[site]=Random::ChooseUniformState (Pij[ancstate], seed, true, randomGsl,  rngBoost);
                else// paternal;
                    newstate =p->paternalSequence[site]=Random::ChooseUniformState (Pij[ancstate], seed, true, randomGsl, rngBoost);
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
                    allSites[site].isSNV = YES;
                    numMU++;
                }
            }
        }
        SimulateFiniteDNAforSite (p->left,  genome, site,allSites, seed,  rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY, doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk, randomGsl, rngBoost);
        SimulateFiniteDNAforSite (p->right, genome, site,allSites, seed, rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq ,Root,  Cijk, randomGsl, rngBoost);
    }
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



void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures,  double  cumfreq[4],  long double  *triNucFreq, const gsl_rng *randomGsl, boost::mt19937* rngBoost)

{
    
    int         chosenTriNucleotide, nextNucleotide;
    
    int            k, n1, n2, n3, rest, site;
    
    long double          *prob4, sum;
    
    
    
    /* memory allocations */
    
    prob4 = (long double  *) calloc (4, sizeof(long double ));
    
    if (!prob4)
        
    {
        
        fprintf (stderr, "Could not allocate the prob4 vector\n");
        
        exit (-1);
        
    }

    /* choose first trinucleotide */
    
    chosenTriNucleotide = Random::ChooseUniformState(triNucFreq, seed, true, randomGsl, rngBoost);
    
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
        
        
        
        nextNucleotide = Random::ChooseUniformState(prob4, seed, true, randomGsl, rngBoost);
        
        p->maternalSequence[site]=p->paternalSequence[site]= allSites[site].referenceAllele = nextNucleotide;
        
        
        
        /* move downstream one position */
        
        n1 = n2;
        
        n2 = n3;
        
        n3 = nextNucleotide;
        
    }
    
    
    
    free(prob4);
    
    prob4=NULL;
    
}

/********************************** EvolveSitesOnTree ***********************************/
/* Evolves all sites (maternal and paternal genomes) on the given tree
 We assume that a site will be ISM or Mk in both maternal and paternal genome
 */
void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, std::vector<int> &DefaultModelSites, std::vector<int> &AltModelSites,  long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, int  &numMU, long double  cumMij[4][4], int altModel, long double  altModelMutationRate, int doUserTree, int doJC, int doHKY, int doGTR, int doGTnR, long double  freqR, long double  freqY, long double  freqAG, long double  freqCT, double  titv, double  freq[4], double  Mij[4][4],   double  Root[],  double  Cijk[],
                     const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost  )
{
    int i = 0;
    
    if (rateVarAmongSites == YES)
        for (i=0; i<numSites; i++)
            allSites[i].rateMultiplier = Random::RandomGamma (alphaSites, seed, true, rngGsl,   rngBoost) / alphaSites;
    
    if (propAltModelSites == 0)/* only default model (ISM diploid) sites */
    {
        numDefaultModelSites = numSites;
        numAltModelSites = 0;
        for (i=0; i<numSites; i++)
            DefaultModelSites[i] = i;
        SimulateISM (treeRoot, genome, NO, seed,
                     DefaultModelSites, numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij,altModelMutationRate, rngGsl,   rngBoost);
        
        
    }
    else if (propAltModelSites == 1)
    {
        numDefaultModelSites = 0;
        numAltModelSites = numSites;
        for (i=0; i<numSites; i++)
            AltModelSites[i] = i;
        
        if (altModel == ISMhap)
        { std::cout<<  "only non ISM sites" << std::endl;
            SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij,    altModelMutationRate,  rngGsl,   rngBoost);
        }
        //        else if (altModel == Mk2)
        //        {
        //            SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, numMU);
        //        }
        else if (altModel == finiteDNA)
        {
            SimulateFiniteDNA (treeRoot, genome, seed,doJC,  doHKY, doGTR,  doGTnR,  freqR, freqY,  freqAG,  freqCT, titv,  freq, Mij, numAltModelSites, AltModelSites, allSites,   rateVarAmongSites, altModelMutationRate,  numMU,   Root,  Cijk, rngGsl,  rngBoost );
            
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
            if (Random::RandomUniform (seed) < propAltModelSites)
                AltModelSites[numAltModelSites++] = i;
            else
                DefaultModelSites[numDefaultModelSites++] = i;
        }
        //            /* Evolve ISM sites */
        if (numDefaultModelSites > 0)
        {
            SimulateISM (treeRoot, genome, NO, seed,
                         DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij,altModelMutationRate, rngGsl,   rngBoost
                         );
            
        }
        //            /* Evolve non-SIM sites */
        if (numAltModelSites > 0)
        {
            if (altModel == ISMhap)
            {
                SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,  allSites, numMU, cumMij,    altModelMutationRate, rngGsl,   rngBoost);
            }
            //            else if (altModel == Mk2)
            //            { SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, numMU);
            //            }
            else if (altModel == finiteDNA)
            {
                SimulateFiniteDNA (treeRoot, genome, seed,doJC,  doHKY, doGTR,  doGTnR,  freqR, freqY,  freqAG,  freqCT, titv,  freq, Mij, numAltModelSites,AltModelSites,  allSites,    rateVarAmongSites,  altModelMutationRate, numMU,   Root,  Cijk, rngGsl,  rngBoost );
            }
            else
            {
                fprintf (stderr, "\n\nERROR: Sorry, the specified model is unknown for me");
                exit(-1);
            }
        }
    }
    
}
/********************************** EvolveGenotypesOnTree ***********************************/

void EvolveGenotypesOnTree (TreeNode *treeRoot,  std::vector<int> numberOfSitesWithKMutations, int numberVariableSites, long int *seed, int rateVarAmongSites, int numSites, std::vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, std::vector<int> &DefaultModelSites, std::vector<int> &AltModelSites,  long double  totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, int  &numMU, long double  cumMij[4][4], int altModel, long double  mutationRate, int doUserTree, int doJC, int doHKY, int doGTR, int doGTnR, long double  freqR, long double  freqY, long double  freqAG, long double  freqCT, double  titv, double  freq[4], double  Mij[4][4],   double  Root[],  double  Cijk[],
                     const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost  )
{
    
        
      long double  scaledTotalTreeLength = totalTreeLength * mutationRate;
    
        SimulateGenotype (treeRoot,  numberOfSitesWithKMutations,  numberVariableSites, NO, seed,
                                        DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  scaledTotalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij, mutationRate, rngGsl,   rngBoost
                                        );
                
   
    
}
void EvolveCNLOHonTree (TreeNode *p, int genome,
                        int &numISMCNLOH,
                        std::vector<SiteStr> &allSites,
                        int numSites,
                        long int *seed,
                        long double CNLOHrate,
                        long double mutationRate,
                        long double totalTreeLength,
                        long double &cumCNLOHbranchLength,
                        const gsl_rng *rngGsl,
                        boost::random::mt19937 * rngBoost )
    {
    int        i, trials, numCNLOH, CNLOHsoFar;
    double    totalCNLOHbranchSum;

    totalCNLOHbranchSum = numSites * totalTreeLength / mutationRate * CNLOHrate;
    numCNLOH = Random::RandomPoisson (totalCNLOHbranchSum, seed,rngGsl,  rngBoost ); /* the number of CNLOH events will be distributed as a Poisson with parameter totalCNLOHbranchSum */
    
     /* if the number of CNLOH events is bigger than the number of sites quit the program and warn the user about violation of ISM */
     if (numISMCNLOH + numCNLOH > (2*numSites))
        {
        fprintf (stderr, "\n\nERROR: The haploid infinite sites model (ISM) for CNLOH events has been violated. There will be");
        fprintf (stderr, "\nmore CNLOH events (%d existing + %d proposed]) than available sites under this model (%d).", numISMCNLOH, numCNLOH, 2*numSites);
        fprintf (stderr, "\nTry using a smaller deletion rate.");
        fprintf (stderr, "\nIf using a user tree, try introducing smaller branch lengths.\n");
        exit (-1);
        }
    
    numISMCNLOH += numCNLOH;
    trials = 0;
    CNLOHsoFar = 0;
    while (CNLOHsoFar < numCNLOH)
        {
        i = Random::RandomUniformTo(numSites, seed, true,rngGsl,rngBoost);  /* choose a site at random */
        while ((genome == MATERNAL && allSites[i].numCNLOHmaternal != 0) ||
               (genome == PATERNAL && allSites[i].numCNLOHpaternal != 0))
            {
            i = Random::RandomUniformTo(numSites, seed, true,rngGsl,rngBoost  );
            if (trials++ > 100*numSites)
                {
                fprintf (stderr, "\n\n ERROR: after %d trials cannot find a site without CNLOH",100*numSites);
                fprintf (stderr, "\nCNLOH = %d  CNLOH so far = %d\n\n",numCNLOH, CNLOHsoFar);
                for (i=0; i<numSites; i++)
                    fprintf (stderr, "\nsite %d has %d CNLOHs",i+1, allSites[i].numCNLOH);
                 fprintf (stderr, "\n");
                 exit(-1);
                }
            }
    
        SimulateCNLOHforSite (p, genome, i, seed, allSites, CNLOHrate, cumCNLOHbranchLength,  totalTreeLength, mutationRate, numCNLOH, rngGsl,   rngBoost);
            
        CNLOHsoFar++;
            
        #ifdef MYDEBUG
            fprintf (stderr, "\nCNLOH = %d   nCNLOH so far = %d\n",numCNLOH, CNLOHsoFar);
            if (allSites[i]].numDeletions > 1)
                {
                fprintf (stderr, "\n\n ERROR: %d CNLOH in site %d",allSites[i].numCNLOH, i+1);
                exit(-1);
                }
            for (i=0; i<numSites; i++)
                fprintf (stderr, "%2d[%d] ",i+1, allSites[i].numCNLOH);
        #endif
            
        }
    }


/********************************** SimulateCNLOHforSite ***********************************/
/*    Simulates a CNLOH event infinite sites model (ISM) for a given site. The branch
    where this mutation is placed is chosen according to its length.
  */
void SimulateCNLOHforSite (TreeNode *p, int genome, int site,long int *seed,std::vector<SiteStr> &allSites,long double CNLOHrate, long double  &cumCNLOHbranchLength, long double totalTreeLength, long double mutationRate, int  &numCNLOH, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
    {
   // static double    cumCNLOHbranchLength, uniform;
        double uniform=0.0;
    int             cell, anccell;
        
    if (p != NULL)
        {
        cell = p->label;
            
        if (p->anc1 == NULL)
            {
            cumCNLOHbranchLength = 0;
            uniform = Random::RandomUniform(seed) * totalTreeLength / mutationRate * CNLOHrate;;
            }
        else
            {
            anccell = p->anc1->label;
            cumCNLOHbranchLength += p->length / mutationRate * CNLOHrate;
                
             if ((cumCNLOHbranchLength < uniform) || /* => there will be no change at this branch */
                ((genome == MATERNAL) && (allSites[site].numCNLOHmaternal > 0)) ||
                ((genome == PATERNAL) && (allSites[site].numCNLOHpaternal > 0)))
                {
                
                if (genome == MATERNAL)
                    p->maternalSequence[site] = p->anc1->maternalSequence[site];
                else
                    p->paternalSequence[site] = p->anc1->paternalSequence[site];
                }
            else /* => there will be a CN_LOH event */
                {
                if (genome == MATERNAL)
                    {
                    p->maternalSequence[site] = p->paternalSequence[site];
                    allSites[site].numCNLOHmaternal++;
                    }
                else
                    {
                    p->paternalSequence[site] = p->maternalSequence[site];
                    allSites[site].numCNLOHpaternal++;
                    }

                allSites[site].numCNLOH++;
                numCNLOH++;
                }
            }
            SimulateCNLOHforSite (p->left, genome, site, seed,  allSites, CNLOHrate, cumCNLOHbranchLength,  totalTreeLength, mutationRate, numCNLOH, rngGsl,   rngBoost);
            SimulateCNLOHforSite (p->right, genome, site, seed,  allSites, CNLOHrate, cumCNLOHbranchLength,  totalTreeLength, mutationRate, numCNLOH, rngGsl,   rngBoost);
        }
    }
/********************************** EvolveDeletionsOnTree ***********************************/
 /*    Simulates point deletions under an infinite haploid sites model (ISM).
    Haploid here means that site 33 in the maternal genome and site 33 in the paternal
    genome will be considered different sites, so both can mutate.
 
   This is a very simple mode, where deletions are simulated after the sequences have been evolved
   (i.e., SNVs are already in place)
 
    The number of point deletions will be distributed as a Poisson with parameter the sum of
    branch lengths (node length * deletion rate) over all sites.
    We will force it to get at most one deletion per site
*/
void EvolveDeletionsOnTree (TreeNode *p, int genome, std::vector<SiteStr> &allSites, int &numISMdeletions,int numSites,  long double totalTreeLength, long double mutationRate, long double deletionRate,int &numDEL,  long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
    {
    int        i, trials, numDeletions, deletionsSoFar;
    double    totalDeletionBranchSum;
        long double cumDeletionBranchLength;
    /*
    totalDeletionBranchSum = 0;
    for (i=0; i<numSites; i++)
        {
        allSites[i].deletionBranchSum = SumBranches (healthyRoot) / mutationRate * deletionRate;;
        totalDeletionBranchSum += allSites[i].deletionBranchSum;
        }
    */
    
    totalDeletionBranchSum = numSites * totalTreeLength / mutationRate * deletionRate;
    numDeletions = Random::RandomPoisson (totalDeletionBranchSum, seed, rngGsl,  rngBoost ); /* the number of deletions will be distributed as a Poisson with parameter totalDeletionBranchSum */
    
     /* if the number of deletions is bigger than the number of sites quit the program and warn the user about violation of ISM */
     if (numISMdeletions + numDeletions > (2*numSites))
        {
        fprintf (stderr, "\n\nERROR: The haploid infinite sites model (ISM) for deletions has been violated. There will be");
        fprintf (stderr, "\nmore deletions (%d existing + %d proposed]) than available sites under this model (%d).", numISMdeletions, numDeletions, 2*numSites);
        fprintf (stderr, "\nTry using a smaller deletion rate.");
        fprintf (stderr, "\nIf using a user tree, try introducing smaller branch lengths.\n");
        exit (-1);
        }
    
    numISMdeletions += numDeletions;
    trials = 0;
    deletionsSoFar = 0;
    while (deletionsSoFar < numDeletions)
        {
        i = Random::RandomUniformTo(numSites, seed, true,rngGsl,rngBoost);  /* choose a site at random */
        while ((genome == MATERNAL && allSites[i].numDeletionsMaternal != 0) ||
               (genome == PATERNAL && allSites[i].numDeletionsPaternal != 0))
            {
            i = Random::RandomUniformTo(numSites, seed, true,rngGsl,rngBoost);
            if (trials++ > 100*numSites)
                {
                fprintf (stderr, "\n\n ERROR: after %d trials cannot find an undeleted site",100*numSites);
                fprintf (stderr, "\ndeletions = %d   deletions so far = %d\n\n",numDeletions, deletionsSoFar);
                for (i=0; i<numSites; i++)
                    fprintf (stderr, "\nsite %d has %d deletions",i+1, allSites[i].numDeletions);
                 fprintf (stderr, "\n");
                 exit(-1);
                }
            }
    
        SimulateDeletionforSite (p, genome, i, allSites,  totalTreeLength, mutationRate, cumDeletionBranchLength,   deletionRate, numDEL, seed);
            
        deletionsSoFar++;
            
        #ifdef MYDEBUG
            fprintf (stderr, "\deletions = %d   deletions so far = %d\n",numDeletions, deletionsSoFar);
            if (allSites[i]].numDeletions > 1)
                {
                fprintf (stderr, "\n\n ERROR: %d deletions in site %d",allSites[i].numDeletions, i+1);
                exit(-1);
                }
            for (i=0; i<numSites; i++)
                fprintf (stderr, "%2d[%d] ",i+1, allSites[i].numDeletions);
        #endif
            
        }
    }

/********************************** SimulateDeletionforSite ***********************************/
/*    Simulates a point deletion infinite sites model (ISM) for a given site. The branch
    where this mutation is placed is chosen according to its length.
  */
void SimulateDeletionforSite (TreeNode *p, int genome, int site, std::vector<SiteStr> &allSites, long double totalTreeLength, long double mutationRate,long double &cumDeletionBranchLength,  long double deletionRate, int &numDEL, long int *seed)
    {
    long double  uniform= 0.0;
    int             cell, anccell;
        
    if (p != NULL)
        {
        cell = p->label;
            
        if (p->anc1 == NULL)
            {
            cumDeletionBranchLength = 0;
//            uniform = RandomUniform(seed) * allSites[site].deletionBranchSum;
             uniform = Random::RandomUniform(seed) * totalTreeLength / mutationRate * deletionRate;
           }
        else
            {
            anccell = p->anc1->label;
            cumDeletionBranchLength += p->length / mutationRate * deletionRate;
                
             if ((cumDeletionBranchLength < uniform) || /* => there will be no change at this branch */
                ((genome == MATERNAL) && (allSites[site].numDeletionsMaternal > 0)) ||
                ((genome == PATERNAL) && (allSites[site].numDeletionsPaternal > 0)))
                {
                    
                if (genome == MATERNAL)
                    p->maternalSequence[site] = p->anc1->maternalSequence[site];
                else
                    p->paternalSequence[site] = p->anc1->paternalSequence[site];
                }
            else /* => there will be a deletion */
                {
                if (p->anc1->maternalSequence[site] != DELETION)  /* checking all this might be excessive */
                     p->anc1->maternalSequence[site] = DELETION;
                    
                else if (p->anc1->maternalSequence[site] == DELETION)
                    {
                    fprintf (stderr, "\n\n ERROR: site %d in genome %d of cell %d cannot be deleted twice under the ISM model\n", site, genome, cell);
                    exit(-1);
                    }
                else
                    {
                    fprintf (stderr, "\n\n ERROR: site %d in genome %d of cell %d has an unknow state %d under the ISM model\n", site, genome, anccell, p->anc1->maternalSequence[site] );
                    exit(-1);
                    }
                    
                if (genome == MATERNAL)
                    allSites[site].numDeletionsMaternal++;
                else if (genome == PATERNAL)
                    allSites[site].numDeletionsPaternal++;
                allSites[site].numDeletions++;
                numDEL++;
                }
            }
            SimulateDeletionforSite (p->left, genome, site,   allSites,  totalTreeLength, mutationRate, cumDeletionBranchLength, deletionRate,numDEL,seed);
            SimulateDeletionforSite (p->right, genome, site,allSites,  totalTreeLength, mutationRate, cumDeletionBranchLength, deletionRate,numDEL, seed);
        }
    }

/********************* AllelicDropout  ************************/
/*
 Remove alleles from single chromosomes at a given ADO rate per genotype
 We assume that ADO is the product of the allele dropout rate as follows:
 A = genotype ADO; a = allele ADO
 A = 2a - a^2
 a = 1 - sqrt(1-A)
*/

void AllelicDropout (int numCells, std::vector<SiteStr> &allSites, int doADOcell, int doADOsite, int numSites,long double fixedADOrate, long double meanADOcell, long double varADOcell, long double meanADOsite,long double varADOsite, std::vector<TreeNode *> &treeTips,    long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
    {
    int i,j;
    double alleleADOrateMean, alleleADOrateCell, alleleADOrateSite;
   double temp;
    //double **alleleADOrate;

    std::vector<std::vector<double> >   alleleADOrate ;

    // fixed ADO rate
    if (doADOcell == NO && doADOsite == NO)
        {
       
        alleleADOrateMean = 1.0 - sqrt (1.0 - fixedADOrate);
        //alleleADOrateMean = fixedADOrate;
            
            for (i=0; i<treeTips.size(); i++){
                std::vector<double> v;
                for (j=0; j<numSites; j++)
                    v.push_back(alleleADOrateMean);
                alleleADOrate.push_back(v);
                
            }
            
        }
    else // variable ADO rates
        {
        for (i=0; i<treeTips.size(); i++)
            {
            std::vector<double> v;
            if (doADOcell == YES)
                alleleADOrateCell = Random::RandomBetaMeanVar(meanADOcell, varADOcell, seed, true, rngGsl,   rngBoost);
            else
                alleleADOrateCell = 0;
            
            for (j=0; j<numSites; j++)
                {
                if (doADOsite == YES)
                    alleleADOrateSite = Random::RandomBetaMeanVar(meanADOsite, varADOsite, seed,  true, rngGsl,   rngBoost);
                else
                    alleleADOrateSite = 0;
                temp=alleleADOrateCell + alleleADOrateSite - (alleleADOrateCell * alleleADOrateSite);
                v.push_back( 1.0 - sqrt (1.0 -temp)); // at the allele level
                }
               alleleADOrate.push_back(v);
            }
        }
    

    addAllelicDropoutToTree(treeTips, allSites,  numSites, alleleADOrate, seed, rngGsl,  rngBoost );
    

    }

void addAllelicDropoutToTree( std::vector<TreeNode *> &treeTips,std::vector<SiteStr> &allSites, int numSites, std::vector<std::vector<double> > &alleleADOrate, long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost ){
    
    TreeNode *tip;
    long double random;
    for (size_t i=0; i<treeTips.size(); i++)//or numCells+1
    {
        tip = treeTips[i];
        for (size_t  site=0; site<numSites; site++)
               {
                   
               if (tip->maternalSequence[site] != DELETION)
                   {
                       random = Random::randomUniformFromGsl2(rngGsl);
                   if (random < alleleADOrate[i][site])
                       {
                       tip->maternalSequence[site] = ADO;
                       allSites[site].hasADO = YES;
                       }
                   }
               if (tip->paternalSequence[site] != DELETION)
                   {
                   if (Random::randomUniformFromGsl2(rngGsl) < alleleADOrate[i][site])
                       {
                       tip->paternalSequence[site] = ADO;
                       allSites[site].hasADO = YES;
                       }
                   }
               }
    }
}

/********************* GenotypeError  ************************/
/*
 Introduce false positive errors directly in the genotypes
 (as opposed to errors just in the reads) (note Eij=0 when i=j)
 
 We assume that the genotype error is the product of the allele error as follows:
 G = genotype error; a = allele error
 G = 2a - a^2
 a = 1 - sqrt(1-E)
*/

void GenotypeError (std::vector<TreeNode *> &treeTips,std::vector<SiteStr> &allSites, int alphabet, int numSites, int numCells,  double meanGenotypingError,  double varGenotypingError, double genotypingError, double Eij[4][4], long int *seed, const gsl_rng *rngGsl,  boost::random::mt19937 * rngBoost)
    {
    size_t     i,j;
    TreeNode *tip;
    //double    genotypingError;
    double alleleError;
    long double *probs;
    long double  error_prob[4][4];
    FillErrorMatrix (error_prob, Eij);
        
    if (alphabet == DNA)
        {
        for (i=0; i<treeTips.size(); i++)//or numCells+1
            {
             tip = treeTips[i];
            for (j=0; j<numSites; j++)
                {
               // genotypingError = Random::RandomBetaMeanVar(meanGenotypingError, varGenotypingError, seed,  true, rngGsl,   rngBoost);
                alleleError = 1.0 - sqrt (1.0 -  genotypingError);

                if ( tip->maternalSequence[j] != ADO && tip->maternalSequence[j] != DELETION && Random::randomUniformFromGsl2(rngGsl) < alleleError)
                    {
                        probs = error_prob[tip->maternalSequence[j]];
                   tip->maternalSequence[j] = Random::ChooseUniformState (probs, seed, true, rngGsl, rngBoost);
                    allSites[j].hasGenotypeError = YES;
                    }
                if (tip->paternalSequence[j] != ADO && tip->paternalSequence[j] != DELETION && Random::randomUniformFromGsl2(rngGsl) < alleleError)
                    {
                        probs = error_prob[tip->paternalSequence[j]];
                    tip->paternalSequence[j] = Random::ChooseUniformState (probs, seed, true, rngGsl, rngBoost);
                    allSites[j].hasGenotypeError = YES;
                    }
                }
            }
        }
    else  /* for binary data */
        {
        for (i=0; i<treeTips.size(); i++)//or numCells+1
            {
                tip = treeTips[i];
            for (j=0; j<numSites; j++)
                {
                genotypingError = Random::RandomBetaMeanVar(meanGenotypingError, varGenotypingError, seed, true, rngGsl,   rngBoost);
                alleleError = 1.0 - sqrt (1.0 -  genotypingError);
                if (tip->maternalSequence[j] != ADO && tip->maternalSequence[j] != DELETION && Random::RandomUniform(seed) < alleleError)
                    {
                    if (tip->maternalSequence[j] == 0)
                        tip->maternalSequence[j] = 1;
                    else if (tip->maternalSequence[j] == 1)
                        tip->maternalSequence[j] = 0;
                    allSites[j].hasGenotypeError = YES;
                    }
                if (tip->paternalSequence[j] != ADO && tip->paternalSequence[j] != DELETION && Random::RandomUniform(seed) < alleleError)
                    {
                    if (tip->paternalSequence[j] == 0)
                        tip->paternalSequence[j] = 1;
                    else if (tip->paternalSequence[j] == 1)
                        tip->paternalSequence[j] = 0;
                    allSites[j].hasGenotypeError = YES;
                    }
                }
            }
        }
    }

void FillErrorMatrix (long double  error_prob[4][4], double Eij[4][4])
{
    size_t i, j;
    
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
            {
                if (i == j)
                    error_prob[i][j] = Eij[i][j];
                else
                    error_prob[i][j] = Eij[i][j];
            }
    }
}

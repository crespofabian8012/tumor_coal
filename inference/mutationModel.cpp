//
//  mutationModel.cpp
//  simul
//
//  Created by Seong-Hwan Jun on 2019-10-08.
//

#include "mutationModel.h"

#include "eigen.hpp"
#include "random.h"

void SimulateISM (TreeNode *treeRoot, int genome, int doISMhaploid, long int *seed,  vector<int> &DefaultModelSites, int numDefaultModelSites, vector<int> &AltModelSites, int numAltModelSites, double totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet,  vector<SiteStr> &allSites, int  &numMU, double cumMij[4][4], double mutationRate)
{
    int   i, trials, numMutations, mutationsSoFar;
    double  totalBranchSum;
    double cumBranchLength =0;
    double uniform =0;
    int mutationAdded;
    double ran=0;
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
        numMutations = RandomPoisson (totalBranchSum, seed);
    }while( numISMmutations + numMutations > (ploidyFactor * numModelSites));
    
    
    numISMmutations = numISMmutations + numMutations;
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
            SimulateISMDNAforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij, mutationRate, uniform, cumBranchLength,  ran);
        else{
            mutationAdded=NO;
            SimulateISMforSite (treeRoot, genome, i, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij,mutationRate ,  cumBranchLength,  uniform,  mutationAdded);
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

/************************************* SimulateMk2ForSite ***************************************/
/* Simulates the nucleotide substitution process for a given site under Mk2 model (see Lewis 2001)
 with equal rates. 0 is the reference (healthy) allele */

//void SimulateMk2forSite (TreeNode *p, int genome, int site, long int *seed, int doUserTree, int rateVarAmongSites, double altModelMutationRate, vector<SiteStr> &allSites, int ***data, int* numMU )
//{
//    double    probOfChange, uniform, branchLength;
//    int     cell, anccell;
//    
//    if (p != NULL)
//    {
//        if (p->isOutgroup == NO)
//        {
//            cell = p->label;
//            anccell = p->anc1->label;
//            
//            if (doUserTree == YES){
//                branchLength = p->lengthModelUnits;//>branchLength;
//                //                 branchLength = p->length;//>branchLength;
//                
//            }
//            else
//            {
//                if (rateVarAmongSites == YES)
//                    branchLength = altModelMutationRate * p->length * allSites[site].rateMultiplier;
//                else
//                    branchLength = altModelMutationRate * p->length;
//            }
//            
//            probOfChange = 0.5 - 0.5 * exp (-2.0 * branchLength);
//            
//            uniform = RandomUniform(seed);
//            if (uniform >= probOfChange) /* => no change */
//                data[genome][cell][site] = data[genome][anccell][site];
//            else /* => there will be change */
//            {
//                if (data[genome][anccell][site] == 0)
//                    data[genome][cell][site] = 1;
//                else
//                    data[genome][cell][site] = 0;
//                
//                if (genome == MATERNAL)
//                    allSites[site].numMutationsMaternal++;
//                else if (genome == PATERNAL)
//                    allSites[site].numMutationsPaternal++;
//                allSites[site].numMutations++;
//                numMU=numMU+1;
//            }
//        }
//        SimulateMk2forSite (p->left,  genome, site, seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
//        SimulateMk2forSite (p->right, genome, site, seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, data, numMU);
//    }
//}

/************************************* SimulateMk2 **********************************************/
/* Simulates the nucleotide substitution process under the Mk2 model (see Lewis 2001),
 also called Cavender-Farris-Neyman CFN  model or Jukes-Cantor (1969) model for two alleles */

//void SimulateMk2 (TreeNode *p, int genome, long int *seed, vector<int> &AltModelSites, int  numAltModelSites, int doUserTree, int rateVarAmongSites, double altModelMutationRate, vector<SiteStr> &allSites, int &numMU)
//{
//    int     i;
//
//    for (i=0; i<numAltModelSites; i++)
//        SimulateMk2forSite (p, genome, AltModelSites[i], seed,  doUserTree,  rateVarAmongSites,  altModelMutationRate, allSites, numMU);
//}

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

void SimulateISMDNAforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, vector<SiteStr> &allSites, int  &numMU, double cumMij[4][4],double mutationRate, double &uniform, double &cumBranchLength, double &ran )
{
    
    //    static double    cumBranchLength, uniform, ran;
    
    int             j, cell, anccell, ancstate;
    
    
    
    if (p != NULL)
        
    {
        
        cell = p->label;
        
        
        
        if ( p->anc1 == NULL)
            
        {
            
            cumBranchLength = 0;
            
            uniform = RandomUniform(seed) * totalTreeLength;
            
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
                
                ran = RandomUniform(seed) * cumMij[ancstate][3];
                
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
        
        SimulateISMDNAforSite (p->left, genome, site, doISMhaploid, seed,  totalTreeLength, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran);
        
        SimulateISMDNAforSite (p->right, genome, site, doISMhaploid, seed, totalTreeLength, allSites, numMU,cumMij,  mutationRate, uniform, cumBranchLength,  ran);
        
    }
    
}

/********************************** SimulateISMForSite ***********************************/
/*    Simulates a 0/1 mutation under an infinite sites model (ISM) for a given site. The branch
 where this mutation is placed is chosen according to its length.
 0 is the reference (healthy) allele
 */
void SimulateISMforSite (TreeNode *p, int genome, int site, int doISMhaploid, long int *seed, double totalTreeLength, vector<SiteStr> &allSites, int  &numMU, double cumMij[4][4], double mutationRate, double &cumBranchLength, double &uniform, int &mutationAdded)
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
            cumBranchLength = 0;
            //            double rUniform=RandomUniform(seed) * totalTreeLength;
            uniform = RandomUniform(seed) * totalTreeLength;
            
        }
        else
        {
            anccell = p->anc1->label;
            //            *cumBranchLength =*cumBranchLength+ p->length;// ->branchLength;
            // *cumBranchLength =*cumBranchLength+ p->lengthModelUnits;// ->branchLength;
            cumBranchLength = cumBranchLength+ p->length;
            //            if ((*cumBranchLength < *uniform) || /* => there will be no change */
            //                ((doISMhaploid == NO)  && (allSites[site].numMutations > 0))
            //               ||
            //               ((doISMhaploid == YES) && (genome == MATERNAL) && (allSites[site].numMutationsMaternal > 0)) ||
            //              ((doISMhaploid == YES) && (genome == PATERNAL) && (allSites[site].numMutationsPaternal > 0))
            //                )
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
        SimulateISMforSite (p->left, genome, site, doISMhaploid, seed, totalTreeLength, allSites, numMU, cumMij, mutationRate, cumBranchLength, uniform, mutationAdded);
        SimulateISMforSite (p->right, genome, site, doISMhaploid, seed,totalTreeLength, allSites, numMU,cumMij, mutationRate, cumBranchLength,  uniform, mutationAdded);
        
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

void SimulateFiniteDNA (TreeNode *p, int genome, long int *seed, int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4], int numAltModelSites, vector<int> &AltModelSites, vector<SiteStr> &allSites,  int rateVarAmongSites, double altModelMutationRate, int &numMU, double Root[], double Cijk[])
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

/************************************* SimulateFiniteDNAforSite **********************************************/
/* Simulates JC, HKY, GTR or GTRnr for a given site */
void SimulateFiniteDNAforSite (TreeNode *p, int genome, int site,vector<SiteStr> &allSites,  long int *seed, int rateVarAmongSites, double altModelMutationRate, int &numMU, int doJC, int doHKY, int doGTR, int doGTnR, double beta,  double kappa, double freqR, double freqY, double freq[4], double Root[], double Cijk[])
{
    double    branchLength, Pij[4][4];
    int     cell, anccell, ancstate, newstate;
    
    int i,j;
    
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
                
                
                if (genome ==MATERNAL )
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
                    numMU++;
                }
            }
        }
        SimulateFiniteDNAforSite (p->left,  genome, site,allSites, seed,  rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY, doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq,  Root,  Cijk);
        SimulateFiniteDNAforSite (p->right, genome, site,allSites, seed, rateVarAmongSites, altModelMutationRate, numMU,  doJC,  doHKY,  doGTR,  doGTnR,  beta,    kappa,   freqR,   freqY,   freq ,Root,  Cijk);
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



void SimulateTriNucFreqGenome (int cell, long int *seed, TreeNode *p, int alphabet, int doUserGenome, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, double cumfreq[4], double *triNucFreq )

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

/********************************** EvolveSitesOnTree ***********************************/
/* Evolves all sites (maternal and paternal genomes) on the given tree
 We assume that a site will be ISM or Mk in both maternal and paternal genome
 */
void EvolveSitesOnTree (TreeNode *treeRoot, int genome, long int *seed, int rateVarAmongSites, int numSites, vector<SiteStr> &allSites, int doGeneticSignatures, int alphaSites, int  propAltModelSites , int numDefaultModelSites, int numAltModelSites, vector<int> &DefaultModelSites, vector<int> &AltModelSites,  double totalTreeLength , int &numISMmutations, int numFixedMutations, int numSNVmaternal, int doSimulateFixedNumMutations,  int alphabet, int  &numMU, double cumMij[4][4], int altModel, double altModelMutationRate, int doUserTree, int doJC, int doHKY, int doGTR, int doGTnR, double freqR, double freqY, double freqAG, double freqCT, double titv, double freq[4], double Mij[4][4],  double Root[], double Cijk[])
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
                     DefaultModelSites, numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij,altModelMutationRate);
        
        
    }
    else if (propAltModelSites == 1)
    {
        numDefaultModelSites = 0;
        numAltModelSites = numSites;
        for (i=0; i<numSites; i++)
            AltModelSites[i] = i;
        
        if (altModel == ISMhap)
        { fprintf(stderr, "only non ISM sites");
            SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij,    altModelMutationRate);
        }
//        else if (altModel == Mk2)
//        {
//            SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, numMU);
//        }
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
                         DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet, allSites, numMU, cumMij,altModelMutationRate
                         );
            
        }
        //            /* Evolve non-SIM sites */
        if (numAltModelSites > 0)
        {
            if (altModel == ISMhap)
            { fprintf(stderr, "pip");
                SimulateISM (treeRoot, genome, YES, seed,DefaultModelSites,numDefaultModelSites, AltModelSites,  numAltModelSites,  totalTreeLength ,  numISMmutations, numFixedMutations, numSNVmaternal,  doSimulateFixedNumMutations,  alphabet,  allSites, numMU, cumMij,    altModelMutationRate);
            }
//            else if (altModel == Mk2)
//            { SimulateMk2 (treeRoot, genome, seed, AltModelSites,   numAltModelSites, doUserTree, rateVarAmongSites,  altModelMutationRate, allSites, numMU);
//            }
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

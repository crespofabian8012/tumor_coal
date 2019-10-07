//
//  proposal.h
//  TumorCellCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 7/8/19.
//  Copyright © 2019 Fausto Fabian Crespo Fernandez. All rights reserved.
//

#ifndef proposal_h
#define proposal_h
void proposalProportionsVector(double *newProportionsvector, double *currentProportionVector, double tuningParameter, int numClones,  long int *seed);
void proposalScalingMoveProportionsVector(double *proportionsVector, int numClones, double omega, long int *seed );
double AcceptanceRateTheta(double *newProportionsvector, TreeNode  *tree, TreeNode *nodes, Population **populations, int numClones, int totalPopSize, double totalDelta);
//void WilsonBaldingMove(Population ** populations, TreeNode* root, TreeNode** nodes, int numNodes,int numClones, long int *seed, int sampleSize, double currentLogLikelihoodSequences, double mutationRate, char *ObservedCellNames[], pll_msa_t * msa);
//static void WilsonBaldingMove(Chain *chain, long int *seed, int sampleSize, double currentLogLikelihoodSequences, ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa);
void WilsonBaldingMove(Chain *chain, long int *seed, double currentLogLikelihoodSequences, ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa);
void newTotalEffectivePopulationSizeMove(Chain *chain, long int *seed,  ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions * mcmcOptions, int *sampleSizes);
void pllmod_utree_rollback_spr1(pll_tree_rollback_t *rollback_info,pll_utree_rb_t *rb, TreeNode   **nodes,
                                double old_father_pop_effect_pop_size,
                                double new_father_pop_effect_pop_size,
                                TreeNode *container_of_p,
                                TreeNode *container_of_r,
                                TreeNode *container_of_p_prime,//mrca
                                TreeNode *container_of_r_prime,
                                TreeNode *container_of_u,
                                TreeNode *container_of_v
                                );
void newProportionsVectorMove(Chain *chain, long int *seed,  ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions *mcmcOptions, int *sampleSizes);
void newScaledGrowthRateMoveforPopulation(Chain *chain, Population *popI, long int *seed,  ProgramOptions *programOptions, char *ObservedCellNames[], pll_msa_t * msa, MCMCoptions * mcmcOptions, int *sampleSizes);
#endif /* proposal_h */
/* allocate cell structure */
//        AllocateCellStructure( cell, numCells, numSites);

/* allocate memory for site information (equal for maternal and paternal) */
//                    allSites = (SiteStr*) calloc (numSites, sizeof(SiteStr));
//                    if (!allSites)
//                    {
//                        fprintf (stderr, "Could not allocate the allSites structure\n");
//                        exit (-1);
//                    }
//                    for (i=0; i<numSites; i++)
//                    {
//                        allSites[i].alternateAlleles = (int *) calloc (4, sizeof(int));
//                        if (!allSites[i].alternateAlleles)
//                        {
//                            fprintf (stderr, "Could not allocate the allSites[i].alternateAlleles[] structure\n");
//                            exit (-1);
//                        }
//                    }

/* the arrays below keep the index for different types of sites */
//                    SNVsites = (int*) calloc (numSites, sizeof(int));
//                    if (!SNVsites)
//                    {
//                        fprintf (stderr, "Could not allocate the SNVsites structure\n");
//                        exit (-1);
//                    }
//                    SFS = (int*) calloc (numSites, sizeof(int));
//                    if (!SFS)
//                    {
//                        fprintf (stderr, "Could not allocate the SNVsites structure\n");
//                        exit (-1);
//                    }
//                    /* the arrays below keep the index for different types of sites */
//                    variantSites = (int*) calloc (numSites, sizeof(int));
//                    if (!variantSites)
//                    {
//                        fprintf (stderr, "Could not allocate the variantSites structure\n");
//                        exit (-1);
//                    }
//
//                    DefaultModelSites = (int*) calloc (numSites, sizeof(int));
//                    if (!DefaultModelSites)
//                    {
//                        fprintf (stderr, "Could not allocate the DefaultModelSites structure\n");
//                        exit (-1);
//                    }
//
//                    AltModelSites = (int*) calloc (numSites, sizeof(int));
//                    if (!AltModelSites)
//                    {
//                        fprintf (stderr, "Could not allocate the AltModelSites structure\n");
//                        exit (-1);
//                    }
//                    /* allocate space for cellNames */
//                    cellNames = (char **) calloc(numNodes, sizeof(char *));
//                    if (!cellNames)
//                    {
//                        fprintf (stderr, "Could not allocate cellNames (%ld)\n", numNodes * sizeof(char *));
//                        exit (-1);
//                    }
//                    for (i=0; i<numNodes; i++)
//                    {
//                        cellNames[i] = (char *) calloc(MAX_NAME, sizeof(char));
//                        if (!cellNames[i])
//                        {
//                            fprintf (stderr, "Could not allocate cellNames[%d] (%ld)\n", i,  MAX_NAME * sizeof(char));
//                            exit (-1);
//                        }
//                    }
/////////////////////////////////
//                    for (i=0; i<numCells+1; i++)
//                        for (j=0; j<numSites; j++)
//                        {   temp=data[MATERNAL][i][j];
//                            //                                  cell[i].site[j].trueMaternalAllele = temp;
//                            temp=data[PATERNAL][i][j];
//                            // cell[i].site[j].truePaternalAllele = temp;
//                        }
//if (ADOrate > 0)
//    AllelicDropout(&seed);

//InitializeGenomes (healthyRoot, &seed);

/* add germline variation if needed */
// if (SNPrate > 0)
//    AddGermlineVariation (&seed);

/* count how many alleles are at each site and how many SNVs we observe now, after ADO and genotype errors/read generation */
//        if (doNGS == YES)
//            numSNVs = CountAllelesInMLGenotypes();
//        else
//            numSNVs = CountAllelesInObservedGenotypes();
/////////////////////////////////////////////////
//        if (doPrintSNVhaplotypes == YES && numSNVs > 0) /* we only print replicates with variation */
//        {
//            if (doPrintSeparateReplicates == NO)
//                fprintf (fpSNVhaplotypes, "[#%d]\n", dataSetNum+1);
//            PrintSNVHaplotypes(fpSNVhaplotypes, NO);
//        }
//        if (doPrintFullGenotypes == YES )
//        {
//            if (doPrintSeparateReplicates == NO)
//                fprintf (fpFullGenotypes, "[#%d]\n", dataSetNum+1);
//            PrintFullGenotypes(fpFullGenotypes);
//        }
//
//                    /* free memory */
//                    for (i=0; i<ploidy; i++)
//                    {
//                        for (j=0; j<2*numCells+1; j++)
//                            free (data[i][j]);
//                        free (data[i]);
//                    }
//                    free (data);
//                      free (allSites);
//                      free (SNVsites);
//                      free (variantSites);
//                      free (DefaultModelSites);
//                      free (AltModelSites);
////                     free(cellNames);

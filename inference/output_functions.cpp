//
//  output_functions.cpp
//  simul
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//

#include "output_functions.hpp"

#include "utils.hpp"

/***************************** PrintUsage *******************************/
/* Prints a short description of program usage */
void PrintUsage()
{
    fprintf (stderr, "\n\nUsage: %s%s []", PROGRAM_NAME, VERSION_NUMBER);
//    fprintf (stderr, "\n\nUsage: %s%s [-n# -x# -c# (# # # # # #) -u# -o# -ttrees.tre -ktimes.txt -## -y# -? -h]", PROGRAM_NAME, VERSION_NUMBER);
//    fprintf (stderr, "\n-n: number of replicates (e.g. -n1000)");
//    fprintf (stderr, "\n-x: haploid/diploid (1 for haploid, 2 for diploid) (e.g. -x2)");
//    fprintf (stderr, "\n-c: number of clones; for each clone: ID number, sample size, population size, birth rate, death rate, time to origin (e.g. -c5\n\t1 4 70000 0.3 0.1 56\n\t2 5 60000 0.4 0.3 110\n\t3 5 95000 0.4 0.3 117\n\t4 4 15000 0.5 0.4 95\n\t5 3 44000 0.3 0.1 53)");
//    fprintf (stderr, "\n-u: mutation rate (e.g. -u9.1e-6)");
//    fprintf (stderr, "\n-o: branch length to the outgroup (root-outgroup) (e.g. -o0.0325) (default is no outgroup)");
//
//    /* outputs */
//    fprintf (stderr, "\n-t: tree file name (e.g. -ttrees.tre)");
//    fprintf (stderr, "\n-k: times file name (e.g. -ktimes.txt)");
//
//    /* other settings */
//    fprintf (stderr, "\n-#: seed (e.g. -#37864287)");
//    fprintf (stderr, "\n-y: noisy (amount of information printed on the screen) (e.g. -y2)");
//    fprintf (stderr, "\n-? -h: Print help\n");
    
    exit(-1);
}

/**************** PrintTrees ***************/
/*  Print unrooted trees to treefile in Newick format */
void PrintTrees(int replicate, TreeNode *treeRootInit,   FILE  *fpTrees, double mutationRate, int doUseObservedCellNames)
{
    /* there isn´t recombination */
    /*fprintf(fpTrees,"Tree.%05d = ", replicate+1);*/
    //    fprintf(fpTrees, "(");
    WriteTree (treeRootInit, mutationRate, fpTrees, doUseObservedCellNames);
    //    fprintf(fpTrees, ");\n");
    fprintf (fpTrees,");\n");
}

/**************** PrintTrees2 ***************/
/*  Print unrooted trees to treefile in Newick format */
void PrintTrees2(int replicate, TreeNode *treeRootInit,   FILE   *fpTrees2 , double mutationRate,char * ObservedCellNames[],int doUseObservedCellNames)
{
    int indexCurrentCell=0;
    
    /* there isn´t recombination */
    /*fprintf(fpTrees2,"Tree.%05d = ", replicate+1);*/
    //   fprintf(fpTrees2, "(");
    WriteTree2 (treeRootInit, mutationRate, fpTrees2, ObservedCellNames, &indexCurrentCell, doUseObservedCellNames);
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
    //asprintf(&currentNewick, *newickString);
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

/********************* PrintTimes **********************/
/* Prints to timesfile a detailed description of
 the tree: nodes, times, branch lengths */

void PrintTimes(int replicate, FILE   *fpTimes, double mutationRate, vector<TreeNode *> &nodes,  int thereisOutgroup)
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

void PrintTimes2(int replicate, FILE  *fpTimes2, double mutationRate,  vector<TreeNode *> &nodes,  int thereisOutgroup)
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

void ListTimes (int j, double mutationRate, vector<TreeNode *> &nodes, FILE *fpTimes, int thereisOutgroup)
{
    /* It does not list superfluous nodes */
    TreeNode  *p;
    int     i = 0;
    
    do
    {
        p = nodes[i];
        
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

void ListTimes2 (int j,  double mutationRate, vector<TreeNode *> &nodes,  FILE *fpTimes2, int thereisOutgroup)
{
    /* It does not list superfluous nodes */
    TreeNode  *p;
    int     i = 0;
    
    do
    {
        p = nodes[i];
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

/********************* getHealthyTip **********************/
/* getHealthyTip*/
TreeNode *getHealthyTip(TreeNode *treeRootInit)
{
    if (treeRootInit !=NULL && treeRootInit->right!=NULL)
        return treeRootInit->right;
    else
        return NULL;
}

/***************************** PrintTrueFullHaplotypes *******************************/
/* Prints observed/ML haplotypes for all sites (variable + invariable) to a file */
void PrintTrueFullHaplotypes (FILE *fp, vector<TreeNode *> &nodes, TreeNode* treeRoot, int numNodes, int doPrintIUPAChaplotypes, int doPrintAncestors, int numSites, int numCells, int alphabet, int doUserTree , int doNGS,   char **cellNames, CellStr            *cell, int        HEALTHY_ROOT, int TUMOR_ROOT , char *cellnames[], int doUseObservedCellName)
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
                p = nodes[i];
                /* print IUPAC haplotype */
                if (p !=NULL){
                    
                    if (p->left==NULL && p->right ==NULL && p->anc1 !=NULL){
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
                p = nodes[i];
                if (p !=NULL){
                    
                    if (p->left==NULL && p->right ==NULL && p->anc1 !=NULL){
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
                p = nodes[i];
                
                if (p->left==NULL && p->right ==NULL && p->anc1 !=NULL){
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
                p = nodes[i];
                if (p !=NULL){
                    
                    if (p->left==NULL && p->right ==NULL && p->anc1 !=NULL){
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


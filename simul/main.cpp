
#include <iostream>
#include <string>

#include "data_utils.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    // take input
    // number of replicates
    // input file containing
    // 1. time of origins
    // 2. number of clones
    // 3. number of cells for each clone
    // 4. population size of each clone
    // 5. birth rate and death rate for each clone
    // output directory
    
    // default evolution model: JC
    
    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    char *ObservedCellNames[MAX_NAME];
    
    int      *CloneNameBegin, *CloneSampleSizeBegin, *ClonePopSizeBegin;
    double   *CloneBirthRateBegin, *CloneDeathRateBegin, *CloneTimeOriginInput;
    double   *CloneGrowthRateBegin;
    double freq[4];
    double cumfreq[4];
    double Mij[4][4];
    double cumMij[4][4];
    
    double Eij[4][4];
    double cumEij[4][4];
   

    double TotalBirthRate, TotalDeathRate;
    int TotalN, i, j, k;

    FILE    *input_file;
    /* Default settings */
    Initialize( Eij, Mij, freq,  &programOptions );
    //size_t seed = std::stol(argv[0]);
    //string output_path = argv[1];
   char* input_path;
     if (argc <= 2)
          input_path = argv[1];
     else{
         fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
         PrintUsage();
         
     }
    
    // 1. call function to parse the input file

        if ((input_file = freopen(input_path, "r", stdin)) != NULL)
        {
            ReadParametersFromFile(&programOptions, &filePaths, &CloneNameBegin, &CloneSampleSizeBegin, &ClonePopSizeBegin, &CloneBirthRateBegin, &CloneDeathRateBegin, &CloneTimeOriginInput,Mij,freq);
        }
        else
        {
            fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
            PrintUsage();
        }
  
    // 2. create and initialize data structures
    
    
    //allocate memory for the population structs
    
    Population **populations = (Population**)malloc (sizeof(struct Population*)  * programOptions.numClones);
    if (!populations)
    {
        fprintf (stderr, "Could not allocate populations (%lu bytes)\n", (programOptions.numClones)  * (long) sizeof(Population*));
        exit (1);
    }
    InitListClones(populations, programOptions.numClones, programOptions.noisy, CloneNameBegin, CloneSampleSizeBegin, CloneBirthRateBegin,  CloneDeathRateBegin, ClonePopSizeBegin, programOptions.TotalNumSequences);
    InitListClonesTimes(populations, programOptions.numClones,  &programOptions.doEstimateTimesOriginClones,  CloneTimeOriginInput  );
     InitNumberNodes(&TotalBirthRate, &TotalDeathRate, &TotalN, populations, &programOptions);
     ListClonesAccordingTimeToOrigin(populations, programOptions.numClones);

      /* set file dirs and names */
    InitFilesPathsOptions(&filePaths, &programOptions);
    
    
    // 3. call function to simulate the data and 4.output the files
     float start = clock();
     SimulateData(&programOptions, CloneNameBegin, CloneSampleSizeBegin, ClonePopSizeBegin,
                    populations,
                     &filePaths,
                     &files,
                     ObservedCellNames,
                     freq,
                      Mij
                     );
    
    // 5. deallocate the memory
    
    Population *pop;
    for( i = 0 ; i < programOptions.numClones; i++)
    {
        pop=*(populations + i);
        free(pop->idsActiveGametes);
        pop->idsActiveGametes =NULL;
        free(pop->CoalescentEventTimes);
        pop->CoalescentEventTimes = NULL;
        for (j = 1; j < ( pop->order); j++)
        {
            if (pop->order >0){
                //free(pop->immigrantsPopOrderedModelTime[j]);
                pop->immigrantsPopOrderedModelTime[j]=NULL;
            }
        }
        free( pop->migrationTimes);
        pop->migrationTimes=NULL;
        if (pop->order >0)
        {
            free(pop->immigrantsPopOrderedModelTime);
            pop->immigrantsPopOrderedModelTime=NULL;
        }
    }
    free(populations);
    populations = NULL;
    free(CloneNameBegin );
    CloneNameBegin =NULL;
    free(CloneSampleSizeBegin);
    CloneSampleSizeBegin=NULL;
    free(ClonePopSizeBegin);
    ClonePopSizeBegin=NULL;
    free(CloneBirthRateBegin );
    CloneBirthRateBegin=NULL;
    free(CloneDeathRateBegin );
    CloneDeathRateBegin=NULL;
    free(CloneTimeOriginInput);
    CloneTimeOriginInput=NULL;
    
    if(programOptions.doPrintSeparateReplicates == NO){
        fprintf(stderr, "\n\nOutput files are in folder \"Results\":");
        if (programOptions.doPrintTrees == YES  )
        {
            fprintf(stderr, "\n Trees printed to files \"%s\"", filePaths.treeFile);
            fclose(files.fpTrees);
            fclose(files.fpTrees2);
        }
        if (programOptions.doPrintTimes == YES)
        {
            fprintf(stderr, "\n Times printed to files  \"%s\"", filePaths.timesFile);
            fclose(files.fpTimes);
            fclose(files.fpTimes2);
        }
    }
    fprintf(stderr, "\n\n*** Simulations finished ***");
    /* ejecution time */
    double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    fprintf(stderr, "\n\n_________________________________________________________________");
    fprintf(stderr, "\nTime processing: %G seconds\n", secs);
    fprintf(stderr, "\nIf you need help type '-?' in the command line of the program\n");
    return 0;
    return 0;
}

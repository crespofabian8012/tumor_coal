
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
    
    // evolution model: JC
    
    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    char *ObservedCellNames[programOptions.numCells];
    
    int      *CloneNameBegin, *CloneSampleSizeBegin, *ClonePopSizeBegin;
    double   *CloneBirthRateBegin, *CloneDeathRateBegin, *CloneTimeOriginInput;
    double   *CloneGrowthRateBegin;
    double Mij[4][4];
     double freq[4];

    double TotalBirthRate, TotalDeathRate;
    int TotalN, i, j, k;

    FILE    *input_file;
    size_t seed = std::stol(argv[0]);
    string output_path = argv[1];
    string input_path = argv[2];
    
    // 1. call function to parse the input file
    if (argc <= 2)
    {
        if ((input_file = freopen("input_path", "r", stdin)) != NULL)
        {
            ReadParametersFromFile(&programOptions, &filePaths, &CloneNameBegin, &CloneSampleSizeBegin, &ClonePopSizeBegin, &CloneBirthRateBegin, &CloneDeathRateBegin, &CloneTimeOriginInput,Mij,freq);
        }
        else
        {
            fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
            PrintUsage();
        }
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
    
    
    // 3. call function to simulate the data
    
    //SimulateData();
    
    
    // 4. output the files
    
    // 5. deallocate the memory
    
    
    return 0;
}

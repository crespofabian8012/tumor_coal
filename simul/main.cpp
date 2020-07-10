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

/*
 * simulator 
 */
#include <iostream>
#include <string>
#include <vector>

#include <libpll/pll_tree.h>

#include "data_types.hpp"
#include "data_utils.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    // default evolution model: JC
    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    
    vector<int> CloneNameBegin, CloneSampleSizeBegin, ClonePopSizeBegin;
    vector<double> CloneBirthRateBegin, CloneDeathRateBegin, CloneTimeOriginInput;
    double freq[4]; // stationary distribution of CTMC over {A,C,G,T}
    double Mij[4][4]; // transition probability matrix of CTMC over {A,C,G,T}
    double Eij[4][4]; // error probability

    double TotalBirthRate, TotalDeathRate;
    int TotalN;

    FILE    *input_file;
    /* Default settings */
    Initialize( Eij, Mij, freq,  programOptions );
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
        ReadParametersFromFile(programOptions, filePaths, CloneNameBegin, CloneSampleSizeBegin, ClonePopSizeBegin, CloneBirthRateBegin, CloneDeathRateBegin, CloneTimeOriginInput,Mij,freq);
    }
    else
    {
        fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
        PrintUsage();
    }

    // 2. create and initialize data structures
    
    vector<Population *> populations;
    InitListClones(populations, programOptions.numClones, programOptions.noisy, CloneNameBegin, CloneSampleSizeBegin, CloneBirthRateBegin,  CloneDeathRateBegin, ClonePopSizeBegin, CloneTimeOriginInput, programOptions.TotalNumSequences);
    InitNumberNodes(TotalBirthRate, TotalDeathRate, TotalN, populations, programOptions);
    ListClonesAccordingTimeToOrigin(populations, programOptions.numClones);

    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    
    // 3. call function to simulate the data and
    // 4.output the files
    float start = clock();
    SimulateData(programOptions, CloneNameBegin, CloneSampleSizeBegin, ClonePopSizeBegin,
                    populations,
                    filePaths,
                    files,
                    freq,
                    Mij );
    
    // 5. deallocate the memory
    for (auto ptr : populations)
    {
        delete ptr;
    }
    populations.clear();

    // 6. output messages
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
    /* execution time */
    double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    fprintf(stderr, "\n\n_________________________________________________________________");
    fprintf(stderr, "\nTime processing: %G seconds\n", secs);
    fprintf(stderr, "\nIf you need help type '-?' in the command line of the program\n");
    return 0;
}

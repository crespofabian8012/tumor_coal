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

//#include <libpll/pll_tree.h>
extern "C"
{
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
#include <libpll/pll.h>

}

#include "data_types.hpp"
#include "data_utils.hpp"
#include "output_functions.hpp"
#include "random.h"

#include <boost/test/unit_test.hpp>
#include <boost/program_options.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>

#include <boost/random.hpp>
using namespace std;
using namespace boost::program_options;



int main(int argc, char *argv[])
{
    // default evolution model: JC
    ProgramOptions programOptions;
    Files files;
    FilePaths filePaths;
    namespace po = boost::program_options;
    
    
    po::options_description global("Global options");
    global.add_options()
    ("help", "BDCancerCoal.1.0")
    ("n", po::value<int>(), "number of replicates[mandatory]")
    ("s", po::value<int>(), "number of sites[mandatory]")
    ("command", po::value<std::string>(), "command to execute")
    ("subargs", po::value<std::vector<std::string> >(), "Arguments for command");

    
    po::positional_options_description pos;
    pos.add("command", 1).
        add("subargs", -1);

    po::variables_map vm;

    po::parsed_options parsed = po::command_line_parser(argc, argv).
        options(global).
        positional(pos).
        allow_unregistered().
        run();

    po::store(parsed, vm);
     
    
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
        Output::PrintUsage();
    }
    
    // 1. call function to parse the input file

    if ((input_file = freopen(input_path, "r", stdin)) != NULL)
    {
        ReadParametersFromFile(programOptions, filePaths, CloneNameBegin, CloneSampleSizeBegin, ClonePopSizeBegin, CloneBirthRateBegin, CloneDeathRateBegin, CloneTimeOriginInput,Mij,freq);
    }
    else
    {
        fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
        Output::PrintUsage();
    }
    
    //allocate random generators
    std::vector<gsl_rng *> rngGslvector(programOptions.numDataSets);
      // std::vector<const gsl_rng * > randomGenerators(programOptions.numDataSets);
    Random::allocateListRandomNumbersGenerators(rngGslvector);
    std::vector<boost::mt19937 *> rngBoostvector = Random::allocateListRandomNumbersGeneratorsBoost(programOptions.numDataSets);

    // 2. create and initialize data structures
    
    vector<Population *> populations;
    InitListClones(populations, programOptions.numClones, programOptions.noisy, CloneNameBegin, CloneSampleSizeBegin, CloneBirthRateBegin,  CloneDeathRateBegin, ClonePopSizeBegin, CloneTimeOriginInput, programOptions.TotalNumSequences, NO, programOptions,  rngGslvector);
    InitNumberNodes(TotalBirthRate, TotalDeathRate, TotalN, populations, programOptions);
    ListClonesAccordingTimeToOrigin(populations, programOptions.numClones);

    /* set file dirs and names */
    InitFilesPathsOptions(filePaths, programOptions);
    InitFiles(files);
   
    // 3. call function to simulate the data and
    // 4.output the files
    float start = clock();
    SimulateData(programOptions, CloneNameBegin, CloneSampleSizeBegin, ClonePopSizeBegin,
                    populations,
                    filePaths,
                    files,
                    freq,
                    Mij,
                    Eij,
                    rngGslvector,
                    rngBoostvector);
    
    // 5. deallocate the memory
    for (auto ptr : populations)
    {
        delete ptr;
    }
    populations.clear();

    // 6. output messages
    if(programOptions.doPrintSeparateReplicates == NO){
        if (programOptions.doPrintTrees == YES  )
        {
            fclose(files.fpTrees->f);
            fclose(files.fpTrees2->f);
        }
        if (programOptions.doPrintTimes == YES)
        {
            fclose(files.fpTimes->f);
            fclose(files.fpTimes2->f);
        }
    }
    std::cout << "\n\nOutput files are in folder \"Results\":"<< "\""<< std::endl;;
    std::cout << "\n Trees saved to folder \"%s\"" << filePaths.treeDir<< "\""<< std::endl;
    std::cout << "\n Times saved to folder  \"%s\"" << filePaths.timesDir<< "\""<< std::endl;
    std::cout << "\n True haplotypes(IUPAC codes) saved to folder  \"%s\"" << filePaths.trueHaplotypesDir << "\""<< std::endl;
    std::cout << "\n True genotypes saved to folder  \"" << filePaths.trueHaplotypesDir<< "\""<< std::endl;
    std::cout << "\n Full haplotypes(IUPAC codes) with errors saved to folder  \"" << filePaths.fullHaplotypesDir<< "\""<<std::endl;
    
    std::cout <<  "\n\n*** Simulations finished ***"<< std::endl;
    /* execution time */
    double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    std::cout << "\n\n_________________________________________________________________" << std::endl;
    std::cout << "\nTime processing: %G seconds\n" <<  secs<< std::endl;
    std::cout << "\nIf you need help type '-?' in the command line of the program\n" << std::endl;
    return 0;
}

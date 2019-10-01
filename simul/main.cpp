
#include <boost/program_options.hpp>

#include <iostream>
#include <string>

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
    
    size_t seed;
    string output_path;
    string input_path;
    
    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("input,i", po::value<string>(&input_path)->required(), "path to input file.")
    ("output,o", po::value<string>(&output_path)->required(), "path to output.")
    ("seed", po::value<size_t>(&seed)->default_value(1), "random seed.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    
    // 1. call function to parse the input file
    // 2. call function to simulate the data
    // 3. output the files
    
    return 0;
}

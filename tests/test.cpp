//
//  test.cpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 01/09/2021.
//
#include "hcm_test.hpp"
#include "smc_test.hpp"

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <stdio.h>
using namespace std;

int main()
{
    std::cout << "Running tests...!" << std::endl;
    //test_hmc();
    test_smc();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

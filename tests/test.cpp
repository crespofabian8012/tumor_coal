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
#include "gtest/gtest.h"
int main(int argc, char **argv)
{

    ::testing::InitGoogleTest(&argc, argv);
    int res = RUN_ALL_TESTS();
    //test_hmc();
    //test_smc();
   return res;
}

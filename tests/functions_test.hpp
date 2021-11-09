//
//  functions_test.hpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 08/09/2021.
//

#ifndef functions_test_hpp
#define functions_test_hpp


extern "C"
    {
#include <gsl/gsl_rng.h>
    }
#include "gtest/gtest.h"
//#include "gtest/internal/gtest-internal.h"
class Population;
class PopulationTest : public ::testing::Test {
protected:
    PopulationTest();
    ~PopulationTest() override;
    
    void SetUp() override;
    Population *pop;
    double K= 0.8;
    gsl_rng * random;
};


#endif

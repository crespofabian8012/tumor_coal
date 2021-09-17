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
  // You can remove any or all of the following functions if their bodies would
  // be empty.

    PopulationTest();

    ~PopulationTest() override;
     // You can do clean-up work that doesn't throw exceptions here.
  

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

    void SetUp() override;

 // void TearDown() override {
      
 // }

    Population *pop;
    double K= 0.8;
    gsl_rng * random;
};

//void  test_functions();
#endif /* functions_test_hpp */

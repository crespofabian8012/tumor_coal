//
//  functions_test.hpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 08/09/2021.
//

#ifndef functions_test_hpp
#define functions_test_hpp

#include "population.hpp"
#include "random.h"
#include "gtest/gtest.h"
class PopulationTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if their bodies would
  // be empty.

  PopulationTest() {
     // You can do set-up work for each test here.
  }

  ~PopulationTest() override {
     // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  void SetUp() override {
     
     int ord = 0;
     long double  timeOriginInput = 0.0;
     int sampleSize = 0;
     int totalSampleSize =0;
     int totalPopSize =0;
     int popSize =0;
     long double  birthRate = 1.0;
     long double  deathRate = 0.99;
     
     pop = new Population(0, ord, timeOriginInput, sampleSize,totalSampleSize,  popSize,totalPopSize, birthRate, deathRate, NO);
      
      double delta = 100;
      double tOriginModelTime = 0.52;
      sampleSize = 20;
   
      pop->delta =delta;
      pop->sampleSize = sampleSize;
      pop->timeOriginSTD = tOriginModelTime;
      long seed = 14263784;
      random = Random::generateRandomObject( seed);
  }

 // void TearDown() override {
      
 // }

    Population *pop;
    double K= 0.8;
    gsl_rng * random;
};

//void  test_functions();
#endif /* functions_test_hpp */

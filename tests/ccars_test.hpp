//
//  ccars_test.hpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 20/10/2021.
//

#ifndef ccars_test_hpp
#define ccars_test_hpp


#include "ccars.hpp"
extern "C"
    {
#include <gsl/gsl_rng.h>
    }
#include "gtest/gtest.h"
//#include "gtest/internal/gtest-internal.h"
class CCLogDensity;
class CCARS;
class CCARSTest : public ::testing::Test {
 protected:


    CCARSTest();
    ~CCARSTest() override;

    void SetUp() override;
    
    CCLogDensity *log_density;
    CCARS *ccars ;
    gsl_rng * random;
};

#endif /* ccars_test_hpp */

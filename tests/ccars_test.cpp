//
//  ccars_test.cpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 20/10/2021.
//

#include "ccars_test.hpp"
#include "ccars.cpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-spi.h"
using testing::ElementsAre;
CCARSTest::CCARSTest(){
    
    SetUp();
    
}

void CCARSTest::SetUp(){
    
    double delta = 100.0;
    double theta = 1;
    double  currModeltime = 0.0;
    double torigin = 0.039103276999;
    std::vector<double> oneMinusConcave ={0.3};
    std::vector<double> oneMinusConvex ={-0.5};
    
    
    double right_bound = 0.999*torigin;
    double interPoint1 = 0.0120;
    double interPoint2 = 0.0256;
    std::vector<double>  x = {currModeltime, interPoint1,interPoint2, right_bound};
    
    
    Eigen::VectorXd xEigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());
    
    int max_points = 10;
    long seed = 14263784;
    random = Random::generateRandomObject( seed);
    
    log_density =
    new GenotypeJCPairLogProposal(2, 16, torigin, delta, theta, currModeltime,
                                  0.0,
                                  0.0, oneMinusConcave, oneMinusConvex);
    
    ccars = new CCARS(log_density, xEigen, currModeltime, right_bound,  max_points);
    
    
}
CCARSTest::~CCARSTest(){
    
    
}
TEST_F(CCARSTest, checkAllPointsAreIncluded) {
    
    EXPECT_EQ(ccars->upper_hull.x.array().size(),5);
    
    EXPECT_THAT(ccars->upper_hull.x.array(), ElementsAre(0.0, 0.0120, 0.021009046623829714,0.0256,0.039064173722001001));
}
TEST_F(CCARSTest, checkHeightsPointsIncluded) {
    
    EXPECT_EQ(ccars->upper_hull.x.array().size(),5);
    EXPECT_EQ(ccars->lower_hull.x.array().size(),5);
    
    EXPECT_THAT(ccars->upper_hull.hx.array(), ElementsAre(0.086854298602573488, 0.0027216204543210609, -0.060339680459994105, -0.24418133608966253,-0.78318603893695515));
    
    EXPECT_THAT(ccars->lower_hull.hx.array(), ElementsAre(0.048726508541922953, 0.0027216204543210609, -0.12057459702887352, -0.24418133608966253,-11.429466410610463));
}
TEST_F(CCARSTest, checkSlopesPointsIncluded) {
    
    EXPECT_EQ(ccars->upper_hull.x.array().size(),5);
    EXPECT_EQ(ccars->lower_hull.x.array().size(),5);
    
    EXPECT_THAT(ccars->upper_hull.hpx.array(), ElementsAre(-7.0110565123543687, -6.9997751757119939, -40.04433078844, -40.032512501419774,-40.032512501419774));
    
    EXPECT_THAT(ccars->lower_hull.hpx.array(), ElementsAre(-3.8337406739668243, -18.160605858662478, -18.148671388936297, -830.74426292076112,-830.74426292076112));
}

TEST_F(CCARSTest, testVectorRatiosBoundedAreas) {
    
    EXPECT_EQ(ccars->ratioAreas.array().size(),4);
    
    EXPECT_THAT(ccars->ratioAreas.array(), ElementsAre(0.019086131794748606, 0.18719234905408222, -0.13787903499436588,0.57437822955348861));
}
TEST_F(CCARSTest, includeNewPointAtEnd) {
    
    double new_x = 0.0272760822140008;
    double new_hx = log_density->h_concave(new_x)+log_density->h_convex(new_x);
   
    double ratioAreasBefore = exp(logSumExp(ccars->lower_hull.cumulative)-logSumExp(ccars->upper_hull.cumulative));
    ccars->addPoint(new_x);
    
    EXPECT_EQ(ccars->upper_hull.x.array().size(),7);
    
    EXPECT_EQ(ccars->lower_hull.x.array().size(),7);
    
    EXPECT_THAT(ccars->upper_hull.x.array(), ElementsAre(0.0, 0.0120, 0.021009046623829714,0.0256, 0.026487195324058948, new_x, 0.039064173722001001));
    
    EXPECT_THAT(ccars->lower_hull.x.array(), ElementsAre(0.0, 0.0120, 0.018789212785232226,0.0256, 0.02643787413304087, new_x, 0.039064173722001001));
    
    EXPECT_THAT(ccars->upper_hull.hx.array(), ElementsAre(0.086854298602573488, 0.0027216204543210609, -0.060339680459994105, -0.24418133608966253, -0.27970254433950076, new_hx,-0.93903768156798639));
    
    EXPECT_THAT(ccars->lower_hull.hx.array(), ElementsAre(0.048726508541922953, 0.0027216204543210609, -0.12057459702887352, -0.24418133608966253,-0.28261535115917047, new_hx, -11.429447939372604));
    
    EXPECT_THAT(ccars->upper_hull.hpx.array(), ElementsAre(-7.0110565123543687, -6.9997751757119939, -40.04433078844, -40.037641415114351,-52.429460634310352, -52.423602469318844, -52.423602469318844));
    
    EXPECT_THAT(ccars->lower_hull.hpx.array(), ElementsAre(-3.8337406739668243, -18.160605858662478, -18.148671388936297, -45.870869566077403,-45.869406648321572, -942.33951894189704,-942.33951894189704));
    
    EXPECT_THAT(ccars->ratioAreas.array(), ElementsAre(0.019086131794748606, 0.18719234905408222, -0.13787903499436588,0.0030979766534886197, -0.0028372623437868255,0.44677597881797565 ));
   
    double ratioAreasAfter = exp(logSumExp(ccars->lower_hull.cumulative)-logSumExp(ccars->upper_hull.cumulative));
    EXPECT_EQ(ratioAreasAfter, 0.80333808289619424);
    EXPECT_GT(ratioAreasAfter, ratioAreasBefore);
}
TEST_F(CCARSTest, includeNewPointAtBegining) {
    
    double new_x = 0.006578;
   // double new_hx = log_density->h_concave(new_x)+log_density->h_convex(new_x);
 
    double ratioAreasBefore = exp(logSumExp(ccars->lower_hull.cumulative)-logSumExp(ccars->upper_hull.cumulative));
    ccars->addPoint(new_x);
    
    EXPECT_EQ(ccars->upper_hull.x.array().size(),7);
    
    EXPECT_EQ(ccars->lower_hull.x.array().size(),7);
    
    
    double ratioAreasAfter = exp(logSumExp(ccars->lower_hull.cumulative)-logSumExp(ccars->upper_hull.cumulative));
    EXPECT_EQ(ratioAreasAfter, 0.76343761060252147);
    EXPECT_GT(ratioAreasAfter, ratioAreasBefore);
}
TEST_F(CCARSTest, testApproximateIntegral) {

    double log_integral = ccars->approximateLogIntegral(random, 15);
    EXPECT_EQ(exp(log_integral),   0.0306713170724892895117);
    
    
}

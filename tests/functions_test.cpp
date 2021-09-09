//
//  functions_test.cpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 08/09/2021.
//

#include "functions_test.hpp"
#include "utils.hpp"

TEST_F(PopulationTest, logRelativePopulationSizeAtTOriginWithK0) {
  //relativePopulationSize at Torigin
  double logRelativePopulationSize = pop->LogCalculateH(pop->timeOriginSTD, pop->timeOriginSTD,pop->delta, 0.0);
  EXPECT_EQ(logRelativePopulationSize,DOUBLE_NEG_INF);
}

TEST_F(PopulationTest, logRelativePopulationSizeAtTOrigin) {
  //relativePopulationSize at Torigin
  double logRelativePopulationSize = pop->LogCalculateH(pop->timeOriginSTD, pop->timeOriginSTD,pop->delta, K);
  EXPECT_EQ(logRelativePopulationSize,DOUBLE_NEG_INF);
}
TEST_F(PopulationTest, logRelativePopulationSizeAtModelTimeK0) {
  //relativePopulationSize at 0<t<Torigin
  double logRelativePopulationSize = pop->LogCalculateH(pop->timeOriginSTD /2.0, pop->timeOriginSTD,pop->delta, 0.0);
   EXPECT_TRUE( logRelativePopulationSize <= 0);

}
TEST_F(PopulationTest, logRelativePopulationSizeAtModelTime) {
  //relativePopulationSize at 0<t<Torigin
  double logRelativePopulationSize = pop->LogCalculateH(pop->timeOriginSTD /2.0, pop->timeOriginSTD,pop->delta, K);
   EXPECT_TRUE( logRelativePopulationSize <= 0);

}
TEST_F(PopulationTest, proposeWaitingTimeNextCoalEvent) {
  //relativePopulationSize at 0<t<Torigin
    
    pop->currentModelTime = 0.0;
    pop->numActiveGametes = pop->sampleSize;
    double waitingTimeKingman;
    double currentTimeKingman;
    double timeNextCoalEvent;
    while (pop->numActiveGametes >1){
           currentTimeKingman = Population::FmodelTstandard (pop->currentModelTime , pop->timeOriginSTD, pop->delta,   K);//this is current time in Kingman coal time
         waitingTimeKingman= pop->proposeWaitingTimeNextCoalEvent(random, K);
           
        currentTimeKingman = currentTimeKingman+waitingTimeKingman;

            timeNextCoalEvent =   Population::GstandardTmodel(currentTimeKingman, pop->timeOriginSTD, pop->delta, K);
           
          EXPECT_TRUE(timeNextCoalEvent < pop->timeOriginSTD);
        pop->currentModelTime =timeNextCoalEvent;
        pop->numActiveGametes = pop->numActiveGametes-1;
    }
}
TEST_F(PopulationTest, logFModelInverseOfG) {

    double currentModelTimePop = 0.0 ;
    pop->currentModelTime = 0.0;
    pop->numActiveGametes = pop->sampleSize;
    double currentTimeKingman = Population::FmodelTstandard (currentModelTimePop , pop->timeOriginSTD, pop->delta,   K);//this is current time in Kingman coal time
    double currentModelTimePopFromG =   Population::GstandardTmodel(currentTimeKingman, pop->timeOriginSTD, pop->delta, K);
     EXPECT_EQ(currentModelTimePopFromG,currentModelTimePop);
   double waitingTimeKingman= pop->proposeWaitingTimeNextCoalEvent(random, K);
   currentTimeKingman = currentTimeKingman+waitingTimeKingman;

    double timeNextCoalEvent =   Population::GstandardTmodel(currentTimeKingman, pop->timeOriginSTD, pop->delta, K);
    
    double waitingTimeKingmanFromInverse = Population::FmodelTstandard (timeNextCoalEvent-pop->getCurrentModelTime() , pop->timeOriginSTD, pop->delta,   K);
    EXPECT_EQ(waitingTimeKingman,waitingTimeKingmanFromInverse);

}

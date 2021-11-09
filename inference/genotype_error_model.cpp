//
//  genotype_error_model.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 17/03/2021.
//

#include "genotype_error_model.hpp"

#include "constants.hpp"
#include "data_types.hpp"

GenotypeErrorModel::GenotypeErrorModel(std::string name, long double seqErrorRate,long double ADOErrorRate , int states){
    this->seqErrorRate= seqErrorRate;
    this->ADOErrorRate = ADOErrorRate;
    this->name = name;
    this->states= states;
    undefinedState = ((pll_state_t) 1 << states) - 1 ;
}
GenotypeErrorModel::GenotypeErrorModel(GenotypeErrorModel& original){
    this->seqErrorRate= original.seqErrorRate;
    this->ADOErrorRate = original.ADOErrorRate;
    this->name = original.name;
    this->states= original.states;
    undefinedState = original.undefinedState ;
}
long double GenotypeErrorModel::getSeqErrorRate() const{
    return seqErrorRate;
}
 long double GenotypeErrorModel::getADOErrorRate() const{
     return  ADOErrorRate;
     
 }
GenotypeErrorModel &GenotypeErrorModel::operator=(const GenotypeErrorModel &original){
    
    if (this == &original)
        return *this;
    
    seqErrorRate = original.seqErrorRate;
    ADOErrorRate = original.ADOErrorRate;
    name = original.name;
    states = original.states;
    undefinedState= original.undefinedState;
    return *this;
}

 void GenotypeErrorModel::computeStateErrorProbPT19(pll_state_t state,
                            std::vector<double>::iterator &clvp) const{
     
     unsigned int state_id = PLL_STATE_CTZ(state);
     static const double one_3 = 1. / 3.;
     static const double one_6 = 1. / 6.;
     static const double one_8 = 1. / 8.;
     static const double three_8 = 3. / 8.;
     static const double one_12 = 1. / 12.;

     double sum_lh = 0.;

     for (size_t k = 0; k < states; ++k)
     {
       if (state == undefinedState)
         clvp[k] = 1.;
       else
       {
         if (k == state_id)
         {
           /* 0 letters away */
           if (HOMO(state_id))
             clvp[k] = 1. - seqErrorRate + 0.5 * seqErrorRate * ADOErrorRate;
           else
             clvp[k] =  (1. - ADOErrorRate ) * (1. - seqErrorRate) + one_12 * seqErrorRate * ADOErrorRate;
         }
         else if (mut_dist[state_id][k] == 1)
         {
           /* 1 letter away */
           if (HOMO(k))
           {
             clvp[k] = one_12 * seqErrorRate * ADOErrorRate +
                       one_3  * (1. - ADOErrorRate) * seqErrorRate;
           }
           else
           {
             if (HOMO(state_id))
             {
               clvp[k] = 0.5 * ADOErrorRate + one_6 * seqErrorRate -
                         three_8 * seqErrorRate * ADOErrorRate;
             }
             else
             {
               clvp[k] = one_6 * seqErrorRate -
                         one_8 * seqErrorRate * ADOErrorRate;
             }
           }
         }
         else
         {
           /* 2 letters away */
           if (HOMO(state_id))
             clvp[k] = one_12 * seqErrorRate * ADOErrorRate;
           else
             clvp[k] = 0.;
         }

         sum_lh += clvp[k];
       }
     }
     
     
     
 }
double GenotypeErrorModel::computeStateErrorProbPT20(pll_state_t state,
                          std::vector<double>::iterator &clvp) const{
    
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;
    double sum_lh = 0.;
    
     if (state == undefinedState)
     {
       for (size_t k = 0; k < states; ++k)
         clvp[k] = 1.;
         
         sum_lh = 1.0 *states;
     }
     else
     {
       auto tstate = state;
       unsigned int ctz = 0;
       unsigned int state_id = 0;

       for (size_t k = 0; k < states; ++k)
         clvp[k] = 0.;

       while (tstate && (ctz = PLL_STATE_CTZ(tstate)) < states)
       {
         state_id = state_id ?  state_id + ctz + 1 : ctz;
         tstate >>= ctz + 1;
         double sum_lh = 0.;

         for (size_t k = 0; k < states; ++k)
         {
           double lh = 0.;

           if (k == state_id)
           {
             if (HOMO(state_id))
               lh = 1. - seqErrorRate + 0.5 * seqErrorRate * ADOErrorRate;
             else
               lh = 1. - seqErrorRate - ADOErrorRate + seqErrorRate * ADOErrorRate;
           }
           else if (mut_dist16[state_id][k] == 1)
           {
             if (HOMO(k))
               lh = (1. - ADOErrorRate) * seqErrorRate * one_6;
             else
             {
               if (HOMO(state_id))
                 lh = 0.5 * ADOErrorRate + one_6 * seqErrorRate -
                     one_3 * seqErrorRate * ADOErrorRate;
               else
                 lh = (1. - ADOErrorRate) * seqErrorRate * one_6;
             }
           }
           else if (HOMO(state_id))
             lh = one_6 * seqErrorRate * ADOErrorRate;
           else
             lh = 0.;

           clvp[k] += lh;

           sum_lh += lh;
         }
        assert(sum_lh>0);
        //std::cout << sum_lh << std::endl;
       }
     }
    return sum_lh;
}

void GenotypeErrorModel::computeStateErrorProbPT17(pll_state_t state,
                           std::vector<double>::iterator &clvp) const{
    
    unsigned int state_id = PLL_STATE_CTZ(state);
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;

    double sum_lh = 0.;

    for (size_t k = 0; k < states; ++k)
    {
      if (state == undefinedState)
        clvp[k] = 1.;
      else
      {
        if (k == state_id)
        {
          if (HOMO(state_id))
            clvp[k] = 1. - seqErrorRate + 0.5 * seqErrorRate * ADOErrorRate;
          else
            clvp[k] = 1. - seqErrorRate - ADOErrorRate + seqErrorRate * ADOErrorRate;
        }
        else if (mut_dist[state_id][k] == 1)
        {
          if (HOMO(k))
            clvp[k] = (1. - ADOErrorRate) * seqErrorRate * one_3;
          else
          {
            if (HOMO(state_id))
              clvp[k] = 0.5 * ADOErrorRate + one_6 * seqErrorRate -
                  one_3 * seqErrorRate * ADOErrorRate;
            else
              clvp[k] = (1. - ADOErrorRate) * seqErrorRate * one_6;
          }
        }
        else if (HOMO(state_id))
          clvp[k] = one_6 * seqErrorRate * ADOErrorRate;
        else
          clvp[k] = 0.;

        sum_lh += clvp[k];
      }
    }
    assert(sum_lh>0);
   // std::cout << sum_lh << std::endl;

}

void GenotypeErrorModel::computeStateErrorProbPT19(pll_state_t state,
                            double *clvp) const{
     
     unsigned int state_id = PLL_STATE_CTZ(state);
     static const double one_3 = 1. / 3.;
     static const double one_6 = 1. / 6.;
     static const double one_8 = 1. / 8.;
     static const double three_8 = 3. / 8.;
     static const double one_12 = 1. / 12.;

     double sum_lh = 0.;

     for (size_t k = 0; k < states; ++k)
     {
       if (state == undefinedState)
         clvp[k] = 1.;
       else
       {
         if (k == state_id)
         {
           /* 0 letters away */
           if (HOMO(state_id))
             clvp[k] = 1. - seqErrorRate + 0.5 * seqErrorRate * ADOErrorRate;
           else
             clvp[k] =  (1. - ADOErrorRate ) * (1. - seqErrorRate) + one_12 * seqErrorRate * ADOErrorRate;
         }
         else if (mut_dist[state_id][k] == 1)
         {
           /* 1 letter away */
           if (HOMO(k))
           {
             clvp[k] = one_12 * seqErrorRate * ADOErrorRate +
                       one_3  * (1. - ADOErrorRate) * seqErrorRate;
           }
           else
           {
             if (HOMO(state_id))
             {
               clvp[k] = 0.5 * ADOErrorRate + one_6 * seqErrorRate -
                         three_8 * seqErrorRate * ADOErrorRate;
             }
             else
             {
               clvp[k] = one_6 * seqErrorRate -
                         one_8 * seqErrorRate * ADOErrorRate;
             }
           }
         }
         else
         {
           /* 2 letters away */
           if (HOMO(state_id))
             clvp[k] = one_12 * seqErrorRate * ADOErrorRate;
           else
             clvp[k] = 0.;
         }

         sum_lh += clvp[k];
       }
     }
     
     
     
 }
double GenotypeErrorModel::computeStateErrorProbPT20(pll_state_t state,
                           double *clvp) const{
    
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;
    double sum_lh = 0.;

     if (state == undefinedState)
     {
       for (size_t k = 0; k < states; ++k)
         clvp[k] = 1.;
         
         sum_lh = 1.0* states;
     }
     else
     {
       auto tstate = state;
       unsigned int ctz = 0;
       unsigned int state_id = 0;

       for (size_t k = 0; k < states; ++k)
         clvp[k] = 0.;

    
    while (tstate && (ctz = PLL_STATE_CTZ(tstate)) < states)
       {
         state_id = state_id ?  state_id + ctz + 1 : ctz;
        
         //std::cout<< "state_id" << state_id << std::endl;
         tstate >>= ctz + 1;
         sum_lh = 0.;

         for (size_t k = 0; k < states; ++k)
         {
          double lh = 0.;

           if (k == state_id)
           {
             if (HOMO(state_id))
               lh = 1. - seqErrorRate + 0.5 * seqErrorRate * ADOErrorRate;
             else
               lh = 1. - seqErrorRate - ADOErrorRate + seqErrorRate * ADOErrorRate;
           }
           else if (mut_dist16[state_id][k] == 1)
           {
             if (HOMO(k))
               lh = (1. - ADOErrorRate) * seqErrorRate * one_6;
             else
             {
               if (HOMO(state_id))
                 lh = 0.5 * ADOErrorRate + one_6 * seqErrorRate -
                     one_3 * seqErrorRate * ADOErrorRate;
               else
                 lh = (1. - ADOErrorRate) * seqErrorRate * one_6;
             }
           }
           else if (HOMO(state_id))
             lh = one_6 * seqErrorRate * ADOErrorRate;
           else
             lh = 0.;

           clvp[k] += lh;

           sum_lh += lh;
         }
     
        assert(sum_lh>0);
        //std::cout << sum_lh << std::endl;
        
       }
     }
//     printf("sumlh: %lf\n", sum_lh);
//     if (state != undefinedState)
//     {
//      for (size_t k = 0; k < states; ++k)
//           clvp[k] /= sum_lh;
    return sum_lh;
    
}

void GenotypeErrorModel::computeStateErrorProbPT17(pll_state_t state,
                            double *clvp) const{
    
    unsigned int state_id = PLL_STATE_CTZ(state);//id of observed state
    static const double one_3 = 1. / 3.;
    static const double one_6 = 1. / 6.;

    double sum_lh = 0.;

    for (size_t k = 0; k < states; ++k)
    {
      if (state == undefinedState)
        clvp[k] = 1.;
      else
      {
        if (k == state_id)
        {
          if (HOMO(state_id))
            clvp[k] = 1. - seqErrorRate + 0.5 * seqErrorRate * ADOErrorRate;
          else
            clvp[k] = 1. - seqErrorRate - ADOErrorRate + seqErrorRate * ADOErrorRate;
        }
        else if (mut_dist[state_id][k] == 1)
        {
          if (HOMO(k))
            clvp[k] = (1. - ADOErrorRate) * seqErrorRate * one_3;
          else
          {
            if (HOMO(state_id))
              clvp[k] = 0.5 * ADOErrorRate + one_6 * seqErrorRate -
                  one_3 * seqErrorRate * ADOErrorRate;
            else
              clvp[k] = (1. - ADOErrorRate) * seqErrorRate * one_6;
          }
        }
        else if (HOMO(state_id))
          clvp[k] = one_6 * seqErrorRate * ADOErrorRate;
        else
          clvp[k] = 0.;

        sum_lh += clvp[k];
      }
    }
    assert(sum_lh>0);
    //std::cout << sum_lh << std::endl;

}


//
//  struct_coal_model.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#include "bd_coal_model.hpp"
BDCoalModel::BDCoalModel(MSA *inputMSA):inputMSA(inputMSA)
{
}
unsigned long BDCoalModel::num_iterations(){
    unsigned long res =inputMSA->getSize();
    return res;
}
shared_ptr<BDCoalState> BDCoalModel::propose_initial(gsl_rng *random, double &log_w, BDCoalPriorParams &params){
    
    shared_ptr<BDCoalState> result;
    
    
    
    
    return result;
}
shared_ptr<BDCoalState> BDCoalModel::propose_next(gsl_rng *random, unsigned int t, const BDCoalState &curr, double &log_w, BDCoalPriorParams &params){
    
    std::shared_ptr<BDCoalState> result(make_shared<BDCoalState>(std::move(curr)));
   
    
    
    return result;
    
}
 double BDCoalModel::log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<BDCoalState> > &genealogy, const BDCoalPriorParams &params){
     double log_weight = 0.0;
     
     return log_weight;
 }
void BDCoalModel::generate_data(gsl_rng *random, size_t T, BDCoalPriorParams &params, vector<double> &latent, vector<double> &obs){
     
     
     
     
     
}

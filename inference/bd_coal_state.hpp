//
//  bd_coal_state.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 05/10/2021.
//

#ifndef bd_coal_state_hpp
#define bd_coal_state_hpp

#include <vector>

using namespace  std;
class BDCoalState{
    
    size_t num_populations;
    
    std::vector<double> scaled_growth_rates;
    std::vector<double> time_origins;
    std::vector<double> proportions_vector;
    double theta;
    
public:
    BDCoalState(size_t num_populations);
    
//    getters
    double getNumberPopulations() const{return num_populations;};
    double getTheta() const{return theta;};
    std::vector<double> getScaledGrowthRates() const{return scaled_growth_rates;};
    double getScaledGrowthRatesAt(size_t i) const{return scaled_growth_rates[i];};
    std::vector<double> getTimeOrigins() const{return time_origins;};
    double getTimeOriginAt(size_t i) const{return time_origins[i];};
    std::vector<double> getProportionsVector() const{return proportions_vector;};
    double getProportionPopulation(size_t i) const{return proportions_vector[i];};
};
#endif /* bd_coal_state_hpp */

//
//  ccars.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 12/10/2021.
//

#ifndef ccars_hpp
#define ccars_hpp

#include "data_types.hpp"


class GenotypeJCPairLogProposal : public CCLogDensity {
    
public:
    int numSites;
    int numStates;
    double Torigin;
    double delta;
    double theta;
    double pair_creation_time;
    double time_left_child;
    double time_right_child;
    double* left_clv;
    double* right_clv;
    double K = 0.8;
    std::vector<double> oneMinusSumTermConcave;
    std::vector<double> oneMinusSumTermConvex;
    
    GenotypeJCPairLogProposal( int numSites,
       int numStates,double Torigin, double delta, double theta, double pair_creation_time,
                              double time_left_child,
                              double time_right_child, double* left_clv, double* right_clv);
    
    
    std::vector<double> clean(std::vector<double> x);
    
    std::vector<double> h_concave(std::vector<double> x);
    std::vector<double> h_convex(std::vector<double> x);
    std::vector<double> h_prime_concave(std::vector<double> x);
    std::vector<double> h_prime_convex(std::vector<double> x);
    double h_concave(double x);
    double h_convex(double x);
    double h_prime_concave(double x);
    double h_prime_convex(double x);
};

#endif /* ccars_hpp */

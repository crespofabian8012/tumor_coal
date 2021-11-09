//
//  hcm_test.cpp
//  tests
//
//  Created by Fausto Fabian Crespo Fernandez on 31/08/2021.
//


#include "hcm_test.hpp"
#include "hmc_nuts.hpp"

#include <iostream>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <numeric>


void  test_hmc(){
    
    int K = 2;
    int F = 100;
    
    double alpha = 1.;
    double beta  = 1.;
    
    std::default_random_engine generator;
    std::gamma_distribution<float> mygamma(alpha, 1/beta);
    
    
    const Eigen::MatrixXd W = (Eigen::ArrayXXd::Random(F, K)*0.5+0.5)*1.0; //Eigen::MatrixXd::Random(F, K);//random values in [-1,1]
    //W = (W + Eigen::MatrixXd::Constant(F,K,1.)) * 1./2.;//add 1.0 to W to get values in //[0,2] and then divide by 2 to get values in  [0,1]
    Eigen::VectorXi v_n(F);
    Eigen::VectorXd h_n(K);
    h_n = h_n.setConstant(1)*5;
    
    std::cout << "h_n:" << std::endl;
    std::cout << h_n << std::endl;
    
    for(int f=0; f<F; f++){
        float lambda_f = W.row(f)*h_n;
        assert(lambda_f>0);
        std::cout << "lambda " << f << ": " << lambda_f << std::endl;
        std::poisson_distribution<int> mypoisson(lambda_f);
        v_n[f] = mypoisson(generator);
    }
    
    std::cout << "Observation:" << std::endl;
    std::cout << v_n.transpose() << std::endl;
    std::cout << "Dictionary W:" << std::endl;
    std::cout << W << std::endl;
    
    float epsilon = 0.0001;
    int iter = 10000;
    h_n.array()+= 20.f ;
    
    //    Eigen::MatrixXd runNutsSampler(const Eigen::VectorXi v_n,//observed data
    //                                   const Eigen::MatrixXd& W,
    //                                   Eigen::VectorXd current_q,//current state
    //                                   double alpha = 1, double beta = 1,
    //                                   float epsilon = 0.01,
    //                                   int iter=100)
    Eigen::MatrixXd samples = runNutsSampler(v_n, W, h_n, alpha, beta, epsilon, iter);
    
    std::cout << "samples of h_n:" << std::endl;
    std::cout << samples << std::endl;
    
    std::cout << "Real h_n:" << std::endl;
    std::cout << h_n.transpose() << std::endl;
    
}

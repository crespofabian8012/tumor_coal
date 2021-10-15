//
//  model.hpp
//  libTumorCoal
//
//  Created by Fausto Fabian Crespo Fernandez on 30/09/2021.
//

#ifndef model_hpp
#define model_hpp

#include "mcmc_parameter.hpp"
#include "utils.hpp"
#include <vector>
class Parameter{
public:
    std::string name;
    bool isFixed = false;
    bool isDiscrete = false;
    bool isHierarchical = false;
    bool isCMTCRate = false;
    bool isNonNegative = false;
    bool isZeroOne = false;
    bool isConstrained = false;
    double constrainUpper = Utils::DOUBLE_INF;
    double constrainLower = Utils::DOUBLE_NEG_INF;
    
    
};
class Model{
    std::vector<Parameter> parameters;
    std::string model_name();
public:
       
       Model(size_t numParams);

       Model(Model&&) = default;  // support moving
       Model& operator=(Model&&) = default;
       Model(const Model&) = default; // support copying
       Model& operator=(const Model&) = default;
    
       virtual ~Model() {};
       virtual std::string model_name() const = 0;
    
};



class AbstractModel{
    
private:
//    std::vector<Model> models;
//    std::vector<MCMCParameter<int>> fixedIntParameters;
//    std::vector<MCMCParameter<double>> fixedDoubleParameters;
//    std::vector<MCMCParameterWithKernel> paramsWKernel;
public:
    AbstractModel();
    void addFixedIntParam(MCMCParameter<int> &param) {};
       
//        if (!fixedIntParameters.contains(param)) {
//            variables.add(variable);
//            variable.addVariableListener(this);
//        }
//
//        // parameters are also statistics
//        if (variable instanceof Statistic) addStatistic((Statistic) variable);
//    }
//
//    public  void removeFixedIntParam(MCMCParameter<int> param) {
//        variables.remove(variable);
//        variable.removeVariableListener(this);
//
//        // parameters are also statistics
//        if (variable instanceof Statistic) removeStatistic((Statistic) variable);
//    }
    
    virtual ~AbstractModel() = default; // make dtor virtual
    AbstractModel(AbstractModel&&) = default;  // support moving
    AbstractModel& operator=(AbstractModel&&) = default;
    AbstractModel(const AbstractModel&) = default; // support copying
    AbstractModel& operator=(const AbstractModel&) = default;

};
    
class Likelihood{
public:
    virtual ~Likelihood() = default;
    Likelihood(Likelihood&&) = default;  // support moving
    Likelihood& operator=(Likelihood&&) = default;
    Likelihood(const Likelihood&) = default; // support copying
    Likelihood& operator=(const Likelihood&) = default;

   
    virtual double getLikelihood() = 0;
    };

class ModelLikelihood: public AbstractModel, public Likelihood {
    
public:
    ModelLikelihood();
    
    
   ModelLikelihood(ModelLikelihood&&) = default;  // support moving
   ModelLikelihood& operator=(ModelLikelihood&&) = default;
   ModelLikelihood(const ModelLikelihood&) = default; // support copying
   ModelLikelihood& operator=(const ModelLikelihood&) = default;
    
   double getLikelihood();
   ~ModelLikelihood() {}; // make dtor virtual
};

#endif /* model_hpp */

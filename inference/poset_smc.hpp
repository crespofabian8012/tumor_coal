//
//  poset_smc.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 16/03/2021.
//

#ifndef poset_smc_hpp
#define poset_smc_hpp




#include "spf.hpp"
#include "poset_smc_params.hpp"
#include "utils.hpp"


#include <map>
#include <unordered_map>
#include <iterator>
#include <set>
#include <random>
#include <boost/thread/mutex.hpp>

class State;
class PosetSMCParams;
typedef std::pair<int, int> pairs;

using namespace std;

class PosetSMC : public ProblemSpecification<State, PosetSMCParams>
{
    static  std::unordered_map<size_t , std::set<pairs > > sizeCombinationMap;
    static boost::mutex mx;
    //mutable std::mutex m;
    size_t numClones;
    size_t numIter;
    std::vector<bool> doPlotPerIteration;

public:
    PosetSMC(size_t numClones, size_t num_iter, bool doPlots);
    unsigned long num_iterations() override;
    std::shared_ptr<State> propose_initial(gsl_rng *random, double &log_w, PosetSMCParams &params) override;
    std::shared_ptr<State> propose_next(gsl_rng *random, unsigned int t, const State &curr, double &log_w, PosetSMCParams &params) override;
    
    // static void generate_data(gsl_rng *random, size_t T, SMCOptions &params, std::vector<double> &latent, std::vector<double> &obs)override ;
    
    double log_weight(unsigned int t, const shared_ptr<ParticleGenealogy<State> > &genealogy, const PosetSMCParams &p) override;
    
    //void set_particle_population(const vector<shared_ptr<State> > &particles) override;
    enum PosetSMCKernel {
        PRIORPRIOR = 0,
        PRIORPOST = 1,
        POSTPOST1 = 2,
        POSTPOST2 = 3,
        TSMC1     = 4
    };
    
    enum PosetSMCIncrementProposal {
          UNIF = 0,
          EXP = 1,
          LBD_COAL = 2
      };
    PosetSMCIncrementProposal incrementProposalDist ;
     PosetSMCKernel kernelType = PRIORPOST ;
    size_t numIncrementsPOSTPOST = 10;
    
    std::string getPosetSMCKernelName() const{
        std::string result;
         switch(kernelType) {
           case PRIORPRIOR:
                 result = "PRIORPRIOR";
                 break;
           case PRIORPOST:
                 result = "PRIORPOST";
                 break;
           case POSTPOST1:
                 result = "POSTPOST1";
                 break;
           case POSTPOST2:
                 result = "POSTPOST2";
                 break;
          case TSMC1:
                 result = "TSMC1";
                 break;
           default:
              result = "UNDEFINED";
              break;
         }
        return result;
    }
    
    std::string getPosetSMCPriorName() const{
        std::string result;
         switch(incrementProposalDist) {
           case UNIF:
                 result = "Uniform";
                 break;
           case EXP:
                 result = "EXP";
                 break;
           case LBD_COAL:
              result = "LBD_coalescent";
               break;
           default:
              result = "UNDEFINED";
              break;
         }
        return result;
    }
    
   static  std::set<std::pair<int, int> > getCombinations(size_t size){
       
        boost::mutex::scoped_lock scoped_lock(mx);
        if (sizeCombinationMap.find(size) != sizeCombinationMap.end())
            
            return sizeCombinationMap[size];
        else
        {
            std::vector<pairs> allPairs= Utils::allCombinations(size, 2);
            std::set<pairs> s(allPairs.begin(), allPairs.end());
            sizeCombinationMap[size] = s;
        }
        return sizeCombinationMap[size];
    }

    static  std::vector<std::pair<int, int> > getCombinationsRandomOrder(size_t size){
       
       
        std::set<Pair> pairSet = getCombinations(size);
        std::vector<Pair> pairs(pairSet.begin(), pairSet.end());
        return pairs;
    }
   
    ~PosetSMC(){};
};


#endif /* poset_smc_hpp */





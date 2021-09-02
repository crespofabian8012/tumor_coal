//
//  msa.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 26/08/2021.
//

#ifndef msa_hpp
#define msa_hpp

#include "data_types.hpp"
#include <vector>
//#include <unordered_map>
extern "C"
    {
#include "libpll/pll.h"
    }

using namespace std;
typedef std::vector<double> ProbVector;
typedef std::vector<ProbVector> ProbVectorList;
class MSA
{
private:
    size_t length;
    size_t num_sites;
    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    NameIdMap label_id_map;
    std::vector<unsigned int> weights;
    std::vector<unsigned int> site_pattern_map;
    ProbVectorList probs;
    // RangeList local_seq_ranges;
    size_t states;
    mutable pll_msa_t * pll_msa;
    mutable bool dirty;
    
    void update_pll_msa() const;
    
public:
    
    MSA() : length(0), num_sites(0), states(0), pll_msa(NULL), dirty(false) {};
    MSA(const unsigned int num_sites) : length(0), num_sites(num_sites),
    states(0), pll_msa(nullptr), dirty(false) {};
    //MSA(const RangeList& rl);
    
    MSA(const pll_msa_t * pll_msa);
    MSA(MSA&& other);
    MSA(const MSA& other) = delete;
    
    ~MSA();
    
    MSA& operator=(MSA&& other);
    MSA& operator=(const MSA& other) = delete;
    const std::string& label(size_t index) const { return labels.at(index); }
    
    const pll_msa_t * getRawPtr() const{return pll_msa;};
    const std::vector<std::string>& getLabels() const { return
        labels; };
    const std::string& at(const std::string& label) const
    { return sequences.at(label_id_map.at(label)); }
    const std::string& at(size_t index) const { return sequences.at(index); }
    const std::string& operator[](const std::string& label) const { return at(label); }
    const std::string& operator[](size_t index) const { return at(index); }
    std::string& operator[](size_t index) { return sequences.at(index); }
    
    size_t getSize() const { return sequences.size(); }
    size_t getLength() const { return length; }
    
    static bool checkMSA(const MSA& msa);
    void free_pll_msa() noexcept;
    void append(const string& sequence, const string& header);
};
#endif /* msa_hpp */

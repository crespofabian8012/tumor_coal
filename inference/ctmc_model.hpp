//
//  model.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 27/09/2021.
//

#ifndef ctmc_model_hpp
#define ctmc_model_hpp


#include "definitions.hpp"
#include <iostream>
#include <limits>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

extern "C"
    {
#include "libpll/pll.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_common.h"
#include <libpll/pll_msa.h> //for pllmod_msa_empirical_frequencies
#include <libpll/pll_tree.h>
#include <libpll/pllmod_util.h>
#include "libpll/pll_optimize.h"
#include <libpll/pllmod_common.h>
#include <libpll/pllmod_algorithm.h>
    
    }

using namespace std;
typedef std::unordered_map<pll_state_t,std::string> StateNameMap;

class SubstitutionModel
{
public:
  SubstitutionModel(const pllmod_subst_model_t& sm) :
    states(sm.states), name(sm.name)
  {
    if (sm.freqs)
      base_freqs.assign(sm.freqs, sm.freqs + sm.states);
    if (sm.rates)
      subst_rates.assign(sm.rates, sm.rates + sm.states*(sm.states-1)/2);
    if (sm.rate_sym)
      rate_sym.assign(sm.rate_sym, sm.rate_sym + sm.states*(sm.states-1)/2);
    if (sm.freq_sym)
      freq_sym.assign(sm.freq_sym, sm.freq_sym + sm.states);
  };

  // getters
  unsigned int getStates() const;
  std::string getName() const;
  const std::vector<double>& baseFreqs() const { return base_freqs; }
  const std::vector<double>& substRates() const { return subst_rates; }
  const std::vector<int>& rateSym() const { return rate_sym; }
  const std::vector<int>& freqSym() const { return freq_sym; }

  unsigned int numRates() const  { return states*(states-1)/2; }
  unsigned int numUniqRates() const
  {
    if (rate_sym.empty())
      return numRates();
    else
      return *std::max_element(rate_sym.cbegin(), rate_sym.cend()) + 1;
  }

  std::vector<double> uniqSubstRates() const
  {
    if (!rate_sym.empty())
    {
      std::vector<double> uniq_rates(numUniqRates());
      for (size_t i = 0; i < subst_rates.size(); ++i)
        uniq_rates[rate_sym[i]] = subst_rates[i];
      return uniq_rates;
    }
    else
      return subst_rates;
  };


  // setters
  void baseFreqs(const std::vector<double>& v)
  {
//    std::cout << "expected: " << _states << ", got: " << v.size() << std::endl;
    if (v.size() != states)
      throw std::invalid_argument("Invalid size of base_freqs vector!");

    base_freqs = v;
  };

  void substRates(const std::vector<double>& v)
  {
    if (v.size() != numRates())
      throw std::invalid_argument("Invalid size of subst_rates vector!");

    subst_rates = v;
  };

  void uniqSubstRates(const std::vector<double>& v)
  {
    if (!rate_sym.empty())
    {
      if (v.size() != numUniqRates())
        throw std::invalid_argument("Invalid size of subst_rates vector!");

      subst_rates.resize(numRates());
      for (size_t i = 0; i < subst_rates.size(); ++i)
        subst_rates[i] = v[rate_sym[i]];
    }
    else
      substRates(v);
  };

private:
  unsigned int states;
  std::string name;
  std::vector<double> base_freqs;
  std::vector<double> subst_rates;
  std::vector<int> rate_sym;
  std::vector<int> freq_sym;
};

class Model
{
public:
  typedef std::unordered_map<int,ParamValue> ParamModeMap;

  Model (DataType data_type = DataType::autodetect, const std::string &model_string = "GTR");
  Model (const std::string &model_string) : Model(DataType::autodetect, model_string) {};

  Model(const Model&) = default;

  /* getters */
  DataType dataType() const { return data_type; };
  std::string dataTypeName() const;
  unsigned int numStates() const { return num_states; };
  std::string getName() const { return name; };

  const pll_state_t* charmap() const;
  const std::vector<std::string>& stateNames() const;            // non-ambiguous states only, eg A C G T
  const StateNameMap& fullStateNamemap() const; // + ambiguous states, eg A C G T M R W S Y K -

  const SubstitutionModel submodel(size_t i) const { return submodels.at(i); };

  unsigned int ratehetMode() const { return rate_het; };
  unsigned int numRatecats() const { return num_ratecats; };
  unsigned int numSubmodels() const { return num_submodels; };
  const std::vector<double>& ratecatRates() const { return ratecat_rates; };
  const std::vector<double>& ratecatWeights() const { return ratecat_weights; };
  const std::vector<unsigned int>& ratecatSubmodels() const { return ratecat_submodels; };
  int gammaMode() const { return gamma_mode; };

  double getAlpha() const { return alpha; };
  double getPinv() const { return pinv; };
  double brlenScaler() const { return brlen_scaler; };
  const std::vector<double>& baseFreqs(unsigned int i) const { return submodels.at(i).baseFreqs(); };
  const std::vector<double>& substRates(unsigned int i) const { return submodels.at(i).substRates(); };

  std::string toString(bool print_params = false, unsigned int precision = 0) const;
  int params_to_optimize() const;
  const ParamModeMap& paramMode() const { return param_mode; }
  ParamValue paramMode(int param) const { return param_mode.at(param); };
  bool param_estimated(int param) const;

  AscBiasCorrection ascbiasType() const { return ascbias_type; }
  const std::vector<unsigned int> & ascbiasWeights() const { return ascbias_weights; }

  /* per alignment site, given in elements (NOT in bytes) */
  size_t clv_entry_size() const { return num_states * num_ratecats; }

  unsigned int  num_free_params() const;

  /* setters */
  void setAlpha(double value) { alpha = value; };
  void setPinv(double value) { pinv = value; };
  void brlenScaler(double value) { brlen_scaler = value; };
  void baseFreqs(size_t i, const std::vector<double>& value) { submodels.at(i).baseFreqs(value); };
  void substRates(size_t i, const std::vector<double>& value) { submodels.at(i).substRates(value); };
  void baseFreqs(const std::vector<double>& value) { for (SubstitutionModel& s: submodels) s.baseFreqs(value); };
  void substRates(const std::vector<double>& value) { for (SubstitutionModel& s: submodels) s.substRates(value); };
  void ratecatRates(std::vector<double> const& value) { ratecat_rates = value; };
  void ratecatWeights(std::vector<double> const& value) { ratecat_weights = value; };

  void paramMode(int param, ParamValue mode) { param_mode[param] = mode; };
  void set_param_mode_default(int param, ParamValue mode)
  {
    if (paramMode(param) == ParamValue::undefined)
      param_mode[param] = mode;
  };

  /* initialization */
  void init_from_string(const std::string& model_string);

private:
  std::string name;
  DataType data_type;
  unsigned int num_states;

  std::string custom_states;
  std::string custom_gaps;
  bool custom_case_sensitive;
  std::shared_ptr<pll_state_t> custom_charmap;
  mutable std::vector<std::string> state_names;
  mutable StateNameMap full_state_namemap;

  unsigned int rate_het;
  unsigned int num_ratecats;
  unsigned int num_submodels;
  std::vector<double> ratecat_rates;
  std::vector<double> ratecat_weights;
  std::vector<unsigned int> ratecat_submodels;
  int gamma_mode;

  double alpha;
  double pinv;
  double brlen_scaler;

  AscBiasCorrection ascbias_type;
  std::vector<unsigned int> ascbias_weights;

  std::vector<SubstitutionModel> submodels;

  ParamModeMap param_mode;

  void autodetect_data_type(const std::string& model_name);
  pllmod_mixture_model_t * init_mix_model(const std::string& model_name);
  void init_model_opts(const std::string& model_opts, const pllmod_mixture_model_t& mix_model);
  void init_state_names() const;
  void set_user_srates(std::vector<double>& srates, bool normalize = true);
  void set_user_freqs(std::vector<double>& freqs);
};

typedef std::unordered_map<size_t, Model> ModelMap;
typedef std::unordered_map<size_t, Model&> ModelRefMap;
typedef std::unordered_map<size_t, const Model&> ModelCRefMap;

void assign(Model& model, const pll_partition_t * partition);
void assign(pll_partition_t * partition, const Model& model);



#endif /* ctmc_model_hpp */

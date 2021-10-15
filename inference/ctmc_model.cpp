//
//  model.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 27/09/2021.
//

#include "ctmc_model.hpp"
#include "data_types.hpp"

const unordered_map<DataType,unsigned int,EnumClassHash>  DATATYPE_STATES { {DataType::dna, 4},
                                                                            {DataType::protein, 20},
                                                                            {DataType::binary, 2},
                                                                            {DataType::genotype10, 10}
                                                                          };

const unordered_map<DataType,const pll_state_t*,EnumClassHash>  DATATYPE_MAPS {
  {DataType::dna, pll_map_nt},
  {DataType::protein, pll_map_aa},
  {DataType::binary, pll_map_bin},
  {DataType::genotype10, pll_map_gt10}
};

const unordered_map<DataType,string,EnumClassHash>  DATATYPE_PREFIX { {DataType::dna, "DNA"},
                                                                      {DataType::protein, "PROT"},
                                                                      {DataType::binary, "BIN"},
                                                                      {DataType::genotype10, "GT"},
                                                                      {DataType::multistate, "MULTI"},
                                                                      {DataType::autodetect, "AUTO"}
                                                                    };
void assign(Model& model, const pll_partition_t * partition)
{
  if (model.numStates() == partition->states &&
      model.numSubmodels() == partition->rate_matrices)
  {
    model.setPinv(partition->prop_invar[0]);
    for (size_t i = 0; i < model.numSubmodels(); ++i)
    {
      model.baseFreqs(i, std::vector<double>(partition->frequencies[i],
                                       partition->frequencies[i] + partition->states));

      size_t n_subst_rates = pllmod_util_subst_rate_count(partition->states);
      model.substRates(i, std::vector<double>(partition->subst_params[i],
                                        partition->subst_params[i] + n_subst_rates));
    }

    if (partition->rate_cats > 1)
    {
      model.ratecatRates(std::vector<double>(partition->rates,
                                       partition->rates + partition->rate_cats));
      model.ratecatWeights(std::vector<double>(partition->rate_weights,
                                         partition->rate_weights + partition->rate_cats));
    }
  }
  else
    throw runtime_error("incompatible partition!");
}

void assign(pll_partition_t * partition, const Model& model)
{
  if (model.numStates() == partition->states &&
      model.numSubmodels() == partition->rate_matrices)
  {
    /* set rate categories & weights */
    pll_set_category_rates(partition, model.ratecatRates().data());
    pll_set_category_weights(partition, model.ratecatWeights().data());

    /* now iterate over rate matrices and set all params */
    for (size_t i = 0; i < partition->rate_matrices; ++i)
    {
      /* set base frequencies */
      assert(!model.baseFreqs(i).empty());
      pll_set_frequencies(partition, i, model.baseFreqs(i).data());

      /* set substitution rates */
      assert(!model.substRates(i).empty());
      pll_set_subst_params(partition, i, model.substRates(i).data());

      /* set p-inv value */
      pll_update_invariant_sites_proportion (partition, i, model.getPinv());
    }
  }
  else
    throw runtime_error("incompatible partition!");
}

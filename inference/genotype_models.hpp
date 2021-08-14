//
//  genotype_models.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 09/08/2021.
//

#ifndef genotype_models_hpp
#define genotype_models_hpp

#include <stdio.h>

#include <string.h>

extern "C"
    {
#include <libpll/pllmod_util.h>
#include <libpll/pllmod_common.h>
    }


/*                                       AA CC GG TT AC AG AT CG CT GT          */
static const double gt_rates_equal_sm[] = {  0, 0, 0, 1, 1, 1, 0, 0, 0,    /* AA */
                                                0, 0, 1, 0, 0, 1, 1, 0,    /* CC */
                                                   0, 0, 1, 0, 1, 0, 1,    /* GG */
                                                      0, 0, 1, 0, 1, 1,    /* TT */
                                                         1, 1, 1, 1, 0,    /* AC */
                                                            1, 1, 0, 1,    /* AG */
                                                               0, 1, 1,    /* AT */
                                                                  1, 1,    /* CG */
                                                                     1 };  /* CT */

/*                                      AA CC GG TT AC AG AT CG CT GT          */
static const double gt_rates_equal[] =  {   1, 1, 1, 1, 1, 1, 1, 1, 1,    /* AA */
                                               1, 1, 1, 1, 1, 1, 1, 1,    /* CC */
                                                  1, 1, 1, 1, 1, 1, 1,    /* GG */
                                                     1, 1, 1, 1, 1, 1,    /* TT */
                                                        1, 1, 1, 1, 1,    /* AC */
                                                           1, 1, 1, 1,    /* AG */
                                                              1, 1, 1,    /* AT */
                                                                 1, 1,    /* CG */
                                                                    1 };  /* CT */

/*                                    AA CC GG TT AC AG AT CG CT GT CA GA TA GC TC TG          */
static const double gt16_rates_equal[120] =
                                    {     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AA */
                                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* CC */
                                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* GG */
                                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* TT */
                                                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AC */
                                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AG */
                                                            1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AT */
                                                               1, 1, 1, 1, 1, 1, 1, 1,   /* CG */
                                                                  1, 1, 1, 1, 1, 1, 1,   /* CT */
                                                                     1, 1, 1, 1, 1, 1,   /* GT */
                                                                        1, 1, 1, 1, 1,   /* CA */
                                                                           1, 1, 1, 1,   /* GA */
                                                                              1, 1, 1,   /* TA */
                                                                                 1, 1,   /* GC */
                                                                                    1};  /* TC */


/*                                      AA   CC   GG   TT   AC   AG   AT   CG   CT   GT */
static const double gt_freqs_equal[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };

#define ONE_16 1./16
static const double gt16_freqs_equal[16] = {ONE_16, ONE_16, ONE_16, ONE_16, ONE_16, ONE_16,
                                            ONE_16, ONE_16, ONE_16, ONE_16, ONE_16, ONE_16,
                                            ONE_16, ONE_16, ONE_16, ONE_16};

/*                                 A  C  G  T              */
//static int gt_sym_freq_equal[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//static int gt_sym_freq_free[]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

/*                                  AA CC GG TT AC AG AT CG CT GT          */
static int gt_sym_rate_free_sm[] = {    0, 0, 0, 1, 2, 3, 0, 0, 0,    /* AA */
                                           0, 0, 4, 0, 0, 5, 6, 0,    /* CC */
                                              0, 0, 7, 0, 8, 0, 9,    /* GG */
                                                 0, 0,10, 0,11,12,    /* TT */
                                                   13,14,15,16, 0,    /* AC */
                                                      17,18, 0,19,    /* AG */
                                                          0,20,21,    /* AT */
                                                            22,23,    /* CG */
                                                               24 };  /* CT */

/* A-C: 1, A-G: 2, A-T: 3, C-G: 4, C-T: 5, G-T: 6, others: 0 */
/*                               AA CC GG TT AC AG AT CG CT GT          */
static int gt_sym_rate_dna4[] =  {   0, 0, 0, 1, 2, 3, 0, 0, 0,    /* AA */
                                        0, 0, 1, 0, 0, 4, 5, 0,    /* CC */
                                           0, 0, 2, 0, 4, 0, 6,    /* GG */
                                              0, 0, 3, 0, 5, 6,    /* TT */
                                                 4, 5, 2, 3, 0,    /* AC */
                                                    6, 1, 0, 3,    /* AG */
                                                       0, 1, 2,    /* AT */
                                                          6, 5,    /* CG */
                                                             4 };  /* CT */

/* A-C = A-T = C-G = G-T: 1, A-G = C-T: 2, others: 0 */
/*                               AA CC GG TT AC AG AT CG CT GT          */
static int gt_sym_rate_hky4[] =  {   0, 0, 0, 1, 2, 1, 0, 0, 0,    /* AA */
                                        0, 0, 1, 0, 0, 1, 2, 0,    /* CC */
                                           0, 0, 2, 0, 1, 0, 1,    /* GG */
                                              0, 0, 1, 0, 2, 1,    /* TT */
                                                 1, 2, 2, 1, 0,    /* AC */
                                                    1, 1, 0, 1,    /* AG */
                                                       0, 1, 2,    /* AT */
                                                          1, 2,    /* CG */
                                                             1 };  /* CT */


/* A-C: 1, A-G: 2, A-T: 3, C-G: 4, C-T: 5, G-T: 6, others: 0 */
/*                                    AA CC GG TT AC AG AT CG CT GT CA GA TA GC TC TG          */
static int gt16_sym_rate_dna4[] =  {      0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 3, 0, 0, 0,   /* AA */
                                             0, 0, 1, 0, 0, 4, 5, 0, 1, 0, 0, 4, 5, 0,   /* CC */
                                                0, 0, 2, 0, 4, 0, 6, 0, 2, 0, 4, 0, 6,   /* GG */
                                                   0, 0, 3, 0, 5, 6, 0, 0, 3, 0, 5, 6,   /* TT */
                                                      4, 5, 2, 3, 0, 0, 0, 0, 2, 3, 0,   /* AC */
                                                         6, 1, 0, 3, 0, 0, 0, 0, 0, 3,   /* AG */
                                                            0, 1, 2, 0, 0, 0, 0, 0, 0,   /* AT */
                                                               6, 5, 2, 0, 0, 0, 0, 5,   /* CG */
                                                                  4, 3, 0, 0, 0, 0, 0,   /* CT */
                                                                     0, 3, 0, 5, 0, 0,   /* GT */
                                                                        4, 5, 0, 0, 0,   /* CA */
                                                                           6, 1, 0, 0,   /* GA */
                                                                              0, 1, 2,   /* TA */
                                                                                 6, 0,   /* GC */
                                                                                    4    /* TC */
};


static const pllmod_subst_model_t gt_model_list[] =
{
/*  name    states  model rates         model freqs   rate symmetries   freq. sym.           */
  {"GT10",       10, NULL,               NULL,              gt_sym_rate_dna4,    NULL, 0 },
  {"GT10JC-SM",  10, gt_rates_equal_sm,  gt_freqs_equal,    NULL,                NULL, 0 },
  {"GT10JC",     10, gt_rates_equal,     gt_freqs_equal,    NULL,                NULL, 0 },
  {"GT10GTR-SM", 10, NULL,               NULL,              gt_sym_rate_free_sm, NULL, 0 },
  {"GT10HKY",    10, NULL,               NULL,              gt_sym_rate_hky4,    NULL, 0 },
  {"GT10GTR",    10, NULL,               NULL,              NULL,                NULL, 0 },
  {"GT16",       16, NULL,               NULL,              gt16_sym_rate_dna4,  NULL, 0 },
  {"GT16JC",     16, gt16_rates_equal,   gt16_freqs_equal,  NULL,                NULL, 0 },
  {"GT16GTR",    16, NULL,               NULL,              NULL,                NULL, 0 }
};

const int GT_MODELS_COUNT = sizeof(gt_model_list) / sizeof(pllmod_subst_model_t);

//static const pllmod_subst_model_alias_t gt_model_aliases[] =
//{
//  {"GTJC",     "GT10JC"},
//  {"GTJC-SM",  "GT10JC-SM"},
//  {"GTGTR4",   "GT10"},
//  {"GTGTR",    "GT10GTR"},
//  {"GTGTR-SM", "GT10GTR-SM"},
//  {"GTHKY4",   "GT10HKY"},
//  {"GPGTR4",   "GT16"}
//};

#endif /* genotype_models_hpp */

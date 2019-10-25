
extern "C"
{
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
}

#ifndef macros_h
#define macros_h

/* Macros */
#define PROGRAM_NAME    "TumorCoal"
#define VERSION_NUMBER    "2.0"
#define NO          0
#define YES         1
#define INCREMENT_NODES   500

#define	infinity                999999
#define	MATERNAL                0
#define	PATERNAL                1
#define BINARY                  0
#define	DNA                     1
#define    MAX_NAME             120
#define    MAX_LINE                3500
#define    MAX_GENOME              1000000

#define ISMhap					0
#define Mk						1
#define finiteDNA				2

#define A						0
#define C						1
#define G						2
#define T						3
#define N                       4
#define R                       5
#define ADO						9
#define DELETION				7
#define M                       6
#define S                       8
#define Y                       10
#define K                       11


//#define genotypes
#define AA                      0
#define CC                      1
#define GG                      2
#define TT                      3
#define AC                      4
#define CA                      4
#define AG                      5
#define GA                      5
#define AT                      6
#define TA                      6
#define A_                      7
#define _A                      7
#define CG                      8
#define GC                      8
#define CT                      9
#define TC                      9
#define C_                      10
#define _C                      10
#define GT                      11
#define TG                      11
#define G_                      12
#define _G                      12
#define T_                      13
#define _T                      13
#define __                      14
//A/C => M
//A/G => R
//A/T => W
//A/_ => a
//
//C/A => M
//C/C => C
//C/G => S
//C/T => Y
//C/_ => c
//
//G/A => R
//G/C => S
//G/G => G
//G/T => K
//G/_ => g
//
//T/A => W
//T/C => Y
//T/G => K
//T/T => T
//T/_ => t
//
//_/A => a
//_/C => c
//_/G => g
//_/T => t
//_/_ => -


//#define M                       6 //A/C or C/A
//#define W                       8 //A/T or T/A
//#define S                        //C/G or G/C
//#define Y                       14 //C/T or T/C
//#define K                       15 //G/T or T/G
//#define a                       16 //A/_ or _/A
//#define c                       17 //C/_ or _/C
//#define g                       18 //G/_ or _/G
//#define t                       19 //T/_ or _/T

#define CG_AT					0	// C>A or G>T
#define CG_GC					1	// C>G or G>C
#define CG_TA					2	// C>T or G>A
#define TA_AT					3	// T>A or A>T
#define TA_CG					4	// T>C or A>G
#define TA_GC					5	// T>G or A>C

//#define M                       10 //A/C or C/A
////#define R                       11 //A/G or G/A
//#define W                       12 //A/T or T/A
//#define S                       13 //C/G or G/C
//#define Y                       14 //C/T or T/C
//#define K                       15 //G/T or T/G
//#define a                       16 //A/_ or _/A
//#define c                       17 //C/_ or _/C
//#define g                       18 //G/_ or _/G
//#define t                       19 //T/_ or _/T
//#define _                       20 //_/_


#define ISMhap					0
#define Mk2						1
#define finiteDNA				2

#define NUM_SIGNATURES			30

#define trinuc(i,j,k)   (i*16 + j*4 + k)
#define trimut(i,j,k)   (i*24 + j*4 + k)


#define CHECK_MUT_DISTRIBUTION
#undef CHECK_MUT_DISTRIBUTION

#define MYDEBUG
#undef MYDEBUG

#define LOAD_INTERNAL_SIGNATURES

#define PRINT_TRIMUTATIONS
#undef PRINT_TRIMUTATIONS

#define PRINT_TRINUC_GENOME
#undef PRINT_TRINUC_GENOME

#define PRINT_TRIMUTCOUNTER
#undef PRINT_TRIMUTCOUNTER


#define GT_MODEL "GTGTR4"
#define JC_MODEL "GTJC"
#define RATE_CATS 1
#define BRLEN_MIN 1e-6
#define BRLEN_MAX 1e+2

/* ascertainment bias correction */
#define PLL_ATTRIB_AB_LEWIS        (1 << 5)
#define PLL_ATTRIB_AB_FELSENSTEIN  (2 << 5)
#define PLL_ATTRIB_AB_STAMATAKIS   (3 << 5)
#define PLL_ATTRIB_AB_MASK         (7 << 5)
#define PLL_ATTRIB_AB_FLAG         (1 << 8)

#define PLL_ATTRIB_RATE_SCALERS    (1 << 9)

#define HOMO(state)   (state<4)
#define HETERO(state) (state>3)

#endif /* macros_h */
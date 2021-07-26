#ifndef definitions_hpp
#define definitions_hpp

/* Macros */
#define PROGRAM_NAME    "TumorCoal"
#define VERSION_NUMBER    "2.0"
#define NO          0
#define YES         1
#define INCREMENT_NODES   500

#define	infinity2                999999
#define	MATERNAL                0
#define	PATERNAL                1
#define BINARY                  0
#define	DNA                     1
#define    MAX_NAME             200
#define    MAX_LINE                3500
#define    MAX_GENOME              1000000

#define ISMhap					0
#define Mk						1
#define finiteDNA				2
namespace Definitions{
    static const int  A = 0;
    static const int  C = 1;
    static const int  G = 2;
    static const int  T= 3;
    static const int  N = 4;
    static const int  R = 5;
    static const int  ADO = 9;
    static const int  DELETION	= 7;
    static const int  M = 6;
    static const int S = 8;
    static const int Y = 10;
    static const int K = 11;

    static const int AA = 0;
    static const int CC = 1;
    static const int GG = 2;
    static const int TT = 3;
    static const int AC = 4;
    static const int CA = 4;
    static const int AG = 5;
    static const int GA = 5;
    static const int AT = 6;
    static const int TA = 6;
    static const int A_ = 7;
    static const int _A = 7;
    static const int CG = 8;
    static const int GC = 8;
    static const int CT = 9;
    static const int TC = 9;
    static const int C_ = 10;
    static const int _C = 10;
    static const int GT = 11;
    static const int TG = 11;
    static const int G_ = 12;
    static const int _G = 12;
    static const int T_ = 13;
    static const int _T = 13;
    static const int __ = 14;


};



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
//#define STATES 4


#endif /* definitions_h */

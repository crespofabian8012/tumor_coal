/*################################################################################
 ##
 ##   Copyright (C) 2018-2020 Fausto Fabian Crespo Fernandez
 ##
 ##   This file is part of the tumor_coal C++ library.
 ##
 ##   Licensed under the Apache License, Version 2.0 (the "License");
 ##   you may not use this file except in compliance with the License.
 ##   You may obtain a copy of the License at
 ##
 ##       http://www.apache.org/licenses/LICENSE-2.0
 ##
 ##   Unless required by applicable law or agreed to in writing, software
 ##   distributed under the License is distributed on an "AS IS" BASIS,
 ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 ##   See the License for the specific language governing permissions and
 ##   limitations under the License.
 ##
 ################################################################################*/

/*
 * utils functions
 */
#include "utils.hpp"
#include <unistd.h>

extern "C"
{
#include <unistd.h>
#include "kseq.h"
#include <gsl/gsl_errno.h>
#include <gsl_complex.h>
#include <gsl_fft_complex.h>
#include <gsl_fft_real.h>
#include <gsl_fft_halfcomplex.h>
}


#include <algorithm>

//#include <fftw3.h>
#include "output_functions.hpp"
//#include "autocovariance.hpp"


KSEQ_INIT(int, read);

/************************ ReadParametersFromFastaFile ***********************/
/*  ReadParametersFromFastaFile */
void Utils::ReadParametersFromFastaFile(char *fileName, int &numCells, int &TotalNumSequences, int &numSites){
    //read fasta
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    
    
    int max_length=0.0;
    int numberSeq;
    
    
    if ((fastaFile = freopen(fileName, "r", stdin)) != NULL)
    {
        seq = kseq_init(fileno(fastaFile));
        numberSeq=0;
        while ((l1 = kseq_read(seq)) >= 0 )
        {
            
            numberSeq=numberSeq+1;
            if ( l1 > max_length)
            {
                max_length=l1;
            }
        }
        
        numSites=max_length;
        if (numberSeq >=1){
            // *numCells =numberSeq;//not counting the healthy cell
            numCells=numberSeq;
            TotalNumSequences=numberSeq;
        }
        else{
            //   *numCells =0;
           numCells=numberSeq;
           TotalNumSequences=numberSeq;
        }
        kseq_destroy(seq);
        fclose(fastaFile);
    }
    else{
        
        fprintf (stderr, "\nERROR: Can't read fasta file.");
        Output::PrintUsage();
        
    }
}
/************************ ReadFastaFile ***********************/
/*  ReadFastaFile */
void Utils::ReadFastaFile(char *fileName, std::vector<std::vector<int> > &ObservedData,  char **ObservedCellNames, ProgramOptions &programOptions){
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    int current=0;
    char *currentSeq;
    char *currentQual;
    int  seqlength=0;
    char *t;
   
    std::vector<int> v;
    if ((fastaFile = freopen(fileName, "r", stdin)) != NULL){
        seq = kseq_init(fileno(fastaFile));
        while ((l1 = kseq_read(seq)) >= 0 ) {
            ObservedCellNames[current] = (char*) malloc(MAX_NAME);
            if (ObservedCellNames[current] != NULL)
            {
                strcpy(ObservedCellNames[current], seq->name.s);
                
            }
            
            seqlength=0;
            currentSeq=seq->seq.s;
            // for ( t= currentSeq; *t != '\0'; t++) {
            std::vector<int> v;
            for ( int index= 0; index < l1; index++) {
                t= currentSeq+index;
                if (programOptions.doUseGenotypes == NO) // use sequences
                    v.push_back( WhichNucChar(*t));
                //ObservedData[current][index]= WhichNucChar(*t);
                else // use genotypypes
                    v.push_back(WhichGenotypeChar(*t));
                seqlength++;
            }
            ObservedData.push_back(v);
            current++;
            if (seq->qual.l){
                // printf("qual: %s\n", seq->qual.s);
                currentQual=seq->qual.s;
            }
        }
        //printf("return value: %d\n", l1);
        kseq_destroy(seq);
        //  gzclose(fastaFile);
        fclose(fastaFile);
    }
}


/********************* WhichIUPAC ************************/
/* Returns the IUPAC representation of the genotype */
/*
 UPAC nucleotide code    Base
 A    Adenine
 C    Cytosine
 G    Guanine
 T (or U)    Thymine (or Uracil)
 R    A or G
 Y    C or T
 S    G or C
 W    A or T
 K    G or T
 M    A or C
 B    C or G or T
 D    A or G or T
 H    A or C or T
 V    A or C or G
 N    unknown state
 . or -    gap
 
 This is what we do:
 
 A/A => A
 A/C => M
 A/G => R
 A/T => W
 A/_ => a
 
 C/A => M
 C/C => C
 C/G => S
 C/T => Y
 C/_ => c
 
 G/A => R
 G/C => S
 G/G => G
 G/T => K
 G/_ => g
 
 T/A => W
 T/C => Y
 T/G => K
 T/T => T
 T/_ => t
 
 _/A => a
 _/C => c
 _/G => g
 _/T => t
 _/_ => -
 
 */

char Utils::WhichIUPAC (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //AA
            return ('A');
        else if (allele2 == 1)    //AC
            return ('M');
        else if (allele2 == 2)    //AG
            return ('R');
        else if (allele2 == 3)    //AT
            return ('W');
        else if (allele2 == Definitions::ADO)    //A?
            return ('a');
        else if (allele2 == Definitions::DELETION)    //A–
            return ('a');
        else
            return ('N');
    }
    else if (allele1 == 1)
    {
        if (allele2 == 0)        //CA
            return ('M');
        else if (allele2 == 1)    //CC
            return ('C');
        else if (allele2 == 2)    //CG
            return ('S');
        else if (allele2 == 3)    //CT
            return ('Y');
        else if (allele2 == Definitions::ADO)    //C?
            return ('c');
        else if (allele2 == Definitions::DELETION)    //C–
            return ('c');
        else
            return ('N');
    }
    else if (allele1 == 2)
    {
        if (allele2 == 0)        //GA
            return ('R');
        else if (allele2 == 1)    //GC
            return ('S');
        else if (allele2 == 2)    //GG
            return ('G');
        else if (allele2 == 3)    //GT
            return ('K');
        else if (allele2 == Definitions::ADO)    //G?
            return ('g');
        else if (allele2 == Definitions::DELETION)    //G–
            return ('g');
        else
            return ('N');
    }
    else if (allele1 == 3)
    {
        if (allele2 == 0)        //TA
            return ('W');
        else if (allele2 == 1)    //TC
            return ('Y');
        else if (allele2 == 2)    //TG
            return ('K');
        else if (allele2 == 3)    //TT
            return ('T');
        else if (allele2 == Definitions::ADO)    //T?
            return ('t');
        else if (allele2 == Definitions::DELETION)    //T–
            return ('t');
        else
            return ('N');
    }
    else if (allele1 == Definitions::ADO)
    {
        if (allele2 == 0)        //?A
            return ('a');
        else if (allele2 == 1)    //?C
            return ('c');
        else if (allele2 == 2)    //?G
            return ('g');
        else if (allele2 == 3)    //?T
            return ('t');
        else if (allele2 == Definitions::ADO)    //??
            return ('-');
        else if (allele2 == Definitions::DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == Definitions::DELETION)
    {
        if (allele2 == 0)        //-A
            return ('a');
        else if (allele2 == 1)    //-C
            return ('c');
        else if (allele2 == 2)    //-G
            return ('g');
        else if (allele2 == 3)    //-T
            return ('t');
        else if (allele2 == Definitions::ADO)    //-?
            return ('-');
        else if (allele2 == Definitions::DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}

/********************* WhichMut ************************/
/* Returns character representation for binary data */

char Utils::WhichMut (int state)
{
    if (state == 0)
        return ('0');
    else if (state == 1)
        return ('1');
    else if (state == Definitions::ADO)
        return ('?');
    else if (state == Definitions::DELETION)
        return ('-');
    else
        return ('N');
}
/********************* WhichConsensusBinary ************************/
/* Returns a consensus representation of the binary genotype */
/*
 0/0 => 0
 0/1 => 1
 1/0 => 1
 1/1 => 2
 
 0/_ => 0
 _/0 => 0
 
 1/_ => 2
 _/1 => 2
 
 _/_ => -
 */

char Utils::WhichConsensusBinary (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //00
            return ('0');
        else if (allele2 == 1)    //01
            return ('1');
        else if (allele2 == Definitions::ADO)    //0?
            return ('0');
        else if (allele2 == Definitions::DELETION)    //0-
            return ('0');
        else
            return ('N');
    }
    else if (allele1 == 1)
    {
        if (allele2 == 0)        //10
            return ('1');
        else if (allele2 == 1)    //11
            return ('2');
        else if (allele2 == Definitions::ADO)    //1?
            return ('2');
        else if (allele2 == Definitions::DELETION)    //0-
            return ('2');
        else
            return ('N');
    }
    else if (allele1 == Definitions::ADO)
    {
        if (allele2 == 0)        //?0
            return ('0');
        else if (allele2 == 1)    //?1
            return ('2');
        else if (allele2 == Definitions::ADO)    //??
            return ('-');
        else if (allele2 == Definitions::DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == Definitions::DELETION)
    {
        if (allele2 == 0)        //-0
            return ('0');
        else if (allele2 == 1)    //-1
            return ('2');
        else if (allele2 == Definitions::ADO)    //-?
            return ('-');
        else if (allele2 == Definitions::DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}

/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char Utils::WhichNuc (int nucleotide)
{
    if (nucleotide == Definitions::A)
        return ('A');
    else if (nucleotide == Definitions::C)
        return ('C');
    else if (nucleotide == Definitions::G)
        return ('G');
    else if (nucleotide == Definitions::T)
        return ('T');
    else if (nucleotide == Definitions::ADO)
        return ('?');
    else if (nucleotide == Definitions::DELETION)
        return ('-');
    else
        return ('N');
}

/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int Utils::CheckMatrixSymmetry(double matrix[4][4])
{
    int i,j;
    
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            if(matrix[i][j] != matrix[j][i])
                return NO;
    return YES;
}

/********************* WhichNucChar ************************/
/* Returns integer representation for character nucleotudes */

int Utils::WhichNucChar (char nucleotide)
{
    if (nucleotide == 'A')
        return (Definitions::A);
    else if (nucleotide == 'C')
        return (Definitions::C);
    else if (nucleotide == 'G')
        return (Definitions::G);
    else if (nucleotide == 'T')
        return (Definitions::T);
    else if (nucleotide == '?')
        return (Definitions::ADO);
    else if (nucleotide == '-')
        return (Definitions::DELETION);
    else if (nucleotide == 'N')
        return (Definitions::N);
//    else if (nucleotide == 'R')
//        return (Definitions::R);
    else
    {
        fprintf (stderr, "\nERROR in WhichNucChar: nucleotide = %c\n",  nucleotide);
        exit(-1);
    }
}

/********************* WhichGenotypeChar ************************/
/* Returns integer representation for character nucleotudes */

int Utils::WhichGenotypeChar (char nucleotide)
{

    if (nucleotide == 'A')
        return (Definitions::AA);
    else if (nucleotide == 'C')
        return (Definitions::CC);
    else if (nucleotide == 'G')
        return (Definitions::GG);
    else if (nucleotide == 'T')
        return (Definitions::TT);
    else if (nucleotide == '?')
        return (Definitions::__);
    else if (nucleotide == '-')
        return (Definitions::__);
    else if (nucleotide == 'N')
        return (Definitions::N);
    else if (nucleotide == 'R')
        return (Definitions::AG);
    else if (nucleotide == 'M')
        return (Definitions::AC);
    else if (nucleotide == 'W')
        return (Definitions::AT);
    else if (nucleotide == 'S')
        return (Definitions::CG);
    else if (nucleotide == 'Y')
        return (Definitions::CT);
    else if (nucleotide == 'K')
        return (Definitions::GT);
    else if (nucleotide == 'a')
        return (Definitions::A_);
    else if (nucleotide == 'c')
        return (Definitions::C_);
    else if (nucleotide == 'g')
        return (Definitions::G_);
    else if (nucleotide == 't')
        return (Definitions::T_);
    else
    {
        fprintf (stderr, "\nERROR in WhichGenotypeChar: nucleotide = %c\n",  nucleotide);
        exit(-1);
    }
}

/********************* parameterMultiplierMCMCmove ************************/


double Utils::parameterMultiplierMCMCmove (double lengthInterval)
{
    double result = 0.0;
    
    result = (lengthInterval)/2.0 + sqrt(1 + (lengthInterval * lengthInterval)* 0.25);
    return result;
}
/********************* SequenceToIntegers ************************/

std::vector<int> Utils::SequenceToIntegers(char ** sequence, int length)
{
    std::vector<int> result(length);
    for(unsigned int i=0; i< length; i++)
    {
        result[i]=WhichNucChar(*sequence[i]);
        
    }
    return(result);
}
/********************* GenotypesToIntegers ************************/

std::vector<int> Utils::GenotypesToIntegers(char ** sequence, int length)
{
    std::vector<int> result(length);
    char * c;
    for(unsigned int i=0; i< length; i++)
    {   c=sequence[i];
        result[i]=WhichGenotypeChar (*c);
        
    }
    return(result);
}
/********************* cb_serialize ************************/

char * Utils::cb_serialize(const pll_rnode_t * node){
    char * newick=NULL;
    if (!node) return NULL;
    if (node->left ==NULL && node->right==NULL)
       asprintf(&newick, "%s:%.15f", node->label , node->length);
    else
       asprintf(&newick, "%s:%.15f", "" , node->length);
    return(newick);
}
/********************* appendCharToCharArray ************************/
char* Utils::appendCharToCharArray(char* array, char a)
{
    size_t len = strlen(array);
    
    char* ret = new char[len+2];
    
    strcpy(ret, array);
    ret[len] = a;
    ret[len+1] = '\0';
    
    return ret;
}
void Utils::init_to_empty_str(char str[MAX_NAME])
{
    for (size_t i = 0; i < MAX_NAME; i++) {
        str[i] = 0;
    }
}
std::vector<pll_rnode_t*>  Utils::filterHealthyTip(pll_rtree_t * trueTree, pll_rnode_t** doublePointer, size_t count, std::string &healthyTipLabel )
{
    std::vector<pll_rnode_t*> result;
    pll_rnode_t* current;
    for(size_t j=0; j <count;++j)
    {
        
        current =trueTree->nodes[j];
        
        if  (current==NULL)
          fprintf (stderr, "\nERROR: A node is null \n");
        
        if  (current!=NULL && current->label != 0 && current->label  ){
            
            if (std::string(current->label).compare(healthyTipLabel)!=0)
                result.push_back(doublePointer[j]);
            
        }
        else
           result.push_back(doublePointer[j]);
     }
return result;
}
//vector<long double >    Utils::potentialScaleReductionLongDouble(int inner_size,int n, int m,  vector<vector<long double>> &means, vector<vector<long double>> &variances )
long double     Utils::potentialScaleReductionLongDouble(int inner_size,int n, int m,  std::vector<long double>& means, std::vector<long double> &variances )
{
    std::vector<long double> elements;
    std::vector<long double> varianceElements;
    std::vector<long double> results;

    long double B;
    long double W;
    long double potentialScaleReduction ;
    
    B = n * Utils::variance(means) ;
    W = (1.0 /m)* Utils::mean(variances) ;
    potentialScaleReduction= (1.0-1.0/n)* W+ (1.0/n)*B;
    potentialScaleReduction = sqrt(potentialScaleReduction / W);

    return potentialScaleReduction;
    
//
//    for(int j=0; j < inner_size;j++)
//    {
//
//        std::transform( means.begin(), means.end(), elements.begin(), [&](vector<long double> &vec){
//            return vec.at(j);
//        });
//        std::transform( variances.begin(), variances.end(), varianceElements.begin(), [&](vector<long double> &vec){
//            return vec.at(j);
//        });
//
//        meanOfMeans = Utils::mean( elements);
//        B = n * Utils::variance(elements) ;
//        W = (1.0 /m)*Utils::mean(varianceElements) ;
//        potentialScaleReduction= (1.0-1.0/n)* W+ (1.0/n)*B;
//        results.push_back(potentialScaleReduction);
//    }
//    return results;
    }
std::vector<long double >    Utils::potentialScaleReductionArray(int numClones,int n, int m,  std::vector<std::vector<long double>> &means, std::vector<std::vector<long double>> &variances ){
    std::vector<long double> elements;
    std::vector<long double> varianceElements;
    std::vector<long double> results;
  
    long double potentialScaleReduction ;
    std::vector<long double > temp;

    for(int j=0; j < numClones;j++)
        {
            elements.clear();
            varianceElements.clear();

            for(int i=0; i < means.size();i++){
                temp = means.at(i);
                elements.push_back(temp.at(j));
                temp=variances.at(i);
                varianceElements.push_back(temp.at(j));
            }
            
            potentialScaleReduction=  Utils::potentialScaleReductionLongDouble(numClones, n, m,  elements, varianceElements);
            results.push_back(potentialScaleReduction);
        }
    

    return results;
}
bool Utils::checkAllElementsGreaterThanZero(std::vector<long double> &values, int size){
    
    bool result=true;
    for(int j=0; j < size;j++){
        
        if (values[j]<=0){
            result =false;
            fprintf (stderr, "\nERROR: Coal time %d: %Lf is not positive \n", j,values[j] );
            break;
            
        }
        
    }
    
    return result;
}
 // Following FFT: https://cp-algorithms.com/algebra/fft.html
void Utils::fft (std::vector<linalgebra::complex> & a, bool invert)
{
  int n = (int) a.size();
 const double PI = 3.141592653589793;
    
  for (int i=1, j=0; i<n; ++i) {
    int bit = n >> 1;
    for (; j>=bit; bit>>=1)
      j -= bit;
    j += bit;
    if (i < j)
      SWAP(a[i], a[j]);
  }
 
  for (int len=2; len<=n; len<<=1) {
    double ang = 2*PI/len * (invert ? -1 : 1);
    linalgebra::complex wlen(cos(ang), sin(ang));
    for (int i=0; i<n; i+=len) {
      linalgebra::complex w(1.0, 0.0);
      for (int j=0; j<len/2; ++j) {
          linalgebra::complex u = a[i+j];
        linalgebra::complex v = a[i+j+len/2] * w;
        a[i+j] = u + v;
        a[i+j+len/2] = u -v;
        w *= wlen;
      }
    }
  }
  if (invert)
    for (int i=0; i<n; ++i)
      a[i] /= n;
}
int Utils::next2Power(int n){
    int i = 1;
    while(i < n)
        i = i << 1;
    return i;
}
inline int Utils::pow2i(int x) { return ((x < 0) ? 0 : (1 << x)); };

long double Utils::autoCorrelation2(std::vector<long double> & a){
    
    int n = next2Power(a.size());
    
    double data[2*n];

    for (unsigned int i = 0; i < n; i++)
      {
         REAL(data,i) = 0.0; IMAG(data,i) = 0.0;
      }

    REAL(data,0) = 1.0;

    for (unsigned int i = 1; i <= 10; i++)
      {
         REAL(data,i) = REAL(data,n-i) = 1.0;
      }
    gsl_fft_complex_radix2_forward (data, 1, n);
    long double result =0.0;
    
    return result;
}
long double Utils::fftReal(double  *data, int n){
    
   
    gsl_fft_real_wavetable * real;
    gsl_fft_halfcomplex_wavetable * hc;
    gsl_fft_real_workspace * work;
    
     work = gsl_fft_real_workspace_alloc (n);
     real = gsl_fft_real_wavetable_alloc (n);
    
    
    gsl_fft_real_transform(data, 1, n, real, work);

    
    gsl_fft_real_wavetable_free (real);
    
    hc = gsl_fft_halfcomplex_wavetable_alloc (n);

    gsl_fft_halfcomplex_inverse(data, 1, n, hc, work);
    
    gsl_fft_halfcomplex_wavetable_free (hc);
    
    
     gsl_fft_real_workspace_free (work);
    
    
    long double result =0.0;
    
    return result;
}
void Utils::four1(long double *data, const int n, const int isign)
{
    int nn,mmax,m,j,istep,i;
    long double  wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

    if (n<2 || n&(n-1)) throw("n must be power of 2 in four1");
    nn = n << 1;
    j = 1;
    for (i=1;i<nn;i+=2) {
        if (j > i) {
            SWAP(data[j-1],data[i-1]);
            SWAP(data[j],data[i]);
        }
        m=n;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (nn > mmax) {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=nn;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j-1]-wi*data[j];
                tempi=wr*data[j]+wi*data[j-1];
                data[j-1]=data[i-1]-tempr;
                data[j]=data[i]-tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

void Utils::four1(std::vector<long double> &data, const int isign) {
    four1(&data[0],data.size()/2,isign);
}

//void Utils::four1(std::vector<complex> &data, const int isign) {
//    four1((long double*)(&data[0]),data.size(),isign);
//}
void Utils::realft(std::vector<long double> &data, const int isign)
{
    int i,i1,i2,i3,i4;
    long double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;

    int n=data.size();
    long double PI= 3.141592653589793;
    theta = PI/double(n>>1);
   
    if (isign == 1) {
        c2 = -0.5;
        four1(data,1);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    for (i=1;i<(n>>2);i++) {
        i2=1+(i1=i+i);
        i4=1+(i3=n-i1);
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r= -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4]= -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[0] = (h1r=data[0])+data[1];
        data[1] = h1r-data[1];
    } else {
        data[0]=c1*((h1r=data[0])+data[1]);
        data[1]=c1*(h1r-data[1]);
        four1(data,-1);
    }
}
void Utils::correl(const std::vector<long double> &data1, const std::vector<long double> &data2, std::vector<long double> &ans)
{
    int no2,i;
    long double  tmp;

    int n=data1.size();
    std::vector<long double> temp(n);
    for (i=0;i<n;i++) {
        ans[i]=data1[i];
        temp[i]=data2[i];
    }
    realft(ans,1);
    realft(temp,1);
    no2=n>>1;
    for (i=2;i<n;i+=2) {
        tmp=ans[i];
        ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/no2;
        ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/no2;
    }
    ans[0]=ans[0]*temp[0]/no2;
    ans[1]=ans[1]*temp[1]/no2;
    realft(ans,-1);
}
void Utils::correl2(const std::vector<long double> &data,long double mean, std::vector<long double> &ans)
{
    int no2,i;
    long double  tmp;

    int N=data.size();
    //size_t M = Utils::next2Power( N);
    std::vector<long double> temp(N);
    for (i=0;i<N;i++) {
        ans[i]=data[i]-mean;
        temp[i]=data[i]-mean;
    }
    realft(ans,1);
    realft(temp,1);
    
    no2 = N>>1;
    for (i=2;i<N;i+=2) {
        tmp=ans[i];
        ans[i]=(ans[i]*temp[i]+ans[i+1]*temp[i+1])/no2;
        ans[i+1]=(ans[i+1]*temp[i]-tmp*temp[i+1])/no2;
    }
    ans[0]=ans[0]*ans[0]/no2;
    ans[1]=ans[1]*ans[1]/no2;
    realft(ans,-1);
}
long double   Utils::tune(long double scale, const long double acceptanceRate){
    
    assert(scale >0);
    
    long double result= scale;
    if (acceptanceRate < 0.001){
        // reduce by 90 percent
        result *= 0.1;
    }
    else if (acceptanceRate < 0.05){
       // reduce by 50 percent
            result *= 0.5;
        }
    else if (acceptanceRate < 0.3){
        // reduce by ten percent
        result *= 0.9;
        }
    else if (acceptanceRate > 0.95){
        // increase by factor of ten
        result *= 10.0;
        if (result > 2.0)
            result= 2.0;
    }
    else if (acceptanceRate > 0.75){
        // increase by double
        result *= 2.0;
        if (result > 2.0)
          result= 2.0;
        
    }
    else if (acceptanceRate > 0.5){
        // increase by ten percent
        result *= 1.1;
        if (result > 2.0)
          result= 2.0;
        
    }
    return result;
}
double Utils::normalize(const std::vector<double> &log_weights, std::vector<double> &weights)
{
    // compute the lognorm
    double log_norm = log_add(log_weights);
    normalize(log_weights, weights, log_norm);
    return log_norm;
}

void Utils::normalize(const std::vector<double> &log_weights, std::vector<double> &weights, double log_norm)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < log_weights.size(); i++)
    {
        weights[i] = exp(log_weights[i] - log_norm);
        sum += weights[i];
    }
    if (abs(sum - 1.0) > 1e-4) {
        std::cout << "Error in normalization. Check that log_weights and log_norm are correctly calculated." << std::endl;
        std::cout << ceil(sum*100000)/100000.0 << " != 1.0" << std::endl;
        std::cout << "log_norm: " << log_norm << std::endl;
        exit(-1);
    }
}
double Utils::log_add(double x, double y)
{
    // make x the max
    if (y > x) {
        double temp = x;
        x = y;
        y = temp;
    }
    // now x is bigger
    if (x == DOUBLE_NEG_INF) {
        return x;
    }
    double negDiff = y - x;
    if (negDiff < -20) {
        return x;
    }
    return x + log(1.0 + exp(negDiff));
}
double Utils::log_add(std::vector<double> x)
{
    double max = DOUBLE_NEG_INF;
    double maxIndex = 0;
    for (unsigned int i = 0; i < x.size(); i++)
    {
        if (x[i] > max) {
            max = x[i];
            maxIndex = i;
        }
    }
    if (max == DOUBLE_NEG_INF) return DOUBLE_NEG_INF;
    // compute the negative difference
    double threshold = max - 20;
    double sumNegativeDifferences = 0.0;
    for (unsigned int i = 0; i < x.size(); i++) {
        if (i != maxIndex && x[i] > threshold) {
            sumNegativeDifferences += exp(x[i] - max);
        }
    }
    if (sumNegativeDifferences > 0.0) {
        return max + log(1.0 + sumNegativeDifferences);
    } else {
        return max;
    }
    
}
std::vector<std::pair<int, int>> Utils::allPairCombinations(int vectorSize){
    
    assert(vectorSize >1);
    std::vector<std::pair<int, int>> result;
    std::pair<int, int> pair;
    for(size_t i=0; i < vectorSize -1 ; i++){
        pair.first = i;
        for(size_t j = i+1; j < vectorSize;j++){
              pair.second = j;
              result.push_back(pair);
        }
        
    }

    return result;
}
std::vector<std::pair<int, int>> Utils::allPairs(int vectorSize){
    
    assert(vectorSize >1);
    std::vector<std::pair<int, int>> result;
    std::pair<int, int> pair;
    for(size_t i=0; i < vectorSize; i++){
        pair.first = i;
        for(size_t j = 0; j < vectorSize;j++){
            if (i!=j){
                pair.second = j;
                result.push_back(pair);
            }
        }
    }

    return result;
}
//void Utils::autocorrelationReal(const Eigen::MatrixBase<long double>& data, long double mean, Eigen::MatrixBase<long double>& ac, Eigen::FFT<long double>& fft)
//{
//    int no2,i;
//    long double  tmp;
//    
//    size_t N = data.size();
//    size_t M = Utils::next2Power( N);
//     
//    
//    Eigen::Matrix<long double, Eigen::Dynamic, 1> centered(2*M);
//    centered.setZero();
//    centered.head(N) =  data.array() - data.mean();
//   
//    Eigen::Matrix<std::complex<long double>, Eigen::Dynamic, 1> freqvec(2*M);
//    
//    fft.fwd(freqvec, centered);
//    
//    
//    freqvec = freqvec.cwiseAbs2();
//  
//    Eigen::Matrix<std::complex<long double>, Eigen::Dynamic, 1> ac_tmp(2*M);
//    fft.inv(ac_tmp, freqvec);
//
//    ac = ac_tmp.head(N).real().array() / (N * N * 2);
//    ac /= ac(0);
//}


//void Utils::multiply (vector<int> & a, vector<int> & b, vector<int> & res)
//{
//  
//  vector<complex> fa;
//  vector<complex> fb;
//
//complex temp = complex(0);
// for (size_t i=0; i<a.size(); ++i){
//     temp= complex(a[i]);
//     fa.push_back(temp);
//    }
//      
//    
//  for (size_t i=0; i<b.size(); ++i){
//   temp= complex(b[i]);
//   fb.push_back(temp);
//  }
//    
//  size_t n = 1;
//  int M = max(a.size(),  b.size());
//  
//  while (n < M)  n <<= 1;
//  n <<= 1;
//  fa.resize (n);
//  fb.resize (n);
// 
//  fft (fa, false);
//  fft (fb, false);
//  for (size_t i=0; i<n; ++i)
//    fa[i] *= fb[i];
//  fft (fa, true);
// 
//  res.resize (n);
//  for (size_t i=0; i<n; ++i)
//    res[i] = (int)round(fa[i].re);
//}


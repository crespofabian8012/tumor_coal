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

#include <algorithm>


#include "kseq.h"
#include "output_functions.hpp"

KSEQ_INIT(int, read);

using namespace std;
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
void Utils::ReadFastaFile(char *fileName, vector<vector<int> > &ObservedData,  char **ObservedCellNames, ProgramOptions &programOptions){
    FILE *fastaFile;
    kseq_t *seq;
    int l1;
    int current=0;
    char *currentSeq;
    char *currentQual;
    int  seqlength=0;
    char *t;
   
    vector<int> v;
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
            vector<int> v;
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
    else if (nucleotide == 'R')
        return (Definitions::R);
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
        return (AA);
    else if (nucleotide == 'C')
        return (CC);
    else if (nucleotide == 'G')
        return (GG);
    else if (nucleotide == 'T')
        return (TT);
    else if (nucleotide == '?')
        return (__);
    else if (nucleotide == '-')
        return (__);
    else if (nucleotide == 'N')
        return (Definitions::N);
    else if (nucleotide == 'R')
        return (AG);
    else if (nucleotide == 'M')
        return (AC);
    else if (nucleotide == 'W')
        return (AT);
    else if (nucleotide == 'S')
        return (CG);
    else if (nucleotide == 'Y')
        return (CT);
    else if (nucleotide == 'K')
        return (GT);
    else if (nucleotide == 'a')
        return (A_);
    else if (nucleotide == 'c')
        return (C_);
    else if (nucleotide == 'g')
        return (G_);
    else if (nucleotide == 't')
        return (T_);
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

vector<int> Utils::SequenceToIntegers(char ** sequence, int length)
{
    vector<int> result(length);
    for(unsigned int i=0; i< length; i++)
    {
        result[i]=WhichNucChar(*sequence[i]);
        
    }
    return(result);
}
/********************* GenotypesToIntegers ************************/

vector<int> Utils::GenotypesToIntegers(char ** sequence, int length)
{
    vector<int> result(length);
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
vector<pll_rnode_t*>  Utils::filterHealthyTip(pll_rnode_t** doublePointer, size_t count, string &healthyTipLabel )
{
    vector<pll_rnode_t*> result;
    pll_rnode_t* current;
    for(size_t j=0; j <count;++j)
    {
        
        current =doublePointer[j];
        
        if  (current->label != 0 && current->label  ){
            
            if (std::string(current->label).compare(healthyTipLabel)!=0)
                result.push_back(doublePointer[j]);
            
        }
        else
           result.push_back(doublePointer[j]);
     }
return result;
}
//vector<long double >    Utils::potentialScaleReductionLongDouble(int inner_size,int n, int m,  vector<vector<long double>> &means, vector<vector<long double>> &variances )
long double     Utils::potentialScaleReductionLongDouble(int inner_size,int n, int m,  vector<long double>& means, vector<long double> &variances )
{
    vector<long double> elements;
    vector<long double> varianceElements;
    vector<long double> results;

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
vector<long double >    Utils::potentialScaleReductionArray(int numClones,int n, int m,  vector<vector<long double>> &means, vector<vector<long double>> &variances ){
    vector<long double> elements;
    vector<long double> varianceElements;
    vector<long double> results;
  
    long double potentialScaleReduction ;
    vector<long double > temp;

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
bool Utils::checkAllElementsGreaterThanZero(vector<long double> &values){
    
    bool result=true;
    for(int j=0; j < values.size();j++){
        
        if (values[j]<=0){
            result =false;
            fprintf (stderr, "\nERROR: Coal time %d: %Lf is not positive \n", j,values[j] );
            break;
            
        }
        
    }
    return result;
}

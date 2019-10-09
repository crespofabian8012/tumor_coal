//
//  utils.cpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 5/10/19.
//

#include "utils.hpp"

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

char WhichIUPAC (int allele1, int allele2)
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
        else if (allele2 == ADO)    //A?
            return ('a');
        else if (allele2 == DELETION)    //A–
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
        else if (allele2 == ADO)    //C?
            return ('c');
        else if (allele2 == DELETION)    //C–
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
        else if (allele2 == ADO)    //G?
            return ('g');
        else if (allele2 == DELETION)    //G–
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
        else if (allele2 == ADO)    //T?
            return ('t');
        else if (allele2 == DELETION)    //T–
            return ('t');
        else
            return ('N');
    }
    else if (allele1 == ADO)
    {
        if (allele2 == 0)        //?A
            return ('a');
        else if (allele2 == 1)    //?C
            return ('c');
        else if (allele2 == 2)    //?G
            return ('g');
        else if (allele2 == 3)    //?T
            return ('t');
        else if (allele2 == ADO)    //??
            return ('-');
        else if (allele2 == DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == DELETION)
    {
        if (allele2 == 0)        //-A
            return ('a');
        else if (allele2 == 1)    //-C
            return ('c');
        else if (allele2 == 2)    //-G
            return ('g');
        else if (allele2 == 3)    //-T
            return ('t');
        else if (allele2 == ADO)    //-?
            return ('-');
        else if (allele2 == DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}

/********************* WhichMut ************************/
/* Returns character representation for binary data */

char WhichMut (int state)
{
    if (state == 0)
        return ('0');
    else if (state == 1)
        return ('1');
    else if (state == ADO)
        return ('?');
    else if (state == DELETION)
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

char WhichConsensusBinary (int allele1, int allele2)
{
    if (allele1 == 0)
    {
        if (allele2 == 0)        //00
            return ('0');
        else if (allele2 == 1)    //01
            return ('1');
        else if (allele2 == ADO)    //0?
            return ('0');
        else if (allele2 == DELETION)    //0-
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
        else if (allele2 == ADO)    //1?
            return ('2');
        else if (allele2 == DELETION)    //0-
            return ('2');
        else
            return ('N');
    }
    else if (allele1 == ADO)
    {
        if (allele2 == 0)        //?0
            return ('0');
        else if (allele2 == 1)    //?1
            return ('2');
        else if (allele2 == ADO)    //??
            return ('-');
        else if (allele2 == DELETION)    //?-
            return ('-');
        else
            return ('N');
    }
    else if (allele1 == DELETION)
    {
        if (allele2 == 0)        //-0
            return ('0');
        else if (allele2 == 1)    //-1
            return ('2');
        else if (allele2 == ADO)    //-?
            return ('-');
        else if (allele2 == DELETION)    //--
            return ('-');
        else
            return ('N');
    }
    else
        return ('N');
}

/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char WhichNuc (int nucleotide)
{
    if (nucleotide == A)
        return ('A');
    else if (nucleotide == C)
        return ('C');
    else if (nucleotide == G)
        return ('G');
    else if (nucleotide == T)
        return ('T');
    else if (nucleotide == ADO)
        return ('?');
    else if (nucleotide == DELETION)
        return ('-');
    else
        return ('N');
}

/************************ CheckMatrixSymmetry **************************/
/* Checks whether a given matrix is symmetric */

int CheckMatrixSymmetry(double matrix[4][4])
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

int WhichNucChar (char nucleotide)
{
    if (nucleotide == 'A')
        return (A);
    else if (nucleotide == 'C')
        return (C);
    else if (nucleotide == 'G')
        return (G);
    else if (nucleotide == 'T')
        return (T);
    else if (nucleotide == '?')
        return (ADO);
    else if (nucleotide == '-')
        return (DELETION);
    else if (nucleotide == 'N')
        return (N);
    else if (nucleotide == 'R')
        return (R);
    else
    {
        fprintf (stderr, "\nERROR in WhichNucChar: nucleotide = %c\n",  nucleotide);
        exit(-1);
    }
}

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

//  eigein.h
//  CellCoal
//
//  Created by David Posada on 21/11/2017.
//  Everything below is shamelessly taken from Yang's Paml package */

//

#ifndef eigein_h
#define eigein_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int EigenREV (double mr, double Qij[], double Root[], double Cijk[]);
//double Qij[16], mr;

/* Everything below is shamelessly taken from Yang's Paml package */
int abyx (double a, double x[], int n);
int xtoy (double x[], double y[], int n);
int matinv( double x[], int n, int m, double space[]);
int eigen(int job, double AAA[], int n, double rr[], double ri[], double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi, double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[], double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi, double vr[], double vi[], int work[]);

//typedef struct { double re, im; } complex;
class complex
{
public:
    double re;
    double im;
    complex(double re, double im);
    complex conjj (complex a);
    complex cplus (complex a, complex b);
    complex cminus (complex a, complex b);
    complex cby (complex a, complex b);
    static complex cdiv (complex a,complex b);
    complex cexpp (complex a);
    complex cfactor (complex x, double a);
    int cxtoy (complex x[], complex y[], int n);
    int cmatby (complex a[], complex b[], complex c[], int n,int m,int k);
    int cmatout (FILE * fout, complex x[], int n, int m);
    int cmatinv( complex x[], int n, int m, double space[]);
};
#define csize(a) (fabs(a.re)+fabs(a.im))

#endif /* eigein_h */

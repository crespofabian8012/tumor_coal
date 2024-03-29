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

//  eigen.c
//  CellCoal
//
//  Created by David Posada on 21/11/2017.
//  Everything below is shamelessly taken from Yang's Paml package */
//
#include "eigen.hpp"
namespace linalgebra{
int EigenREV (double mr, double Qij[], double Root[], double Cijk[])
{
    /* freq[] is constant
     */
    int i,j,k;
    double U[16], V[16], T1[16], T2[16];
    
    abyx (1/mr, Qij, 16);
    
    if ((k=eigen (1, Qij, 4, Root, T1, U, V, T2))!=0) {
        fprintf(stderr, "\ncomplex roots in EigenREV");
        exit(0);
    }
    xtoy (U, V, 16);
    matinv (V, 4, 4, T1);
    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
            for (k=0; k<4; k++)
                Cijk[i*4*4+j*4+k] = U[i*4+k]*V[k*4+j];
    return (0);
}

int abyx (double a, double x[], int n)
{
    int i; for (i=0; i<n; x[i]*=a,i++);
    return(0);
}
int xtoy (double x[], double y[], int n)
{
    int i;
    for (i=0; i<n; y[i]=x[i],i++);
    return(0);
}

int matinv( double x[], int n, int m, double space[])
{
    /* x[n*m]  ... m>=n
     */
    // register int i,j,k;
    int i,j,k;
    int *irow=(int*) space;
    double ee=1.0e-20, t,t1,xmax;
    double det=1.0;
    
    for (i=0; i<n; i++)  {
        xmax = 0.;
        for (j=i; j<n; j++) {
            if (xmax < fabs(x[j*m+i]))  {
                xmax = fabs( x[j*m+i] );
                irow[i] = j;
            }
        }
        det *= xmax;
        if (xmax < ee)   {
            fprintf(stderr,"\nDet becomes zero at %3d!\t\n", i+1);
            return(-1);
        }
        if (irow[i] != i) {
            for (j=0; j<m; j++) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i] * m + j];
                x[ irow[i] * m + j] = t;
            }
        }
        t = 1./x[i*m+i];
        for (j=0; j<n; j++) {
            if (j == i) continue;
            t1 = t*x[j*m+i];
            for (k=0; k<n; k++)  x[j*m+k] -= t1*x[i*m+k];
            x[j*m+i] = -t1;
        }
        for (j=0; j<m; j++)   x[i*m+j] *= t;
        x[i*m+i] = t;
    }                            /* i  */
    for (i=n-1; i>=0; i--) {
        if (irow[i] == i) continue;
        for (j=0; j<n; j++)  {
            t = x[j*m+i];
            x[j*m+i] = x[ j*m + irow[i] ];
            x[ j*m + irow[i] ] = t;
        }
    }
    return (0);
}

/***********************************************************
 *  This eigen() works for eigenvalue/vector analysis
 *         for real general square matrix A
 *         A will be destroyed
 *         rr,ri are vectors containing eigenvalues
 *         vr,vi are matrices containing (right) eigenvectors
 *
 *              A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
 *
 *  Algorithm: Handbook for Automatic Computation, vol 2
 *             by Wilkinson and Reinsch, 1971
 *             most of source codes were taken from a public domain
 *             solftware called MATCALC.
 *  Credits:   to the authors of MATCALC
 *
 *  return     -1 not converged
 *              0 no complex eigenvalues/vectors
 *              1 complex eigenvalues/vectors
 *  Tianlin Wang at University of Illinois
 *  Thu May  6 15:22:31 CDT 1993
 ***************************************************************/

#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     53    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */

#define pos(i,j,n)      ((i)*(n)+(j))

int eigen(int job, double AAA[], int n, double rr[], double ri[], double vr[], double vi[], double work[])
{
    /*  double work[n*2]: working space
     */
    int low,hi,i,j,k, it, istate=0;
    double tiny=sqrt(pow((double)BASE,(double)(1-DIGITS))), t;
    
    balance(AAA,n,&low,&hi,work);
    elemhess(job,AAA,n,low,hi,vr,vi, (int*)(work+n));
    if (-1 == realeig(job,AAA,n,low,hi,rr,ri,vr,vi)) return (-1);
    if (job) unbalance(n,vr,vi,low,hi,work);
    
    /* sort, added by Z. Yang */
    for (i=0; i<n; i++) {
        for (j=i+1,it=i,t=rr[i]; j<n; j++)
            if (t<rr[j]) { t=rr[j]; it=j; }
        rr[it]=rr[i];   rr[i]=t;
        t=ri[it];       ri[it]=ri[i];  ri[i]=t;
        for (k=0; k<n; k++) {
            t=vr[k*n+it];  vr[k*n+it]=vr[k*n+i];  vr[k*n+i]=t;
            t=vi[k*n+it];  vi[k*n+it]=vi[k*n+i];  vi[k*n+i]=t;
        }
        if (fabs(ri[i])>tiny) istate=1;
    }
    
    return (istate) ;
}

/* complex funcctions
 */

complex::complex(double re,double im)
{
    this->re = re;
    this->im = im;
}

complex::complex(int re){
    this->re = (double) (1.0 *re);
    this->im = (double)0.0;
}
complex complex::conjj (complex a)
{
    a.im = -a.im;
    return(a);
}

#define csize(a) (fabs(a.re)+fabs(a.im))

complex linalgebra::complex::cplus (complex a, complex b)
{
    double re = a.re+b.re;
    double im = a.im+b.im;
    complex c(re, im);
    return (c);
}

complex linalgebra::complex::cminus (complex a, complex b)
{
    double re = a.re-b.re;
    double im = a.im-b.im;
    complex c(re, im);
    return (c);
}

complex linalgebra::complex::cby (complex a, complex b)
{
    double re = a.re*b.re-a.im*b.im ;
    double im = a.re*b.im+a.im*b.re ;
    complex c(re, im);
    return (c);
}

complex linalgebra::complex::cdiv (complex a,complex b)
{
    double ratio, den;
    complex c(0, 0);
    
    if (fabs(b.re) <= fabs(b.im)) {
        ratio = b.re / b.im;
        den = b.im * (1 + ratio * ratio);
        c.re = (a.re * ratio + a.im) / den;
        c.im = (a.im * ratio - a.re) / den;
    }
    else {
        ratio = b.im / b.re;
        den = b.re * (1 + ratio * ratio);
        c.re = (a.re + a.im * ratio) / den;
        c.im = (a.im - a.re * ratio) / den;
    }
    return(c);
}

complex linalgebra::complex::cexpp (complex a)
{
    complex c(0, 0);
    c.re = exp(a.re);
    if (fabs(a.im)==0) c.im = 0;
    else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
    return (c);
}

complex linalgebra::complex::cfactor (complex x, double a)
{
    complex c(0,0);
    c.re = a*x.re;
    c.im = a*x.im;
    return (c);
}

int linalgebra::complex::cxtoy (complex x[], complex y[], int n)
{
    int i;
    FOR (i,n) y[i]=x[i];
    return (0);
}

int linalgebra::complex::cmatby (complex a[], complex b[], complex c[], int n,int m,int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
 */
{
    int i,j,i1;
    complex t(0,0);
    
    FOR (i,n)  FOR(j,k) {
        for (i1=0,t=complex(0,0); i1<m; i1++)
            t = cplus (t, cby(a[i*m+i1],b[i1*k+j]));
        c[i*k+j] = t;
    }
    return (0);
}

int linalgebra::complex::cmatout (FILE * fout, complex x[], int n, int m)
{
    int i,j;
    for (i=0,FPN(fout); i<n; i++,FPN(fout))
        FOR(j,m) fprintf(fout, "%7.3f%7.3f  ", x[i*m+j].re, x[i*m+j].im);
    return (0);
}

int linalgebra::complex::cmatinv( complex x[], int n, int m, double space[])
{
    /* x[n*m]  ... m>=n
     */
    int i,j,k, *irow=(int*) space;
    double xmaxsize, ee=1e-20;
    complex xmax(0,0), t(0,0), t1(0,0);
    
    FOR(i,n)  {
        xmaxsize = 0.;
        for (j=i; j<n; j++) {
            if ( xmaxsize < csize (x[j*m+i]))  {
                xmaxsize = csize (x[j*m+i]);
                xmax = x[j*m+i];
                irow[i] = j;
            }
        }
        if (xmaxsize < ee)   {
            fprintf(stderr,"\nDet goes to zero at %8d!\t\n", i+1);
            return(-1);
        }
        if (irow[i] != i) {
            FOR(j,m) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i]*m+j];
                x[ irow[i]*m+j] = t;
            }
        }
        t = cdiv (complex(1,0), x[i*m+i]);
        FOR(j,n) {
            if (j == i) continue;
            t1 = cby (t,x[j*m+i]);
            FOR(k,m)  x[j*m+k] = cminus (x[j*m+k], cby(t1,x[i*m+k]));
            x[j*m+i] = cfactor (t1, -1);
        }
        FOR(j,m)   x[i*m+j] = cby (x[i*m+j], t);
        x[i*m+i] = t;
    }
    for (i=n-1; i>=0; i--) {
        if (irow[i] == i) continue;
        FOR(j,n)  {
            t = x[j*m+i];
            x[j*m+i] = x[j*m+irow[i]];
            x[ j*m+irow[i]] = t;
        }
    }
    return (0);
}


void balance(double mat[], int n,int *low, int *hi, double scale[])
{
    /* Balance a matrix for calculation of eigenvalues and eigenvectors
     */
    double c,f,g,r,s;
    int i,j,k,l,done;
    /* search for rows isolating an eigenvalue and push them down */
    for (k = n - 1; k >= 0; k--) {
        for (j = k; j >= 0; j--) {
            for (i = 0; i <= k; i++) {
                if (i != j && fabs(mat[pos(j,i,n)]) != 0) break;
            }
            
            if (i > k) {
                scale[k] = j;
                
                if (j != k) {
                    for (i = 0; i <= k; i++) {
                        c = mat[pos(i,j,n)];
                        mat[pos(i,j,n)] = mat[pos(i,k,n)];
                        mat[pos(i,k,n)] = c;
                    }
                    
                    for (i = 0; i < n; i++) {
                        c = mat[pos(j,i,n)];
                        mat[pos(j,i,n)] = mat[pos(k,i,n)];
                        mat[pos(k,i,n)] = c;
                    }
                }
                break;
            }
        }
        if (j < 0) break;
    }
    
    /* search for columns isolating an eigenvalue and push them left */
    
    for (l = 0; l <= k; l++) {
        for (j = l; j <= k; j++) {
            for (i = l; i <= k; i++) {
                if (i != j && fabs(mat[pos(i,j,n)]) != 0) break;
            }
            if (i > k) {
                scale[l] = j;
                if (j != l) {
                    for (i = 0; i <= k; i++) {
                        c = mat[pos(i,j,n)];
                        mat[pos(i,j,n)] = mat[pos(i,l,n)];
                        mat[pos(i,l,n)] = c;
                    }
                    
                    for (i = l; i < n; i++) {
                        c = mat[pos(j,i,n)];
                        mat[pos(j,i,n)] = mat[pos(l,i,n)];
                        mat[pos(l,i,n)] = c;
                    }
                }
                
                break;
            }
        }
        
        if (j > k) break;
    }
    
    *hi = k;
    *low = l;
    
    /* balance the submatrix in rows l through k */
    
    for (i = l; i <= k; i++) {
        scale[i] = 1;
    }
    
    do {
        for (done = 1,i = l; i <= k; i++) {
            for (c = 0,r = 0,j = l; j <= k; j++) {
                if (j != i) {
                    c += fabs(mat[pos(j,i,n)]);
                    r += fabs(mat[pos(i,j,n)]);
                }
            }
            
            if (c != 0 && r != 0) {
                g = r / BASE;
                f = 1;
                s = c + r;
                
                while (c < g) {
                    f *= BASE;
                    c *= BASE * BASE;
                }
                
                g = r * BASE;
                
                while (c >= g) {
                    f /= BASE;
                    c /= BASE * BASE;
                }
                
                if ((c + r) / f < 0.95 * s) {
                    done = 0;
                    g = 1 / f;
                    scale[i] *= f;
                    
                    for (j = l; j < n; j++) {
                        mat[pos(i,j,n)] *= g;
                    }
                    
                    for (j = 0; j <= k; j++) {
                        mat[pos(j,i,n)] *= f;
                    }
                }
            }
        }
    } while (!done);
}


/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void unbalance(int n,double vr[],double vi[], int low, int hi, double scale[])
{
    int i,j,k;
    double tmp;
    
    for (i = low; i <= hi; i++) {
        for (j = 0; j < n; j++) {
            vr[pos(i,j,n)] *= scale[i];
            vi[pos(i,j,n)] *= scale[i];
        }
    }
    
    for (i = low - 1; i >= 0; i--) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[pos(i,j,n)];
                vr[pos(i,j,n)] = vr[pos(k,j,n)];
                vr[pos(k,j,n)] = tmp;
                
                tmp = vi[pos(i,j,n)];
                vi[pos(i,j,n)] = vi[pos(k,j,n)];
                vi[pos(k,j,n)] = tmp;
            }
        }
    }
    
    for (i = hi + 1; i < n; i++) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[pos(i,j,n)];
                vr[pos(i,j,n)] = vr[pos(k,j,n)];
                vr[pos(k,j,n)] = tmp;
                
                tmp = vi[pos(i,j,n)];
                vi[pos(i,j,n)] = vi[pos(k,j,n)];
                vi[pos(k,j,n)] = tmp;
            }
        }
    }
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void elemhess(int job,double mat[],int n,int low,int hi, double vr[],
              double vi[], int work[])
{
    /* work[n] */
    int i,j,m;
    double x,y;
    
    for (m = low + 1; m < hi; m++) {
        for (x = 0,i = m,j = m; j <= hi; j++) {
            if (fabs(mat[pos(j,m-1,n)]) > fabs(x)) {
                x = mat[pos(j,m-1,n)];
                i = j;
            }
        }
        
        if ((work[m] = i) != m) {
            for (j = m - 1; j < n; j++) {
                y = mat[pos(i,j,n)];
                mat[pos(i,j,n)] = mat[pos(m,j,n)];
                mat[pos(m,j,n)] = y;
            }
            
            for (j = 0; j <= hi; j++) {
                y = mat[pos(j,i,n)];
                mat[pos(j,i,n)] = mat[pos(j,m,n)];
                mat[pos(j,m,n)] = y;
            }
        }
        
        if (x != 0) {
            for (i = m + 1; i <= hi; i++) {
                if ((y = mat[pos(i,m-1,n)]) != 0) {
                    y = mat[pos(i,m-1,n)] = y / x;
                    
                    for (j = m; j < n; j++) {
                        mat[pos(i,j,n)] -= y * mat[pos(m,j,n)];
                    }
                    
                    for (j = 0; j <= hi; j++) {
                        mat[pos(j,m,n)] += y * mat[pos(j,i,n)];
                    }
                }
            }
        }
    }
    if (job) {
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                vr[pos(i,j,n)] = 0.0; vi[pos(i,j,n)] = 0.0;
            }
            vr[pos(i,i,n)] = 1.0;
        }
        
        for (m = hi - 1; m > low; m--) {
            for (i = m + 1; i <= hi; i++) {
                vr[pos(i,m,n)] = mat[pos(i,m-1,n)];
            }
            
            if ((i = work[m]) != m) {
                for (j = m; j <= hi; j++) {
                    vr[pos(m,j,n)] = vr[pos(i,j,n)];
                    vr[pos(i,j,n)] = 0.0;
                }
                vr[pos(i,m,n)] = 1.0;
            }
        }
    }
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */

int realeig(int job,double mat[],int n,int low, int hi, double valr[],
            double vali[], double vr[],double vi[])
{
    complex v(0,0);
    double p=0,q=0,r=0,s=0,t,w,x,y,z=0,ra,sa,norm,eps;
    int niter,en,i,j,k,l,m;
    double precision  = pow((double)BASE,(double)(1-DIGITS));
    
    eps = precision;
    for (i=0; i<n; i++) {
        valr[i]=0.0;
        vali[i]=0.0;
    }
    /* store isolated roots and calculate norm */
    for (norm = 0,i = 0; i < n; i++) {
        for (j = max(0,i-1); j < n; j++) {
            norm += fabs(mat[pos(i,j,n)]);
        }
        if (i < low || i > hi) valr[i] = mat[pos(i,i,n)];
    }
    t = 0;
    en = hi;
    
    while (en >= low) {
        niter = 0;
        for (;;) {
            
            /* look for single small subdiagonal element */
            
            for (l = en; l > low; l--) {
                s = fabs(mat[pos(l-1,l-1,n)]) + fabs(mat[pos(l,l,n)]);
                if (s == 0) s = norm;
                if (fabs(mat[pos(l,l-1,n)]) <= eps * s) break;
            }
            
            /* form shift */
            
            x = mat[pos(en,en,n)];
            
            if (l == en) {             /* one root found */
                valr[en] = x + t;
                if (job) mat[pos(en,en,n)] = x + t;
                en--;
                break;
            }
            
            y = mat[pos(en-1,en-1,n)];
            w = mat[pos(en,en-1,n)] * mat[pos(en-1,en,n)];
            
            if (l == en - 1) {                /* two roots found */
                p = (y - x) / 2;
                q = p * p + w;
                z = sqrt(fabs(q));
                x += t;
                if (job) {
                    mat[pos(en,en,n)] = x;
                    mat[pos(en-1,en-1,n)] = y + t;
                }
                if (q < 0) {                /* complex pair */
                    valr[en-1] = x+p;
                    vali[en-1] = z;
                    valr[en] = x+p;
                    vali[en] = -z;
                }
                else {                      /* real pair */
                    z = (p < 0) ? p - z : p + z;
                    valr[en-1] = x + z;
                    valr[en] = (z == 0) ? x + z : x - w / z;
                    if (job) {
                        x = mat[pos(en,en-1,n)];
                        s = fabs(x) + fabs(z);
                        p = x / s;
                        q = z / s;
                        r = sqrt(p*p+q*q);
                        p /= r;
                        q /= r;
                        for (j = en - 1; j < n; j++) {
                            z = mat[pos(en-1,j,n)];
                            mat[pos(en-1,j,n)] = q * z + p *
                            mat[pos(en,j,n)];
                            mat[pos(en,j,n)] = q * mat[pos(en,j,n)] - p*z;
                        }
                        for (i = 0; i <= en; i++) {
                            z = mat[pos(i,en-1,n)];
                            mat[pos(i,en-1,n)] = q * z + p * mat[pos(i,en,n)];
                            mat[pos(i,en,n)] = q * mat[pos(i,en,n)] - p*z;
                        }
                        for (i = low; i <= hi; i++) {
                            z = vr[pos(i,en-1,n)];
                            vr[pos(i,en-1,n)] = q*z + p*vr[pos(i,en,n)];
                            vr[pos(i,en,n)] = q*vr[pos(i,en,n)] - p*z;
                        }
                    }
                }
                en -= 2;
                break;
            }
            if (niter == MAXITER) return(-1);
            if (niter != 0 && niter % 10 == 0) {
                t += x;
                for (i = low; i <= en; i++) mat[pos(i,i,n)] -= x;
                s = fabs(mat[pos(en,en-1,n)]) + fabs(mat[pos(en-1,en-2,n)]);
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }
            niter++;
            /* look for two consecutive small subdiagonal elements */
            for (m = en - 2; m >= l; m--) {
                z = mat[pos(m,m,n)];
                r = x - z;
                s = y - z;
                p = (r * s - w) / mat[pos(m+1,m,n)] + mat[pos(m,m+1,n)];
                q = mat[pos(m+1,m+1,n)] - z - r - s;
                r = mat[pos(m+2,m+1,n)];
                s = fabs(p) + fabs(q) + fabs(r);
                p /= s;
                q /= s;
                r /= s;
                if (m == l || fabs(mat[pos(m,m-1,n)]) * (fabs(q)+fabs(r)) <=
                    eps * (fabs(mat[pos(m-1,m-1,n)]) + fabs(z) +
                           fabs(mat[pos(m+1,m+1,n)])) * fabs(p)) break;
            }
            for (i = m + 2; i <= en; i++) mat[pos(i,i-2,n)] = 0;
            for (i = m + 3; i <= en; i++) mat[pos(i,i-3,n)] = 0;
            /* double QR step involving rows l to en and columns m to en */
            for (k = m; k < en; k++) {
                if (k != m) {
                    p = mat[pos(k,k-1,n)];
                    q = mat[pos(k+1,k-1,n)];
                    r = (k == en - 1) ? 0 : mat[pos(k+2,k-1,n)];
                    if ((x = fabs(p) + fabs(q) + fabs(r)) == 0) continue;
                    p /= x;
                    q /= x;
                    r /= x;
                }
                s = sqrt(p*p+q*q+r*r);
                if (p < 0) s = -s;
                if (k != m) {
                    mat[pos(k,k-1,n)] = -s * x;
                }
                else if (l != m) {
                    mat[pos(k,k-1,n)] = -mat[pos(k,k-1,n)];
                }
                p += s;
                x = p / s;
                y = q / s;
                z = r / s;
                q /= p;
                r /= p;
                /* row modification */
                for (j = k; j <= (!job ? en : n-1); j++){
                    p = mat[pos(k,j,n)] + q * mat[pos(k+1,j,n)];
                    if (k != en - 1) {
                        p += r * mat[pos(k+2,j,n)];
                        mat[pos(k+2,j,n)] -= p * z;
                    }
                    mat[pos(k+1,j,n)] -= p * y;
                    mat[pos(k,j,n)] -= p * x;
                }
                j = min(en,k+3);
                /* column modification */
                for (i = (!job ? l : 0); i <= j; i++) {
                    p = x * mat[pos(i,k,n)] + y * mat[pos(i,k+1,n)];
                    if (k != en - 1) {
                        p += z * mat[pos(i,k+2,n)];
                        mat[pos(i,k+2,n)] -= p*r;
                    }
                    mat[pos(i,k+1,n)] -= p*q;
                    mat[pos(i,k,n)] -= p;
                }
                if (job) {             /* accumulate transformations */
                    for (i = low; i <= hi; i++) {
                        p = x * vr[pos(i,k,n)] + y * vr[pos(i,k+1,n)];
                        if (k != en - 1) {
                            p += z * vr[pos(i,k+2,n)];
                            vr[pos(i,k+2,n)] -= p*r;
                        }
                        vr[pos(i,k+1,n)] -= p*q;
                        vr[pos(i,k,n)] -= p;
                    }
                }
            }
        }
    }
    
    if (!job) return(0);
    if (norm != 0) {
        /* back substitute to find vectors of upper triangular form */
        for (en = n-1; en >= 0; en--) {
            p = valr[en];
            if ((q = vali[en]) < 0) {            /* complex vector */
                m = en - 1;
                if (fabs(mat[pos(en,en-1,n)]) > fabs(mat[pos(en-1,en,n)])) {
                    mat[pos(en-1,en-1,n)] = q / mat[pos(en,en-1,n)];
                    mat[pos(en-1,en,n)] = (p - mat[pos(en,en,n)]) /
                    mat[pos(en,en-1,n)];
                }
                else {
                    v = complex::cdiv(complex(0.0,-mat[pos(en-1,en,n)]),
                                      complex(mat[pos(en-1,en-1,n)]-p,q));
                    mat[pos(en-1,en-1,n)] = v.re;
                    mat[pos(en-1,en,n)] = v.im;
                }
                mat[pos(en,en-1,n)] = 0;
                mat[pos(en,en,n)] = 1;
                for (i = en - 2; i >= 0; i--) {
                    w = mat[pos(i,i,n)] - p;
                    ra = 0;
                    sa = mat[pos(i,en,n)];
                    for (j = m; j < en; j++) {
                        ra += mat[pos(i,j,n)] * mat[pos(j,en-1,n)];
                        sa += mat[pos(i,j,n)] * mat[pos(j,en,n)];
                    }
                    if (vali[i] < 0) {
                        z = w;
                        r = ra;
                        s = sa;
                    }
                    else {
                        m = i;
                        if (vali[i] == 0) {
                            v = complex::cdiv(complex(-ra,-sa),complex(w,q));
                            mat[pos(i,en-1,n)] = v.re;
                            mat[pos(i,en,n)] = v.im;
                        }
                        else {                      /* solve complex equations */
                            x = mat[pos(i,i+1,n)];
                            y = mat[pos(i+1,i,n)];
                            v.re = (valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q;
                            v.im = (valr[i] - p)*2*q;
                            if ((fabs(v.re) + fabs(v.im)) == 0) {
                                v.re = eps * norm * (fabs(w) +
                                                     fabs(q) + fabs(x) + fabs(y) + fabs(z));
                            }
                            v = complex::cdiv(complex(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
                            mat[pos(i,en-1,n)] = v.re;
                            mat[pos(i,en,n)] = v.im;
                            if (fabs(x) > fabs(z) + fabs(q)) {
                                mat[pos(i+1,en-1,n)] =
                                (-ra - w * mat[pos(i,en-1,n)] +
                                 q * mat[pos(i,en,n)]) / x;
                                mat[pos(i+1,en,n)] = (-sa - w * mat[pos(i,en,n)] -
                                                      q * mat[pos(i,en-1,n)]) / x;
                            }
                            else {
                                v = complex::cdiv(complex(-r-y*mat[pos(i,en-1,n)],
                                                          -s-y*mat[pos(i,en,n)]),complex(z,q));
                                mat[pos(i+1,en-1,n)] = v.re;
                                mat[pos(i+1,en,n)] = v.im;
                            }
                        }
                    }
                }
            }
            else if (q == 0) {                             /* real vector */
                m = en;
                mat[pos(en,en,n)] = 1;
                for (i = en - 1; i >= 0; i--) {
                    w = mat[pos(i,i,n)] - p;
                    r = mat[pos(i,en,n)];
                    for (j = m; j < en; j++) {
                        r += mat[pos(i,j,n)] * mat[pos(j,en,n)];
                    }
                    if (vali[i] < 0) {
                        z = w;
                        s = r;
                    }
                    else {
                        m = i;
                        if (vali[i] == 0) {
                            if ((t = w) == 0) t = eps * norm;
                            mat[pos(i,en,n)] = -r / t;
                        }
                        else {            /* solve real equations */
                            x = mat[pos(i,i+1,n)];
                            y = mat[pos(i+1,i,n)];
                            q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
                            t = (x * s - z * r) / q;
                            mat[pos(i,en,n)] = t;
                            if (fabs(x) <= fabs(z)) {
                                mat[pos(i+1,en,n)] = (-s - y * t) / z;
                            }
                            else {
                                mat[pos(i+1,en,n)] = (-r - w * t) / x;
                            }
                        }
                    }
                }
            }
        }
        /* vectors of isolated roots */
        for (i = 0; i < n; i++) {
            if (i < low || i > hi) {
                for (j = i; j < n; j++) {
                    vr[pos(i,j,n)] = mat[pos(i,j,n)];
                }
            }
        }
        /* multiply by transformation matrix */
        
        for (j = n-1; j >= low; j--) {
            m = min(j,hi);
            for (i = low; i <= hi; i++) {
                for (z = 0,k = low; k <= m; k++) {
                    z += vr[pos(i,k,n)] * mat[pos(k,j,n)];
                }
                vr[pos(i,j,n)] = z;
            }
        }
    }
    /* rearrange complex eigenvectors */
    for (j = 0; j < n; j++) {
        if (vali[j] != 0) {
            for (i = 0; i < n; i++) {
                vi[pos(i,j,n)] = vr[pos(i,j+1,n)];
                vr[pos(i,j+1,n)] = vr[pos(i,j,n)];
                vi[pos(i,j+1,n)] = -vi[pos(i,j,n)];
            }
            j++;
        }
    }
    return(0);
}
}


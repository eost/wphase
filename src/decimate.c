/***************************************************************************
*
*                     W phase source inversion package              
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Universite de Strasbourg / CNRS 
*                                  April 2013
*
*    Neither the name of the California Institute of Technology (Caltech) 
*    nor the names of its contributors may be used to endorse or promote 
*    products derived from this software without specific prior written 
*    permission
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
****************************************************************************/

/*      Decimation      */

#include<stdio.h>
#include<stdlib.h>
#include"proto_alloc.h"
#define  PARYTY_EVEN 0

typedef struct
{
    int     n1,n2  ;
    double *coeffs ;
} FIR_filter ;

typedef struct
{
    int counter     ;
    double *segment ;
} channel ;

typedef struct{
    int   ratio, n, *facs;
} Cascade;

void error(const char* p,const char* p2)
{
    fprintf(stderr,"%s %s\n",p,p2);
    exit(1);
}

double dot(double *a,int na,double *b,int nb)
{
    int j;
    if (na != nb) error("Error: length mismatch in 'dot'\n", "");
    double val = 0.;
    for (j=0; j<na; j++)
        val += a[j]*b[j];
    return val;
}

int init_FIR(double *coeffs,int Ncoeffs,FIR_filter *FIR)
{
    int    i,j,N;
    double anorm=0.;
    N = Ncoeffs*2-1;
    if (PARYTY_EVEN)
        N += 1;
    FIR->coeffs = double_alloc(N);
    for(i=0;i<Ncoeffs;i++)
    {
        FIR->coeffs[i] = coeffs[i];
        anorm += coeffs[i];
    }
    FIR->n2  = Ncoeffs;
    anorm   *= 2.;
    anorm   -= FIR->coeffs[Ncoeffs-1];
    FIR->n1  = FIR->n2 - 1;
    if (PARYTY_EVEN)
    {
        FIR->n1 += 1;
        anorm   += FIR->coeffs[Ncoeffs-1];
    }
    for(j=0; j < FIR->n2; j++)
        FIR->coeffs[j] /= anorm;
    for(j=0; j < FIR->n1; j++)
        FIR->coeffs[j+FIR->n2] = FIR->coeffs[FIR->n1-1-j];
    return 0;
}  

void init_chan(channel *chan,FIR_filter *FIR)
{
    chan->counter = 0;
    chan->segment = double_calloc(FIR->n1+FIR->n2);
}

void push_back(double *v, int nv, double back)
{
    int i;
    for (i=0;i<nv-1;i++)
        v[i] = v[i+1];
    v[nv-1] = back;
}

int decimate(FIR_filter *FIR, int dec_fac, double *yi, int ni, double *yo, int *no)
{
    int i,j,N;
    double sample_value ;
    channel chan ;

    if ( dec_fac == 1 ) return 0;
    init_chan(&chan,FIR);
    j = 0 ;
    N = FIR->n1+FIR->n2 ;
    for(i=FIR->n2;i<N;i++)
        chan.segment[i]=yi[i-FIR->n2];
    for(i=FIR->n1;i<ni+FIR->n1;i++)
    { 
        if (i<ni)
            sample_value = yi[i] ;
        else
            sample_value = 0. ;
        push_back(chan.segment,N,sample_value);
        if(chan.counter >= 0 && chan.counter % dec_fac == 0)
        {
            yo[j] = dot(FIR->coeffs,N,chan.segment,N);
            chan.counter = 0;
            j++;
        }
        ++(chan.counter);
    }
    (*no)=j;
    return 0;
}

void init_casc(Cascade *cascade, int N, int *facs)
{
//  const int N=3;
    int k;
    (*cascade).n = N;
    (*cascade).facs = (int *)calloc(N, sizeof(int));
    (*cascade).ratio = 1;
    for (k=0; k<N; k++)
    {
        (*cascade).facs[k] = facs[k];
        (*cascade).ratio *= facs[k];
    }
}

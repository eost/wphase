/****************************************************************
*	W phase package - Decimation
*                                           
*       History
*         JAN 2010  Original Coding in C++    Luis Rivera
*         OCT 2010  Further Updates           Zacharie Duputel
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include"proto_alloc.h"
#define  PARYTY_EVEN 0

typedef struct{
  int     n1,n2  ;
  double *coeffs ;
} FIR_filter ;

typedef struct{
  int counter     ;
  double *segment ;
} channel ;

void 
error(const char* p,const char* p2)
{
  fprintf(stderr,"%s %s\n",p,p2);
  exit(1);
}

double 
dot(double *a,int na,double *b,int nb)
{
  int j;
  if (na != nb) error("Error: length mismatch in 'dot'\n", "");
  double val = 0.;
  for (j=0; j<na; j++)
	val += a[j]*b[j];
  return val;
}

int 
init_FIR(double *coeffs,int Ncoeffs,FIR_filter *FIR)
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
    FIR->coeffs[j+FIR->n2] = coeffs[FIR->n1-1-j];
  return 0;
}  

void
init_chan(channel *chan,FIR_filter *FIR)
{
  chan->counter = 0;
  chan->segment = double_calloc(FIR->n1+FIR->n2);
}

void
push_back(double *v, int nv, double back)
{
  int i;
  for (i=0;i<nv-1;i++)
    v[i] = v[i+1];
  v[nv-1] = back;
}

int 
decimate(FIR_filter *FIR, int dec_fac, double *yi, int ni, double *yo, int *no)
{
  int i,j,N;
  double sample_value ;
  channel chan ;
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



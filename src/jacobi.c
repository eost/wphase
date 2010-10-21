/****************************************************************
*	W phase package - jacobi 
*      (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1. 
*                         
*       History
*             2010   adapted for double precision by Z.Duputel   
*       License
*             Distributed under the same License as the Numerical
*             Recipes source and binaries, used here under permission 
*             by Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proto_alloc.h"  /* subroutines in proto_alloc.c */
#define ROTATE(a,i,j,k,l) g=a[i][j] ;\
                          h=a[k][l] ;\
                          a[i][j]=g-s*(h+g*tau) ;\
                          a[k][l]=h+s*(g-h*tau) ;

void 
jacobi(double **a, int n, int np, double *d, double **v, int *nrot)
{
  int    NMAX = 500, i, ip, iq, j;
  double c,g,h,s,sm,t,tau,theta,tresh, *b, *z ;

  b = double_alloc(NMAX) ;
  z = double_calloc(NMAX) ;


  /* Set Diagonal matrix */
  for (ip = 0 ; ip<n ; ip++)
    {
      for (iq=0 ; iq<n ; iq++)
	v[ip][iq] = 0. ;
      v[ip][ip] = 1.   ;
    }
  for(ip=0 ; ip<n ; ip++)
    {
      b[ip] = d[ip] = a[ip][ip] ;
      z[ip] = 0.0 ;
    }
  
  *nrot=0 ;
  for(i=1 ; i<=50 ; i++) 
    {
      sm = 0.0 ;
      for (ip=0 ; ip<n-1 ; ip++)
	for (iq=ip+1 ; iq<n ; iq++)
	  sm = sm + fabs(a[ip][iq]) ;
      if (sm == 0.0)
	{
	  free((void *)z) ;
	  free((void *)b) ;
	  return;
	}
      if(i < 4)
	tresh = 0.2 * sm / (n*n) ;
      else
	tresh=0.0 ;
      for (ip = 0 ; ip<n-1 ; ip++)//22
	{
	  for (iq=ip+1 ; iq<n ; iq++) //21
	    {
	      g = 100.0 * fabs(a[ip][iq]) ;
	      if((i>4) && (float)(fabs(d[ip])+g) == (float)(fabs(d[ip]))
		 && (float)(fabs(d[iq])+g) == (float)(fabs(d[iq])))
		a[ip][iq] = 0.0 ;
	      else if(fabs(a[ip][iq]) > tresh)
		{
		  h = d[iq]-d[ip]  ;
		  if((float)(fabs(h)+g) == (float)(fabs(h)))
		    t = a[ip][iq]/h  ;
		  else
		    {
		      theta = 0.5 * h/(a[ip][iq]) ;
		      t     = 1./(fabs(theta) + sqrt(1.0+theta*theta));
		      if(theta < 0.) t = -t  ;
		    }
		  c   = 1./sqrt(1+t*t) ;
		  s   = t * c       ;
		  tau = s/(1.0+c)   ;
		  h   = t*a[ip][iq] ;
		  z[ip] -= h ;
		  z[iq] += h ;
		  d[ip] -= h ;
		  d[iq] += h ;
		  a[ip][iq] = 0.0 ;
		  for (j=0 ; j<=ip-1 ; j++)
		    {
		      ROTATE(a,j,ip,j,iq)
		    }
		  for (j=ip+1 ; j<=iq-1 ; j++)
		    {
		      ROTATE(a,ip,j,j,iq)
		    }
		  for(j = iq+1 ; j<n ; j++)
		    {
		      ROTATE(a,ip,j,iq,j)
		    }
		  for (j=0 ; j<n ; j++)
		    {
		      ROTATE(v,j,ip,j,iq)
		    }
		  ++(*nrot);
		}
	    }
	}
      for(ip=0 ; ip<n ; ip++)
	{
          b[ip] += z[ip] ;
          d[ip]  = b[ip] ;
          z[ip]  = 0.0   ; 
	}
    }
  fprintf(stderr, "too many iterations in jacobi\n");
}


void 
eigsrt(double *d, double **v, int n)
{
  int    i, j, k ;
  double p       ;
  for (i=0 ; i < n-1 ; i++)
    {
      p = d[k=i] ;
      for (j=i+1 ; j<n ; j++)
	if (d[j] >= p) p = d[k=j] ;
      if (k != i)
	{
	  d[k] = d[i] ;
	  d[i] = p    ;
	  for (j=0 ; j<n ; j++)
	    {
	      p       = v[j][i] ;
	      v[j][i] = v[j][k] ;
	      v[j][k] = p       ;
	    }
	}
    }
}

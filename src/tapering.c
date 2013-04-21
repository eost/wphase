#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rwsacs.h"
 
void 
taper_syn(double *R, double *T, double *P, sachdr *hdr)
{
  int    n1, n2, j;
  double fac, blr, blt, blp;

  if (hdr->t[0]==-12345.0 || hdr->delta == -12345.0 || 
      hdr->b==-12345.0|| hdr->npts == -12345)
    {
      fprintf(stderr,"*** ERROR - P arrival time not specified in sac header T[0] ***\n");
      exit(1);
    }
  n1 = 0;
  n2 = (int)((hdr->t[0]-hdr->b)/hdr->delta);
  if (n2 > hdr->npts) n2 = hdr->npts;  
  if(n2 > n1+1)
    {
      blr = 0.;
      blt = 0.;
      blp = 0.;
      for (j=n1; j<n2; j++)  
	{
	  blr += R[j] ;
	  blt += T[j] ;
	  blp += P[j] ;        
	}
      blr /= (float)(n2-n1) ;
      blt /= (float)(n2-n1) ;
      blp /= (float)(n2-n1) ;
      for (j=0; j<hdr->npts; j++)
	{
	  R[j] -= blr ;
	  T[j] -= blt ;
	  P[j] -= blp ;        
	}
      
      for (j=n1; j<n2; j++)  
	{
	  fac   = (1.-cos((float)(j-n1)*M_PI/(float)(n2-n1-1)))/2.;
	  R[j] *= fac ;
	  T[j] *= fac ;
	  P[j] *= fac ;        
	}
    }
  return;
}

void 
taper_syn_only_Z(double *Z, sachdr *hdr)
{
  int    n1, n2,j;
  double fac, bl;
  
  if (hdr->t[0]==-12345.0 || hdr->delta == -12345.0 || 
      hdr->b==-12345.0|| hdr->npts == -12345)
    {
      fprintf(stderr,"*** ERROR - P arrival time not specified in sac header T[0] ***\n");
      exit(1);
    }

  n1 = 0;
  n2 = (int)((hdr->t[0]-hdr->b)/hdr->delta);
  if (n2 > hdr->npts) n2 = hdr->npts;    
  if(n2 > n1+1)
    {
      bl = 0.;
      for (j=n1; j<n2; j++)
	{
	  bl += Z[j];
	}
      bl /= (float)(n2-n1);
      for (j=0; j<hdr->npts; j++)
	{
	  Z[j] -= bl;
	}
      
      for (j=n1; j<n2; j++)
	{
	  fac   = (1.-cos((float)(j-n1)*M_PI/(float)(n2-n1-1)))/2.;
	  Z[j] *= fac;
	}
    }
  return;
}

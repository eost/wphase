#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c          */
#include "rwsacs.h"       /* rwsacs.c               */
#include "butterworth.h"  /* butterworth_sin.c */


void 
make_stf(int lh, char *type, double *h)
{
  int    i, LH;
  double al0 ;

  LH = lh+1 ;
  al0=((double) lh)/1.51743;
  if (type[0] == 'g')      /* gaussian */
    for (i = 0; i <= lh ; i++)
      h[i] = exp(- pow( ((double) i)/al0 ,2)) ;
  else if (type[0] == 'q') /* parabolic */
    for (i = 0; i <= lh ; i++)
      h[i] = 1. - pow( ((double) i)/((double)lh) ,2) ;
  else if (type[0] == 'l') /* triangle */
    for (i = 0; i <= lh ; i++)
      h[i] = 1. - ((double) i)/((double)lh) ;
  else if (type[0] == 'b') /* box car */
    for (i = 0; i <= lh ; i++)
      h[i] = 1. ;
  else if (type[0] == 'c') /* cosine */
    for (i = 0; i <= lh ; i++)
      h[i] = 1. + cos(M_PI*((double) i)/((double)lh));
  else
    {
      fprintf (stderr, "ERROR: Invalid STF type (arg 4)\n") ;
      fprintf (stderr, "       arg 4 must be one of these letters:\n\t g :gaussian \n\t q : parabolic \n\t l : triangle \n\t b : box car \n\t c : cosine\n") ;
      exit (1) ;
    }

  al0 = 0.;
  for ( i=0; i<LH; i++ )
    al0 = al0 + h[i];
  al0 = 2. * al0 - h[0] ;
  for (i = 0; i <= lh ; i++)
    h[i] = h[i]/al0;
}



void 
runave(double *x, double *y, int npts, int lh, char *type)
{
  int    LH, end, pos;
  int    i, j, il, ir;
  double cum, *h;

  
  h = double_alloc(1000) ; 

  /* Set STF */
  make_stf(lh,type,h);
  
  /* Ponderated Average in the inner section of the time series*/
  LH = lh + 1 ; /* half nb of samples (lh is the half nb of intervals)*/
  end = npts - LH ;
  pos = lh;
  for (i = lh ; i <= end ; i++)
    {
      il = i;
      ir = i;
      cum = h[0] * x[i] ; 
      for (j = 1; j < LH; j++)
	{
	  il = il - 1;
	  ir = ir + 1;
	  cum = cum + h[j] * ( x[il] + x[ir] ) ;
	}
      y[pos] = cum ;
      pos++;
    }
  /* Set average values for borders */
  for( i=0 ; i < lh ; i++ )
    {
      y[i] = y[lh] ;
      y[npts - lh + i] = y[npts - LH];
    }
  free((void*)h);
}



void 
conv_by_stf(double *delay, double *half_width, char *itype, sachdr *hdr, double *x_in, double *x_conv)
{
  int lh;

  /* Remove baseline */
  dtrd (x_in, hdr->npts, 60) ;
  printf("baseline removed\n");
  hdr->b = hdr->b + (*delay);

  /* Perform convolution */
  lh = (int) ((*half_width)/hdr->delta + 0.1) ;
  printf("b, nlh, dt, npts=%f %d %f %d\n", hdr->b, lh, hdr->delta, hdr->npts);
  runave(x_in,x_conv,hdr->npts,lh,itype) ;

  /* Remove baseline again */
  dtrd (x_conv, hdr->npts, 60) ;
}

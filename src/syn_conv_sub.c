/***************************************************************************
*
*	              W phase source inversion package 	            
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

/*	STF convolution        */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c          */
#include "rwsacs.h"       /* rwsacs.c               */


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
  int i, j, il, ir, LH;
  double cum, *h, xl, xr;

  /* STF */
  LH = lh + 1 ; /* #samples/2 (lh is #intervals/2) */
  h = double_alloc(LH) ;
  make_stf(lh,type,h);  
  /* Weighted Average */
  for (i = 0 ; i < npts ; i++)
    {
      il = i;
      ir = i;
      cum = h[0] * x[i] ;
      for (j = 1; j < LH; j++)
	{
	  il--;
	  ir++;
	  xl = (il< 0   ) ? 0. : x[il];
	  xr = (ir>=npts) ? 0. : x[ir];
	  cum = cum + h[j] * ( xl + xr ) ;
	}
      y[i] = cum ;
    }
  free((void*)h);
}

/************************************************************************/
/*                      left_taper(x, npts, start, end)                      */
/************************************************************************/
/*    Tapering                                                          */
void left_taper(double *x, int npts, int nbeg, int nend)
{
  int i;
  double ang, av = 0.;
  if ( nend > nbeg)
	{
	  for (i = nbeg; i <= nend ; i++)
		av += x[i];
	  av /= (nend-nbeg+1);
	  
	  for (i = 0; i < nbeg; i++)
		x[i] = 0.;
	  for (i = nbeg; i < npts ; i++)
		x[i] -= av;
	  
	  ang = M_PI / ((double)(nend-nbeg)) ;
	  for (i = nbeg; i <= nend ; i++)
		x[i] *= (1.0 - cos((double)(i-nbeg)*ang))/2.0 ;
	}
  else
	  printf("Warning 'nend <= nbeg'; not applying taper\n");
}

void
conv_by_stf(double delay,double half_width,char *itype,sachdr *hdr,double *x_in,double *x_conv)
{
  int lh, nbeg, nend;
  const int max_base     = 59;   /*  SAMPLES  */
  const int pre_P_safety = 10;   /*  SAMPLES  */

  nbeg = 0;
  nend =  (int)floor((hdr->t[0] - hdr->b - pre_P_safety)/hdr->delta);
  if ( nend - nbeg > max_base ) nbeg = nend - max_base;

  left_taper(x_in, hdr->npts, nbeg, nend);
  lh = (int)floor(half_width/hdr->delta + 0.5);
  runave(x_in,x_conv,hdr->npts,lh,itype);
  hdr->b    += delay;
}

/****************************************************************
*	W phase package - Butterworth filtering
*                                           
*       History
*             2008  Original Coding     Luis Rivera and Hiroo Kanamori
*             2010  Further Updates     Zacharie Duputel
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

/********************************************************************************/
/* ROUTINES IN C FOR FILTERING DATA USING A DIGITAL BANDPASS BUTTERWORTH FILTER */
/********************************************************************************/
/* Main routines :                                                              */
/* butter_sos__ permits to design a second order section of the filter          */
/* filter_data__ filters data from a second order section                       */
/********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "proto_alloc.h"    /* proto_alloc.c  */
#include "rwtextfiles.h"    /* rwtextfiles.c  */
#include "butterworth.h"    


/****************************************************/
/*               m=calc_mean(x,n)                   */
/****************************************************/
/* > Calculate mean (m) of the n first samples of x */
double 
calc_mean(double *x,int *nbmoy)
{
  int    i ;
  double s1, avex ;
  s1 = 0.;
  for ( i=0; i<(*nbmoy); i++ )
    s1 = s1 + x[i];
  avex = s1 /( (double)(*nbmoy) ); /* mean */
  return avex;
}

/**************************************************/
/*             dtrd(x, npts, y, NDEG)             */
/**************************************************/
/*  > Removes trend and/or mean from signal       */
/*  input : x : input signal                      */
/*          npts : nb of data points              */
/*          NDEC : flag to specify which of mean  */
/*                 or trend as to be removed      */
/*                                                */
/*  output signal : x                             */
void 
dtrd(double *x, int npts, int nbmoy)
{
  double an, s1, s2, slope, avei, avex;
  int i;
  an = npts;
  s1 = 0.0 ;
  s2 = 0.0 ;
  if(nbmoy != npts)
    {
      if (nbmoy > npts)
	nbmoy = npts;
      avex = calc_mean(x,&nbmoy) ;
      for ( i=0; i<npts; i++ )
	x[i] = x[i] - avex;          /* rm mean */
      avex = calc_mean(x,&nbmoy) ;
    }
  else
    {
      for ( i=0; i<npts; i++ )
	{
	  s1 = s1 + x[i] ;
	  s2 = s2 + s1 ;
	}
      avei = 0.5 * ((double)npts + 1.0);
      avex = s1 / ( (double)npts );
      slope = - 12.0 * (s2 - avei * s1) / (npts*(pow(npts,2) - 1.0)) ;
      for ( i=0; i<npts; i++ )
	x[i] = x[i] - avex - slope * ((double)i - avei) ; /* rm trend + mean */
    }
}

void 
taper(double *x, int npts, double start, double end)
{
  int i, m1, m2, m3;
  double ang, xi, cs;


  m1 = (int) ( npts * start + 0.5 ) ; /* left side -> nb of samples to be tapered */ 
  
  if (m1 > 0)
    {
      ang = M_PI / ( (double)m1 ) ;
      for (i=0; i < m1 ; i++)
	{
	  xi = (double)(i + 1) ;
	  cs = (1.0 - cos(xi * ang))/2.0 ;
	  x[i] = x[i] * cs ; /* tapering x by a sin2 */
	} 
    }
  m2 = (int)(npts * ((int)end) + 0.5) ; /* right side -> nb of samples to be tapered */ 
  if (m2 > 0)
    {
      m3 = npts - m2 ;
      ang = M_PI / ( (double)m2) ;
      for (i=m3 ; i < npts ; i++)
	{
	  xi = (double)(i - npts) ;
	  cs = ( 1.0 - cos(xi * ang)/2 ) ;
	  x[i] = x[i] * cs ; /* tapering x by a sin2 */
	}
    }
}






/*************************************************************************/
/*                    out=decrec (b1, b2, a1, k, fi, fo)                 */
/*************************************************************************/
/*    Recusive function to apply a recursive filtering in the special    */
/*    case where fo(1)=fo(2)=0 and a2 =0.                                */
/*    >denomonator polynomial is z**2 + a1*z       (!!a2 == 0!!)         */
/*    >numerator polynomial is z**2 + b1*z + b2                          */
/*         fi  = input array                                             */
/*         fo  = output array                                            */   
/*	   out = fo[k]                                                   */
/*    Remark : this function appear to be less effective than filt       */
double 
decrec (double b1, double b2, double a1, int k, double *fi, double *fo)
{
  if ( k == 1 )
    {
      fo[0] = 0 ;
      fo[1] = 0 ;
    }
  else
    {
      fo[k] = fi[k] + b1*fi[k-1] + b2*fi[k-2] ; 
      fo[k] = fo[k] - a1*decrec(b1, b2, a1, k-1, fi, fo) ;
    }
  return fo[k];
}


/*******************************************************************/
/*        out=filtrec (b1, b2, a1, a2, k, fi, fo)                  */
/*******************************************************************/
/*    Recusive function to apply a second order recursive          */
/*    filter to the data                                           */
/*    >denomonator polynomial is z**2 + a1*z + a2                  */
/*    >numerator polynomial is z**2 + b1*z + b2                    */
/*         fi  = input array                                       */
/*         fo  = output array                                      */   
/*	   out = fo[k]                                             */
/*    Remark : this function appear to be less effective than filt */
double 
filtrec (double b1, double b2, double a1, double a2, int k, double *fi, double *fo)
{
  if ( k==0 )
    fo[k] = fi[k] ;
  else if ( k == 1 )
    {
      fo[k] = fi[k] + b1*fi[k-1] ;
      fo[k] = fo[k] - a1*filtrec(b1, b2, a1, a2, k-1, fi, fo) ;
    }
  else
    {
      fo[k] = fi[k] + b1*fi[k-1] + b2*fi[k-2] ; 
      fo[k] = fo[k] - a1*filtrec(b1, b2, a1, a2, k-1, fi, fo) - a2*fo[k-2] ;
    }
  return fo[k];
}


/**************************************************************************/
/*                   bandpass (fl,fh,dt,n,p,b)                        */
/**************************************************************************/
/*
   Routine to compute bandpass filter poles for Butterworth filter
       fc = desired cutoff frequency
       dt = sample rate in seconds
       n  = number of poles (odd or even number is accepted)
       p  = pole locations  (RETURNED)
       z  = zero locations  (RETURNED)
       b  = gain factor     (RETURNED)
*/
/*   Program calculates a continuous Butterworth low pass IIRs with required */
/*    cut off frequency.                                                     */
/*   Then a discrete filter is calculated utilizing the bilinear transform   */
void 
bandpass    (double fl,  double fh,  double dt, int n, complex *p, double *b)
{
   double      wl, wh, wc, wc2, wb, b0 ;
   int         i, i1, inc, iinf, isup, ierror=1;
   complex     add_c(), mul_c(), div_c(), cmul_c(), sub_c() ;
   complex     conj_c() ;
   complex     one, x, y, sq_discr, tmp, pole, zc ;

/*             Frequency pre-warping, initialize variables         */
   wl  = tan(fl * dt * M_PI) ;
   wh  = tan(fh * dt * M_PI) ;
   wb  = wh - wl;
   wc2 = wh*wl;
   one = cmplx(1., 0.) ;
   tmp = cmplx(4.*wc2, 0.);

   if ( !(n%2) ) /* odd case   */
     {
       inc  = 2 ;
       iinf = 0 ;
       isup = n ;
     }
   else          /* even case */
     {
       inc  = 1 ;
       iinf = 1 ;
       isup = (n+1)/2 ;
       pole = cmplx(-wb,0.) ; /* pole at -1 for the lp prototype*/
       sq_discr  =  sqrt_c(sub_c(mul_c(pole, pole), tmp )) ; 
       p[0]      =  cmul_c(.5, add_c(pole, sq_discr)) ; /*bp*/
       p[1]      =  cmul_c(.5, add_c(pole, sq_discr)) ;      
       for( i=0 ; i< 2; i++)
	 {
	   x           = add_c(one,p[i]) ; /* analog -> digital using the */
	   y           = sub_c(one,p[i]) ; /*   bilinear transformation   */
	   p[i]        = div_c(x,y,&ierror) ;
	 }
     }

   for (i=iinf ; i < isup ; i += inc)
	 {
/*             Calculate position of poles for continuous filter   */
	   i1 = i + inc - 1 ;
	   pole.real = -wb*cos((float)i1*M_PI/(float)(inc*n)) ;
	   pole.imag =  wb*sin((float)i1*M_PI/(float)(inc*n)) ;   /* lp-prototype */
	   sq_discr  =  sqrt_c(sub_c(mul_c(pole, pole), tmp )) ; 
	   p[2*i]    =  cmul_c(.5, add_c(pole, sq_discr)) ;     /* bp */
	   p[2*i+2]  =  cmul_c(.5, sub_c(pole, sq_discr)) ;     /* bp */
	 }
   
   for ( i=2*iinf ; i<2*n ; i += 2)
   {
/*             Calculate position of poles for discrete    */
/*             filter using the bilinear transformation    */
     x           = add_c(one,p[i]) ;
     y           = sub_c(one,p[i]) ;
     p[i]        = div_c(x,y,&ierror) ;
     p[i+1]      = conj_c(p[i]) ;
   }

/*             Calculate bp filter gain   */
   wc           = 2.*M_PI*sqrt(fl*fh) ;
   zc           = cmplx(cos(wc*dt),sin(wc*dt));
   tmp          = sub_c(mul_c(zc,zc), one) ;
   b0           = sqrt(tmp.real*tmp.real+tmp.imag*tmp.imag) ;
   b0           = pow(b0, (double)n) ;
   for (i=0 ; i<2*n ; i +=2)
   {
       x        = sub_c(zc,p[i]) ;
       y        = sub_c(zc,p[i+1]) ;
       x        = mul_c(x,y) ;
       b0      /= sqrt(x.real*x.real+x.imag*x.imag) ;
   }
   b0 = 1./b0 ;
   *b = b0 ;
}


/**************************************************************************/
/*                    poles2sos(pb, n, b1, b2, a1, a2)                    */
/**************************************************************************/
/*
   Conversion from poles (in the array pb) to the second order section
   coefficients in the arrays b1, b2, a1 and a2 for Butterworth filter
   of order n.                                                             */
void 
poles2sos(complex *pb, int order, double *b1, double *b2, double *a1, double *a2)
{
  int     i ;
  for ( i=0 ; i<2*order ; i +=2)
    {
        b1[i/2] =  0. ;
	b2[i/2] = -1. ;
        a1[i/2] = -2.*pb[i].real ;
	a2[i/2] =  pb[i].real*pb[i].real + pb[i].imag*pb[i].imag ;
	
    }
}



/**************************************************************************/
/*           butter_sos__(fl, fh, n, sr, b1, b2, a1, a2, scale)           */
/**************************************************************************/
/*
  Creation of a second ordinary section for a digital butterworth filter
  
  input : fl : low cut off frequency
          fh : high cut off frequency
          n  : order of the filter (even or odd)
	  sr : sampling period (in seconds)
         
  output : b1, b2, a1, a2 pointer to an proto_allocd array containing the
           second order section 
	   scale : scale factor of the second order section               */
void 
butter_sos(double *lowCutoff, double *highCutoff, int *filterOrder, double *sampPeriod, double **b1, double **b2, double **a1, double **a2, double *gain) 
{
  complex *pb ;
  double   dt;
  pb = complex_alloc((*filterOrder)*2);

  /* Get low pass filter parameters    */
  dt = *sampPeriod;
 
  /* Get band pass filter poles if necessary     */
  bandpass(*lowCutoff,*highCutoff,dt,*filterOrder,pb,gain) ;

  /* Get second ordinary section from pole & zeros */
  (*b1) = double_alloc( *filterOrder ) ;
  (*b2) = double_alloc( *filterOrder ) ;
  (*a1) = double_alloc( *filterOrder ) ;
  (*a2) = double_alloc( *filterOrder ) ;
  poles2sos ( pb, *filterOrder, *b1, *b2, *a1, *a2 ) ;
  free((void *)pb);
  return;
}


/*************************************************************/
/*    filter_data__ (di, np, n, b1, b2, a1, a2, cst)         */
/*************************************************************/
/*     Routine to apply a second order section to the data
       denomonator polynomial is z**2 + a1*z + a2
       numerator polynomial is z**2 + b1*z + b2
           di  = input array
           do  = output array
           np  = number of points
	   cst = scale factor of the sos
	   n   = number of rows of the matrix sos
	   opt = optionnal parameter for the deconvolution
                  > if opt == 1 the two first samples of the 
		                deconvolved signal are zero.
		  > if opt == 0 the convolution is computed
		                as the remaining sos. 
*/
void 
filter_data(double *data, int *numSamp, int *nfilt, double *b1, double *b2, double *a1, double *a2, double *scale, int *opt)
{
  int    i ;
  double  tmp;

  if ((*opt)<0 || (*opt)>1)
    {
      fprintf(stderr,"ERROR : invalid parameter opt given to filter_data\n");
      exit(1);
    }
  
  if (*opt)
    {
      for (i=0; i<3;i++)
	data[i] = 0.;
      tmp = decrec (b1[0], b2[0], a1[0], (*numSamp)-1, data, data) ;
    }

  for ( i=(*opt) ; i < (*nfilt) ; i ++)
    tmp=filtrec (b1[i], b2[i], a1[i], a2[i], (*numSamp)-1, data, data) ;

  for ( i=0 ; i<*numSamp ; i++) 
    data[i] =(*scale)*data[i] ;

  return;
}



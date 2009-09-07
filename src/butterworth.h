/* Header of butterworth_sin.c subroutines */

#include "complex.h"      /* complex.c      */


/****************************************************/
/*               m=calc_mean(x,n)                   */
/****************************************************/
/* > Calculate mean (m) of the n first samples of x */
double calc_mean(double *x,int *nbmoy);

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
void dtrd(double *x, int npts, int nbmoy);


void taper(double *x, int npts, double start, double end);



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
double decrec (double b1, double b2, double a1, int k, double *fi, double *fo);


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
double filtrec (double b1, double b2, double a1, double a2, int k, double *fi, double *fo);


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
void bandpass (double fl,  double fh,  double dt, int n, complex *p, double *b);


void  poles2sos(complex *pb, int order, double *b1, double *b2, double *a1, double *a2);


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
                                                                          
void butter_sos(double *lowCutoff, double *highCutoff, int *filterOrder, double *sampPeriod, double **b1, double **b2, double **a1, double **a2, double *gain);


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
void  filter_data(double *data, int *numSamp, int *nfilt, double *b1, double *b2, double *a1, double *a2, double *scale, int *opt);

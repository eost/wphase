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

/*      Butterworth filtering                                   */
/* WARNING: Everything is coded assuming even filter order      */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "proto_alloc.h"


/**************************************************/
/*             dtrd(x, npts, y, NDEG)             */
/**************************************************/
/* Removes trend and/or mean from signal          */
void dtrd(double *x, int npts)
{
    double s1, s2, slope, avei, avex;
    int i;
    s1 = 0.0;
    s2 = 0.0;  
    for ( i=0; i<npts; i++ )
    {
        s1 = s1 + x[i];
        s2 = s2 + s1;
    }
    avei = 0.5 * ((double)npts + 1.0);
    avex = s1 / ( (double)npts );
    slope = - 12.0 * (s2 - avei * s1) / (npts*(pow(npts,2) - 1.0));
    for ( i=0; i<npts; i++ )
        x[i] = x[i] - avex - slope * ((double)i - avei); /* rm trend + mean */
}

/****************************************************/
/*             rmean(x, npts, jbeg, jend)           */
/****************************************************/
/*  Substracts a constant from the input vector     */
/*  The constant is defined as the average          */
/*  of a segment of the input signal.               */
void rmean (double *x, const int npts, int *jbeg, int *jend)
{
    double ave;
    int    j;

    if ( *jbeg < 0 ) 
        *jbeg = 0;
    if ( *jend  > npts-1 )  
        *jend = npts-1;

    ave = 0.0;
    for (j=*jbeg; j<=*jend; j++)
        ave += x[j];
    ave /= (double)(*jend - *jbeg + 1);

    for (j=0; j<npts; j++)
        x[j] -= ave;
}

/************************************************************************/
/*                      taper(x, npts, start, end)                      */
/************************************************************************/
/*    Tapering                                                          */
void taper(double *x, int npts, int nbeg, int nend)
{
    int i, iend;
    double ang, xi, cs;
    ang = M_PI / ((double)nbeg);
    for (i=0; i<nbeg ; i++)
    {
        xi = (double)(i + 1);
        cs = (1.0 - cos(xi*ang))/2.0;
        x[i] = x[i] * cs; /* taper with a sin2 */
    }   
    iend = npts-nend; 
    ang = M_PI / ( (double)nend);
    for (i=iend ; i<npts ; i++)
    {
        xi = (double)(i - npts);
        cs = ( 1.0 - cos(xi * ang)/2 );
        x[i] = x[i] * cs; /* taper with a sin2 */
    }
}

/************************************************************************/
/*               app_sos( b1, b2, a1, a2, npts, sig)                    */
/************************************************************************/
/* Apply a second order section                                         */
void app_sos(double b1, double b2, double a1, double a2, int npts, double *sig)
{
    int    k;
    double x0, x1, x2, y0, y1, y2;
    x1 = 0., x2 = 0.,  y1 = 0., y2 = 0.;
    for (k=0; k<npts; k++)
    {
        x0 = sig[k];
        y0 = x0 + b1*x1 + b2*x2 - a1*y1 - a2*y2;
        x2 = x1; x1 = x0;
        y2 = y1; y1 = y0;
        sig[k] = y0;
    }
}

/**************************************************************************/
/*                   lpbu2sos (fc,dt,n,b1,b2,a1,a2,g)                     */
/**************************************************************************/
/* Compute lowpass filter poles for Butterworth filter                    */
int lpbu2sos(double fc, double dt, int n, double *g, double *b1, double *b2, double *a1, double *a2)
{
    double   ang, wc, gain;
    int      j, k;
    complex  pole, cone, swp, *p;

    p = (complex*)calloc(n, sizeof(complex));
    cone = 1.;
    wc   = tan(fc * dt * M_PI);
    for (j=0; j<n; j++)
    {
        ang  = (double)(2*j+1)*M_PI/(double)(2*n);
        pole =  wc*(-sin(ang) + cos(ang)*_Complex_I);
        p[j] = (cone+pole)/(cone-pole);
    }

    for(j=0; j<n; j++)
        for(k=0; k<n; k++)
            if (carg(p[j]) < carg(p[k]))
            {
                swp  = p[j];
                p[j] = p[k];
                p[k] = swp;
            }
    for(j=0; j<n/2; j++)
    {
        b1[j] =  2.;
        b2[j] =  1.;
        a1[j] = -2.*creal(p[j]);
        a2[j] =  creal(p[j])*creal(p[j]) + cimag(p[j])*cimag(p[j]);
    }

    /*  Calculate lp filter gain   */
    gain = 1. ;
    for (k=0 ; k<n ; k ++)
        gain  /= cabs(cone-p[k])/2. ;
    *g = 1./gain ; 

    free(p);
    return n/2;
}

/**************************************************************************/
/*                   bpbu2sos (fl,fh,dt,n,g,b1,b2,a1,a2)                    */
/**************************************************************************/
/* Compute bandpass filter sos for Butterworth filter                   */
int bpbu2sos(double fl, double fh, double dt, int n, double *g, double *b1, double *b2, double *a1, double *a2)
{
    double   wl, wh, wc, wc2, wb, ang, gain;
    int      j, k;
    complex  sq_discr, pole, ctwo, swp, zc, *p;
  
    p = (complex*)calloc(2*n, sizeof(complex));

    ctwo = 2.;
    wl   = tan(fl * dt * M_PI);
    wh   = tan(fh * dt * M_PI);
    wb   = wh - wl;
    wc2  = wh*wl;
  
    for (j=0; j<n; j++)
    {
        ang      = (double)(2*j+1)*M_PI/(double)(2*n);
        pole     =  wb*(-sin(ang) + cos(ang)*_Complex_I);
        sq_discr = csqrt(pole*pole-4.*wc2) ;
        p[j]     = (ctwo+pole+sq_discr)/(ctwo-pole-sq_discr);
        p[n+j]   = (ctwo+pole-sq_discr)/(ctwo-pole+sq_discr);
    }
 
    for(j=0; j<2*n; j++)
        for(k=0; k<2*n; k++)
            if (carg(p[j]) < carg(p[k]))
            {
                swp  = p[j];
                p[j] = p[k];
                p[k] = swp;
            }

    for(j=0; j<n; j++)
    {
        b1[j] =  0.;
        b2[j] = -1.;
        a1[j] = -2.*creal(p[j]);
        a2[j] =  creal(p[j])*creal(p[j]) + cimag(p[j])*cimag(p[j]);
    }

    /*  Calculate bp filter gain   */
    wc = 2.*M_PI*sqrt(fl*fh);
    zc = cos(wc*dt)+sin(wc*dt)*_Complex_I;
    gain = pow(cabs(zc*zc - 1.0),n);
    for (k=0 ; k<2*n ; k ++)
        gain  /= cabs(zc-p[k]);
    *g = 1./gain;

    free(p);
    return n;
}


/****************************************************/
/*   filter_with_sos (cst,b1,b2,a1,a2,n,sig,npts)   */
/****************************************************/
/* Apply a second order section to sig              */
void filter_with_sos(double g, double *b1, double *b2, double *a1, double *a2, int nsects, double *sig, int npts)
{
    int j;
    for (j=0; j<nsects; j++)
        app_sos(b1[j], b2[j], a1[j], a2[j], npts, sig);
    for (j=0; j<npts; j++)
        sig[j] *= g;
}


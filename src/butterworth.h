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

/*  W phase package - Butterworth filtering headers            */
/* WARNING: Everything is coded assuming even filter order     */

/****************************************************/
/*                  dtrd(x, npts)                   */
/****************************************************/
/*  > Removes trend and/or mean from signal         */
/*  input : x : input signal                        */
/*          npts : nb of data points                */
/*  output signal : x                               */
void dtrd(double *x, int npts);

/****************************************************/
/*             rmean(x, npts, jbeg, jend)           */
/****************************************************/
/*  Substracts a constant from the input vector     */
/*  The constant is defined as the average          */
/*  of a segment of the input signal.               */
/*                                                  */
/*  input : x        : input signal                 */
/*          npts     : length of x                  */
/*          jbeg,jend: defines the segment of x     */
/*                     whose average is to be       */
/*                     substracted from x.          */
/*                     Both x[jbeg] and x[jend]     */
/*                     are included in the segment  */
/*  output signal    : x                            */
/*  Note:                                           */
/*  if jbeg < 0 or jend > npts-1 they are modified  */
/****************************************************/
void rmean (double *x, const int npts, int *jbeg, int *jend);

/****************************************************/
/*               taper(x, npts, p1, p2)             */
/****************************************************/
/* Taper of x from samples 0 to nbeg (left side)    */
/*        and from npts-nend to npts-1 (right side) */
void taper(double *x, int npts, int nbeg, int nend);

/**************************************************************************/
/*                   lpbu2sos (fc,dt,n,b1,b2,a1,a2,g)                     */
/**************************************************************************/
/* Compute lowpass filter poles for Butterworth filter                    */
int lpbu2sos(double fc, double dt, int n, double *g, double *b1, double *b2, double *a1, double *a2);

/**************************************************************************/
/*                   bpbu2sos (fl,fh,dt,n,g,b1,b2,a1,a2)                    */
/**************************************************************************/
/* Compute bandpass filter sos for Butterworth filter                   */
int bpbu2sos(double fl, double fh, double dt, int n, double *g, double *b1, double *b2, double *a1, double *a2);

/****************************************************/
/*   filter_with_sos (cst,b1,b2,a1,a2,n,sig,npts)   */
/****************************************************/
/* Apply a second order section to the signal sig   */
/*     denomonator polynomial is z**2 + a1*z + a2   */
/*     numerator polynomial is z**2 + b1*z + b2     */
/*     cst         : scale factor of the sos        */
/*     b1,b2,a1,a2 : sos                            */
/*     n           : number of rows in the sos      */
/*     sig         : input/output signal            */
/*      npts        : number of sample in sig       */
void filter_with_sos(double g, double *b1, double *b2, double *a1, double *a2, int nsects, 
		    double *sig, int npts);



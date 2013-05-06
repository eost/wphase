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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rwsacs.h"
#include "proto_alloc.h"


/* Internal subroutine */
void get_name(char *A, char *B, char *C, char *D, char *name);

/* Compute vertical displacement from GF and moment components */
void 
summ_up(double *M2_cmt, char *segm1, char *segm2, double **GFs, 
	double *R, double *T, double *P, sachdr *hdr)
{
  int  j, k, kmin, max, ierror=1;
  char *sac_GF;

  max = hdr->npts ;
  sac_GF = char_alloc(FSIZE);

  /* Vertical     */
  get_name(segm1, "/RR/", segm2, "Z.SAC", sac_GF);
  rhdrsac(sac_GF, hdr, &ierror) ; /* Read header */
  if (hdr->npts != max) hdr->npts = max ;
  rdatsac(sac_GF, hdr, GFs[0], &ierror) ; /* Read samples */

  get_name(segm1, "/TT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[1], &ierror) ; /* Read samples */

  get_name(segm1, "/PP/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[2], &ierror) ; /* Read samples */

  get_name(segm1, "/RT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[3], &ierror) ; /* Read samples */


  /* Longitudinal */
  get_name(segm1, "/RR/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[4], &ierror) ; /* Read samples */

  get_name(segm1, "/TT/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[5], &ierror) ; /* Read samples */

  get_name(segm1, "/PP/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[6], &ierror) ; /* Read samples */

  get_name(segm1, "/RT/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[7], &ierror) ; /* Read samples */

  /* Transversal  */
  get_name(segm1, "/RP/", segm2, "T.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[8], &ierror) ; /* Read samples */

  get_name(segm1, "/TP/", segm2, "T.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[9], &ierror) ; /* Read samples */

  kmin=4;
  for(k=0; k<4; k++)
    if(fabs(M2_cmt[k]) > 0.)
      {
	kmin = k;
	break; 
      }
  
  for(k=kmin; k<4; k++)
    if(fabs(M2_cmt[k]) > 0.)
      for(j=0; j < hdr->npts; j++) {
	R[j] += ((double)M2_cmt[k])*GFs[  k][j];
	T[j] += ((double)M2_cmt[k])*GFs[4+k][j]; }
  
  kmin=6;
  for(k=4; k<6; k++)
    if(fabs(M2_cmt[k]) > 0.) {
      kmin = k;
      break; }
  
  for(k=kmin; k<6; k++)
    if(fabs(M2_cmt[k]) > 0.) 
      for(j=0; j < hdr->npts; j++)
	P[j] += M2_cmt[k]*GFs[4+k][j];


  free((void *)sac_GF);
  return;
}

/* Compute vertical displacement from GF and moment components */
void 
summ_up_only_Z(double *M2_cmt, char *segm1, char *segm2, double **GFs, double *R, sachdr *hdr)
{
  int  j, k, kmin, max, ierror=1;
  char *sac_GF;

  max = hdr->npts ;
  sac_GF = char_alloc(FSIZE);

  /* Read GF for Vertical displacement */
  get_name(segm1, "/RR/", segm2, "Z.SAC", sac_GF);
  rhdrsac(sac_GF, hdr, &ierror) ;
  if (hdr->npts != max) hdr->npts = max ;
  rdatsac(sac_GF, hdr, GFs[0], &ierror) ;

  get_name(segm1, "/TT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[1], &ierror) ;

  get_name(segm1, "/PP/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[2], &ierror) ;

  get_name(segm1, "/RT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdr, GFs[3], &ierror) ;

  kmin=4;
  for(k=0; k<4; k++)
    if(fabs(M2_cmt[k]) > 0.) 
      {
	kmin = k;
	break;
      }
  
  for(k=kmin; k<4; k++)
    if(fabs(M2_cmt[k]) > 0.)
      for(j=0; j < hdr->npts; j++)
	R[j] += ((double)M2_cmt[k])*GFs[  k][j];

  free((void *)sac_GF);
  return;
}


void 
get_name(char *A, char *B, char *C, char *D, char *name)
{
    strncpy(name, A, FSIZE);
    strcat( name, B);
    strcat( name, C);
    strcat( name, D);
    return;
}

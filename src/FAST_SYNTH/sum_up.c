#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../proto_alloc.h"
#include "../complex.h"
#include "../rwsacs.h"

/* Internal subroutine */
void get_name(char *A, char *B, char *C, char *D, char *name);

/* Compute vertical displacement from GF and moment components */
void 
summ_up(M2_cmt, segm1, segm2, GFs, R, T, P, hdr, hdrg)
     double *M2_cmt, ** GFs ;
     double *R, *T, *P      ;
     char   *segm1, *segm2   ;
     sachdr *hdr, *hdrg     ;
{
  int  j, k, kmin, max, ierror=1;
  char *sac_GF;

  max = hdr->npts ;
  
  sac_GF = char_alloc(FSIZE);

  /* Vertical     */
  get_name(segm1, "/RR/", segm2, "Z.SAC", sac_GF);
  rhdrsac(sac_GF, hdrg, &ierror) ; /* Read header */
  if (hdrg->npts != max) hdrg->npts = max ;
  rdatsac(sac_GF, hdrg, GFs[0], &ierror) ; /* Read samples */

  get_name(segm1, "/TT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[1], &ierror) ; /* Read samples */

  get_name(segm1, "/PP/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[2], &ierror) ; /* Read samples */

  get_name(segm1, "/RT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[3], &ierror) ; /* Read samples */


  /* Longitudinal */
  get_name(segm1, "/RR/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[4], &ierror) ; /* Read samples */

  get_name(segm1, "/TT/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[5], &ierror) ; /* Read samples */

  get_name(segm1, "/PP/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[6], &ierror) ; /* Read samples */

  get_name(segm1, "/RT/", segm2, "L.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[7], &ierror) ; /* Read samples */

  /* Transversal  */
  get_name(segm1, "/RP/", segm2, "T.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[8], &ierror) ; /* Read samples */

  get_name(segm1, "/TP/", segm2, "T.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[9], &ierror) ; /* Read samples */

  /* Set header variables */
  hdr->delta = hdrg->delta ;
  hdr->b     = hdrg->b     ;

  kmin=4;
  for(k=0; k<4; k++)
    if(fabs(M2_cmt[k]) > 0.){
      kmin = k;
      break; }
  
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

void 
get_name(char *A, char *B, char *C, char *D, char *name)
{
    strncpy(name, A, FSIZE);
    strcat( name, B);
    strcat( name, C);
    strcat( name, D);
    return;
}

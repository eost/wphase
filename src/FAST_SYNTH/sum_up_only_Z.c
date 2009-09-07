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
summ_up_only_Z(M2_cmt, segm1, segm2, GFs, R, hdr, hdrg)
     double *M2_cmt, ** GFs, *R ;
     char   *segm1, *segm2      ;
     sachdr *hdr, *hdrg         ;
{
  int  j, k, kmin, max, ierror=1;
  char *sac_GF;

  max = hdr->npts ;
  
  sac_GF = char_alloc(FSIZE);

  /* Read GF for Vertical displacement */
  get_name(segm1, "/RR/", segm2, "Z.SAC", sac_GF);
  rhdrsac(sac_GF, hdrg, &ierror) ;
  if (hdrg->npts != max) hdrg->npts = max ;
  rdatsac(sac_GF, hdrg, GFs[0], &ierror) ;

  get_name(segm1, "/TT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[1], &ierror) ;

  get_name(segm1, "/PP/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[2], &ierror) ;

  get_name(segm1, "/RT/", segm2, "Z.SAC", sac_GF);
  rdatsac(sac_GF, hdrg, GFs[3], &ierror) ;

  /* Set header variables */
  hdr->delta = hdrg->delta;
  hdr->b = hdrg->b;

  kmin=4;
  for(k=0; k<4; k++)
    {
      if(fabs(M2_cmt[k]) > 0.) 
	{
	  kmin = k;
	  break;
	}
    }
  
  for(k=kmin; k<4; k++)
    {
      if(fabs(M2_cmt[k]) > 0.)
	{
	  for(j=0; j < hdr->npts; j++)
	    {
	      R[j] += ((double)M2_cmt[k])*GFs[  k][j];
	    }
	}
    }

  free((void *)sac_GF);
  return;
}

void 
get_name(char *A, char *B, char *C, char *D, char *name)
{
    strncpy(name, A, 128);
	strcat( name, B);
	strcat( name, C);
	strcat( name, D);
	return;
}

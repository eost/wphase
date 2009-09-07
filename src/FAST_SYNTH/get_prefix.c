#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../proto_alloc.h"
#define NDEPTHS 200

/* External subroutine */
void get_depths(char *path, float *depths, int *nd);

/******************************************************************/
/*     get_prefix(dep, xdeg, bsegm1, bsegm2, bdepth, bdist )      */
/******************************************************************/
/* Select the best depth and dist for Green's function extraction */
/* Input : dep  : focal depths (from pde)                         */
/*         xdeg : epicentral distance                             */
/* Ouput : bsegm1 : path to best depths GF                        */
/*         bsegm2 : sac file containing the best dist GF          */
/*         bdepths: best focal depths                             */
/*         bdist  : best focal distance                           */
void 
get_prefix(float cmt_dep, float xdeg, char *best_segm1, char *best_segm2, 
				float *best_depth, float *best_dist )
{
  char *gf_path;
  float depths[NDEPTHS];
  int   jd, nd, idist;

/* ############################################################### */
/* This section should be updated if the distance grid is modified */
  idist       = 2*(int)floor(5.*xdeg)+1;
  if (idist <   1) 
    {
      fprintf(stderr, "Warning distance too short!; using 0.1 deg\n");
      idist = 1;
    }
  if (idist > 899) 
    {
      fprintf(stderr, "Warning distance too big!;  using 89.9 deg\n");
      idist = 899;
    } 
/* ############################################################### */

  gf_path = char_alloc(FSIZE);
  
  if(getenv("GF_PATH"))
    strncpy(gf_path, getenv("GF_PATH"), FSIZE);
  else
    {
      fprintf(stderr, "Error: GF_PATH environment variable not defined\n");
      exit(1);
    }
  
  get_depths(gf_path, depths, &nd);

  /* Finding the best depth */
  *best_depth = depths[0];
  for (jd=1; jd < nd; jd++)
    {
      if(fabs(depths[jd]-cmt_dep) < fabs(*best_depth-cmt_dep)) 
	*best_depth = depths[jd];
      //printf("%6.1f %6.1f %6.1f\n", cmt_dep, depths[jd], *best_depth);
    }
  strncpy(best_segm1, gf_path, FSIZE);
  sprintf(best_segm2, "/H%05.1f",  *best_depth);
  strcat (best_segm1, best_segm2);
  sprintf(best_segm2, "GF.%04d.SY.LH", idist); 
  free((void *)gf_path);
  *best_dist = (float)idist/10.;
  //fprintf(stdout, "Using gf's %s/*/%sZ.SAC\n", best_segm1,best_segm2);
  return;
}

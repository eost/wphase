#include <stdio.h>
#include <string.h>
#include <math.h>
#include "rwsacs.h"
#include "proto_alloc.h"
#include "read_i_files.h"
#include "travel_times.h" 

int  yyyymmdd2jjj(int year, int month, int day) ;

void rotate_cmt(double* M_cmt, double* N_cmt, double az);

void get_prefix(double cmt_dep, double xdegs, char* segm1, char* segm2, 
		double* best_dep, double* best_dist);

long int get_length(char *segm1, char *segm2);

void summ_up(double* M2_cmt, char* segm1, char* segm2, double** GFs, 
	     double* R, double *T, double *P, sachdr *hdr) ;

void summ_up_only_Z(double *M2_cmt, char *segm1, char *segm2, double **GFs, 
		    double *R, sachdr *hdr) ;

/**************************************************************/
/* Computing synthetic seismograms for a given moment tensor  */
/**************************************************************/
/* Input params :                                             */
/*                - az   = azimuth (deg)                      */
/*                - xdeg = epicentral distance (deg)          */
/*                - eq.vm[1] = input moment tensor            */
/*                - eq.evdp = centroid depth                  */
/*                - GFs : Green function vectors              */
/*                - Z,TH,PH,E,N : vectors of __LEN_SIG__      */
/*                                                            */
/* Output params : - hdr  = synthetic SAC header information  */
/*                 - Z,TH,PH,E,N : synthetics                 */
void 
fast_synth_sub(double az, double baz, double xdeg, double *tv, double *dv, int nd, 
	       str_quake_params *eq, sachdr *hdr, double **GFs, double *Z, double *TH, double *PH)
{
  int i, ierror=1 ;
  char  segm1[FSIZE], segm2[FSIZE] ;
  double   best_depth, best_dist, P_tt ;  

  rotate_cmt(eq->vm[1], eq->vm[0], az) ;  /* Rotating the moment tensor   */
  get_prefix(eq->evdp, xdeg, segm1, segm2, &best_depth, &best_dist) ; /* Get the best GF filename     */
  hdr->npts  = get_length(segm1, segm2) ; /* GF sample number   */
  if (hdr->npts > __LEN_SIG__)            
    hdr->npts = __LEN_SIG__ ;                
  for(i=0; i<hdr->npts; i++) /* Set vectors to zero */
    {
      Z[i]  = 0. ;
      TH[i] = 0. ;
      PH[i] = 0. ;
    }
  summ_up(eq->vm[0], segm1, segm2, GFs, Z, TH, PH, hdr) ; /* Computing synthetics Z,TH,PH */
  trav_time(xdeg,tv,dv,nd,&P_tt,&ierror);
  hdr->o      = 0.   ; 
  hdr->az     = az   ;
  hdr->baz    = baz  ;
  hdr->gcarc  = xdeg ;
  hdr->t[0]   = P_tt ;
  strcpy(hdr->kt[0],"P       ");
  hdr->nzyear = eq->ot_ye ;
  hdr->nzhour = eq->ot_ho ;
  hdr->nzmin  = eq->ot_mi ;
  hdr->evla   = (float)eq->evla  ;
  hdr->evlo   = (float)eq->evlo  ;
  hdr->evdp   = (float)eq->evdp  ;
  hdr->nzsec  = floor(eq->ot_se) ;
  hdr->nzmsec = (int) (1000.*(eq->ot_se - hdr->nzsec) + .5) ;
  hdr->nzjday = yyyymmdd2jjj(hdr->nzyear, eq->ot_mo, eq->ot_dm);
}


void 
fast_synth_only_Z_sub(double az, double baz, double xdeg, double *tv, double *dv, int nd, 
		      str_quake_params *eq, sachdr *hdr, double **GFs, double *Z)
{
  int i ;
  char  segm1[FSIZE], segm2[FSIZE] ;
  double   best_depth, best_dist   ;  

  rotate_cmt(eq->vm[1], eq->vm[0], az) ;  /* Rotating the moment tensor   */
  get_prefix(eq->evdp, xdeg, segm1, segm2, &best_depth, &best_dist) ; /* Get the best GF filename     */
  hdr->npts  = get_length(segm1, segm2) ; /* GF sample number   */
  if (hdr->npts > __LEN_SIG__)            
    hdr->npts = __LEN_SIG__ ;                
  for(i=0; i<hdr->npts; i++) /* Set vectors to zero */
    Z[i]  = 0. ;
  summ_up_only_Z(eq->vm[0], segm1, segm2, GFs, Z, hdr) ; /* Computing synthetics Z */
  hdr->o      = 0.   ; 
  hdr->az     = az   ;
  hdr->baz    = baz  ;
  hdr->gcarc  = xdeg ;
  hdr->nzyear = eq->ot_ye ;
  hdr->nzhour = eq->ot_ho ;
  hdr->nzmin  = eq->ot_mi ;
  hdr->evla   = (float)eq->evla  ;
  hdr->evlo   = (float)eq->evlo  ;
  hdr->evdp   = (float)eq->evdp  ;
  hdr->nzsec  = floor(eq->ot_se) ;
  hdr->nzmsec = (int) (1000.*(eq->ot_se - hdr->nzsec) + .5) ;
  hdr->nzjday = yyyymmdd2jjj(hdr->nzyear, eq->ot_mo, eq->ot_dm);
}

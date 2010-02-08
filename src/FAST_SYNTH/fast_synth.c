#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../rwsacs.h"
#include "../proto_alloc.h"
#include "../travel_times.h"
#include "../read_i_files.h"

#define  YES  1
#define  NON  0


/* External routines */
void read_cmt(str_quake_params *eq, sachdr *hdr);

int  read_stats(char* file, char*** stats, char*** nets, float** stlats, 
		float** stlons);

void distaz(double cmt_lat, double cmt_lon, float* stlats, float* stlons, 
	    int nstat, float* dists, float* azs, float* bazs, float* xdegs,
	    long int* nerr);  

void rotate_cmt(double* M_cmt, double* N_cmt, double az);

void get_prefix(float cmt_dep, float xdegs, char* segm1, char* segm2, 
		float* best_dep, float* best_dist);

void summ_up(double* M2_cmt, char* segm1, char* segm2, double** GFs, 
	     double* R, double *T, double *P, sachdr* hdr, sachdr* hdrg);


void save_sac(char *stnm, char *netwk, char *chan, float *lat, float *lon, 
	      sachdr *hdr,  double*   depval);

long int get_length(char *segm1, char *segm2);

/**********************************************************************************************/
/* The variables and functions declared hereafter are used for tapering the non causal signal */
void init_tapering(int *nh, int *nd, float *h, double **dv, double **tv);
void taper(double *R, double *T, double *P, int *npts, float *dt, float *d, 
	   double *tv, double *dv, int *nd);
/**********************************************************************************************/


int 
main(int argc, char **argv)
{

  int j ;
  /* CMT variables */
  char   stat_file[FSIZE] ; 
  double *M1_cmt, *M2_cmt ; 
  str_quake_params eq     ;  
  /* Stations variables */
  char  **stats, **nets             ;
  float *stlats,*stlons             ;
  float *dists, *azs, *bazs, *xdegs ;
  int   jstat, nstat                ;
  /* Sac header variables */
  long int nerr                       ;
  char     segm1[FSIZE], segm2[FSIZE] ;
  float    best_depth, best_dist      ;  
  sachdr   hdr, hdrg                  ; 
  /* Seismograms */
  double **GFs ;
  double *Z, *TH, *PH ;
  /* Tapering */
  int    tapering = NON ;
  int    nh=NDEPTHS, nd=NDISTAS ;
  double *tv, *dv    ;

  /* Input parameters */
  if(argc != 3 && argc != 4)
    {
      fprintf(stderr,"Error: syntax: fast_synth cmt_file stats_file [-t]\n") ;
      exit(1) ;
    }

  /* Double allocations */
  GFs    = double_alloc2(10,__LEN_SIG__);/* GFs: Rrr, Rtt, Rpp, Rrt */
  M1_cmt = double_alloc(6)   ; /* Moment tensor           */  
  M2_cmt = double_alloc(6)   ; /* Rotated moment tensor   */
  eq.vm  = double_alloc2p(2) ;
  eq.vm[1] = &M1_cmt[0] ;  /* eq.vm[1] is extracted from cmtfile */
  eq.vm[0] = &M2_cmt[0] ;  
 
  /* Sac header allocations */
  hdr_alloc(&hdr)  ; /* header for synthetics      */
  hdr_alloc(&hdrg) ; /* header for green functions */

  /* Tapering before the P arrival */
  if(argc == 4)
    {
      if (strcmp(argv[3], "-t") == 0){
	fprintf(stdout, "Tapering between Origin and P_arrival times\n");
	tapering = YES;}
      else
	fprintf(stdout, "WARNING: ignoring unknown option: %s\n", argv[3]);
    }

  /* Reading the CMT file */
  strcpy(eq.cmtfile, argv[1]);
  read_cmt(&eq, &hdr);
   
 
  if(tapering == YES) 
    init_tapering(&nh, &nd, &hdr.evdp, &dv, &tv);

  /* Reading the stations file */
  strncpy(stat_file, argv[2], FSIZE);
  nstat = read_stats(stat_file, &stats, &nets, &stlats, &stlons);
  
  /* Float allocations */
  dists = float_calloc(nstat);
  azs   = float_calloc(nstat);
  bazs  = float_calloc(nstat);
  xdegs = float_calloc(nstat);

  /* Distance, azimuth, back-azimuth, etc */
  distaz((double)hdr.evla, (double)hdr.evlo, stlats, stlons, nstat, 
         dists, azs, bazs, xdegs, &nerr);  
  hdr.o = 0.; /* otime */
    
  /* Treating each station */
  for (jstat=0; jstat < nstat; jstat++)
    {
         printf("%-5s", stats[jstat]);
      /* Rotating the moment tensor */
         rotate_cmt(M1_cmt, M2_cmt, (double)azs[jstat]) ;
         get_prefix(hdr.evdp, xdegs[jstat], segm1, segm2, &best_depth, &best_dist) ;
      /* Computing seismograms */
         hdr.npts  = get_length(segm1, segm2) ;
         if (hdr.npts > __LEN_SIG__) hdr.npts = __LEN_SIG__ ;
         Z  = double_calloc(hdr.npts) ; /* Memory allocation */
	 TH = double_calloc(hdr.npts) ; /* Memory allocation */
	 PH = double_calloc(hdr.npts) ; /* Memory allocation */
         summ_up(M2_cmt, segm1, segm2, GFs, Z, TH, PH, &hdr, &hdrg) ;
      /* true distance correction (group velocity ~ 9km/s)  */
      /*      hdr.b += (xdegs[jstat] - best_dist)*111.1/9.0;*/
         if (tapering == YES)
	   taper(Z, TH, PH, &hdr.npts, &hdr.delta, &(xdegs[jstat]), tv, dv, &nd) ;
      /* Writing the output sac files */
	 hdr.gcarc = xdegs[jstat];
	 hdr.az = azs[jstat];
	 hdr.baz = bazs[jstat];
         save_sac(stats[jstat], nets[jstat], "LHZ", &stlats[jstat] , &stlons[jstat], &hdr, Z) ;
         save_sac(stats[jstat], nets[jstat], "LHL", &stlats[jstat] , &stlons[jstat], &hdr,TH) ;
         save_sac(stats[jstat], nets[jstat], "LHT", &stlats[jstat] , &stlons[jstat], &hdr,PH) ; 
      /* Memory Freeing */
         free((void *)  Z) ; 
	 free((void *) PH) ; 
	 free((void *) TH) ; 
    }

  /* Freeing memory */
  for(j=0; j<10; j++)
    free((void *)GFs[j])  ;
  free((void**)GFs)       ;  
  free((void *)M1_cmt)    ;
  free((void *)M2_cmt)    ;
  free((void**)eq.vm)     ;
  free((void *)dists)     ;
  free((void *)azs)       ;
  free((void *)bazs)      ;
  free((void *)xdegs)     ;
  for (j=0;j< nstat;j++)
    {
      free((void *) stats[j]) ;
      free((void *) nets[j])  ;
    }
  free((void **) stats) ;
  free((void **) nets)  ;
  free((void *) stlats) ;
  free((void *) stlons) ;
  if(tapering == YES) 
    {
      free((void *) tv);
      free((void *) dv);
    }
  printf("\n");
  exit(0);
}


void 
init_tapering(int *nh, int *nd, float *h, double **dv, double **tv)
{
  int    ierror = 1;
  double dep;
  dep = (double)(*h) ;
  
  /* Allocate memory */
  *tv = double_alloc(*nd); /* travel times */
  *dv = double_alloc(*nd); /* distances    */

  /* Set travel time table for depth = dep             */
  /*   ierror  = 1 => error if dep > maxtabledep       */
  /*   ierror != 1 => extrapolate if dep > maxtabledep */
  trav_time_init(nh, nd, &dep, *dv, *tv, &ierror);
}
 
void 
taper(double *R, double *T, double *P, int *npts, float *dt, float *d, double *tv, double *dv, int *nd)
{
  int    n1, n2,j,ierror=1;
  double P_tt, fac, blr, blt, blp, delta, gcarc;
  
  delta = (double)(*dt);
  gcarc = (double)(*d);
  
  /* Set travel time for dist = gcarc                     */
  /*   ierror  = 1 => error if  gcarc> maxtabledist       */
  /*   ierror != 1 => extrapolate if gcarc > maxtabledist */
  trav_time( &gcarc, tv, dv, nd, &P_tt, &ierror);

  n1 = 0;
  n2 = (int)(floor((P_tt-SAFETY_DELAY)/delta));
  if (n2 > *npts) n2 = *npts;
  
  if(n2 > n1+1)
    {
      blr = 0.;
      blt = 0.;
      blp = 0.;
      for (j=n1; j<n2; j++)  {
	blr += R[j] ;
	blt += T[j] ;
	blp += P[j] ;        }
      blr /= (float)(n2-n1) ;
      blt /= (float)(n2-n1) ;
      blp /= (float)(n2-n1) ;
      for (j=0; j<*npts; j++){
	R[j] -= blr ;
	T[j] -= blt ;
	P[j] -= blp ;        }
      
      for (j=n1; j<n2; j++)  {
	fac   = (1.-cos((float)(j-n1)*M_PI/(float)(n2-n1-1)))/2.;
	R[j] *= fac ;
	T[j] *= fac ;
	P[j] *= fac ;        }
    }
  return;
}

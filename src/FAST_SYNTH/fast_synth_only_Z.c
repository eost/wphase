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

void summ_up_only_Z(double* M2_cmt, char* segm1, char* segm2, double** GFs, 
		    double* R, sachdr* hdr, sachdr* hdrg);

void rotate_traces(float* T, float* P, float baz, int  npts, float* N, 
		   float* E);

void save_sac(char *stnm, char *netwk, char *chan, float *lat, float *lon, 
	      sachdr *hdr,  double*   depval);

long int get_length(char *segm1, char *segm2);

/**********************************************************************************************/
/* The variables and functions declared hereafter are used for tapering the non causal signal */
void init_tapering(int *nh, int *nd, float *h, double **dv, double **tv);
void taper(double *Z, int *npts, float *delta, float *gcarc, double *tv, double *dv, int *nd);
/**********************************************************************************************/


int 
main(int argc, char **argv)
{

  /* CMT variables */
  char   *stat_file       ; 
  double *M1_cmt, *M2_cmt ; 
  str_quake_params eq     ;  

  /* Travel times variables */
  int    nh=NDEPTHS, nd=NDISTAS ;
  double *tv, *dv               ;

  /* Stations variables */
  char  **stats, **nets             ;
  float *stlats,*stlons             ;
  float *dists, *azs, *bazs, *xdegs ;
  int   jstat, nstat                ;

  /* Sac header variables */
  long int nerr                  ;
  char     *segm1, *segm2        ;
  float    best_depth, best_dist ;  
  sachdr   hdr, hdrg             ; 

  /* GFs */
  double **GFs ;
  int    j     ;

  /* Seismograms */
  double *Z ;

  /* Flag */
  int tapering = NON ;

  /* Input parameters */
  if(argc != 3 && argc != 4)
    {
      fprintf(stderr,"Error: syntax: fast_synth_only_Z  cmt_file stats_file [-t]\n") ;
      exit(1) ;
    }

  /* String allocations */
  stat_file = char_alloc(FSIZE) ;
  segm1     = char_alloc(FSIZE) ;
  segm2     = char_alloc(FSIZE) ;

  /* Double allocations */
  GFs    = double_alloc2(4,__LEN_SIG__);/* GFs: Rrr, Rtt, Rpp, Rrt */
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
         rotate_cmt(M1_cmt, M2_cmt, (double)azs[jstat]);
         get_prefix(hdr.evdp, xdegs[jstat], segm1, segm2, &best_depth, &best_dist);
      /* Computing the vertical seismograms */
         hdr.npts  = get_length(segm1, segm2);
         if (hdr.npts > __LEN_SIG__) hdr.npts = __LEN_SIG__ ;
         Z = double_calloc(hdr.npts); /* Memory allocation */
         summ_up_only_Z(M2_cmt, segm1, segm2, GFs, Z, &hdr, &hdrg);
      /* true distance correction (group velocity ~ 9km/s)  */
      /*      hdr.b += (xdegs[jstat] - best_dist)*111.1/9.0;*/
         if (tapering == YES)  
	   taper(Z, &hdr.npts, &hdr.delta, &(xdegs[jstat]), tv, dv, &nd);
      /* Writing the output sac files */
         save_sac(stats[jstat], nets[jstat], "LHZ", &stlats[jstat] , &stlons[jstat], &hdr,Z); 
      /* Memory Freeing */
         free((void *)Z  ); 
    }

  /* Freeing memory */
  for(j=0; j<4; j++)
    free((void *)GFs[j])  ;
  free((void**)GFs)       ;  
  free((void *)M1_cmt)    ;
  free((void *)M2_cmt)    ;
  free((void**)eq.vm)     ;
  free((void *)stat_file) ;
  free((void *)segm1)     ;
  free((void *)segm2)     ;
  free((void *)dists)     ;
  free((void *)azs)       ;
  free((void *)bazs)      ;
  free((void *)xdegs)     ;

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
taper(double *Z, int *npts, float *dt, float *d, double *tv, double *dv, int *nd)
{
  int    n1, n2,j,ierror=1;
  double P_tt, fac, bl, delta, gcarc;
  
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
      bl = 0.;
      for (j=n1; j<n2; j++)
	{
	  bl += Z[j];
	}
      bl /= (float)(n2-n1);
      for (j=0; j<*npts; j++)
	{
	  Z[j] -= bl;
	}
      
      for (j=n1; j<n2; j++)
	{
	  fac   = (1.-cos((float)(j-n1)*M_PI/(float)(n2-n1-1)))/2.;
	  Z[j] *= fac;
	}
    }
  return;
}

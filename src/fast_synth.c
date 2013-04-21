/***************************************************************************
*
*	              W phase source inversion package 	            
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Université de Strasbourg / CNRS 
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
#include "rwsacs.h"
#include "proto_alloc.h"
#include "read_i_files.h"
#include "travel_times.h"

#define  YES  1
#define  NON  0


/* External routines */
int  read_stats(char* file, char*** stats, char*** nets, float** stlats, 
		float** stlons);

void init_tapering(int nh, int nd, double h, double **dv, double **tv);

void taper_syn(double *R, double *T, double *P, int *npts, float *dt, float *d, 
	   double *tv, double *dv, int nd);

void distaz(double cmt_lat, double cmt_lon, float* stlats, float* stlons, 
	    int nstat, float* dists, float* azs, float* bazs, float* xdegs,
	    long int* nerr);  

void save_sac(char *stnm, char *netwk, char *chan, float *lat, float *lon, 
	      sachdr *hdr,  double*   depval);

void fast_synth_sub(double az, double baz, double xdeg, double *tv, double *dv, int nd, 
		    str_quake_params *eq, sachdr *hdr, double **GFs, double *Z, double *TH, double *PH);


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
  sachdr   hdr                        ; 
  /* Seismograms */
  double **GFs ;
  double *Z, *TH, *PH ;
  /* Tapering */
  int    tapering = NON, ierror = 1 ;
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
  Z   = double_alloc(__LEN_SIG__) ;   /* Vertical components     */
  TH  = double_alloc(__LEN_SIG__) ;   /* Radial components       */
  PH  = double_alloc(__LEN_SIG__) ;   /* Transverse components   */
  M1_cmt = double_alloc(6)   ; /* Moment tensor           */  
  M2_cmt = double_alloc(6)   ; /* Rotated moment tensor   */
  eq.vm  = double_alloc2p(2) ;
  eq.vm[1] = &M1_cmt[0] ;  /* eq.vm[1] is extracted from cmtfile */
  eq.vm[0] = &M2_cmt[0] ;  
  tv = double_alloc(nd); /* travel times */
  dv = double_alloc(nd); /* distances    */
 
  /* Sac header allocation */
  hdr_alloc(&hdr)  ; /* header for synthetics      */

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
  get_cmtf(&eq, 2) ;

  /* Set travel time table for depth = dep */
  trav_time_init(nh,nd,eq.evdp,dv,tv,&ierror);

  /* Reading the stations file */
  strncpy(stat_file, argv[2], FSIZE);
  nstat = read_stats(stat_file, &stats, &nets, &stlats, &stlons);
  
  /* Float allocations */
  dists = float_calloc(nstat);
  azs   = float_calloc(nstat);
  bazs  = float_calloc(nstat);
  xdegs = float_calloc(nstat);

  /* Distance, azimuth, back-azimuth, etc */
  distaz(eq.evla, eq.evlo, stlats, stlons, nstat, dists, azs, bazs, xdegs, &nerr);  
  hdr.o = 0.; /* otime */
    
  /* Treating each station */
  for(jstat=0;jstat< nstat;jstat++)
    {
      printf("%-5s", stats[jstat]) ;
      fast_synth_sub(azs[jstat], bazs[jstat], xdegs[jstat], tv, dv, nd, &eq, &hdr, GFs, Z, TH, PH);
      if (tapering == YES)
	taper_syn(Z,TH,PH,&hdr.npts,&hdr.delta,&(xdegs[jstat]),tv,dv,nd) ;
      save_sac(stats[jstat], nets[jstat], "LHZ", &stlats[jstat] , &stlons[jstat], &hdr, Z) ;
      save_sac(stats[jstat], nets[jstat], "LHL", &stlats[jstat] , &stlons[jstat], &hdr,TH) ;
      save_sac(stats[jstat], nets[jstat], "LHT", &stlats[jstat] , &stlons[jstat], &hdr,PH) ; 
    }

  /* Freeing memory */
  free((void *)  Z) ; 
  free((void *) PH) ; 
  free((void *) TH) ; 
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

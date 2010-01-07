#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <locale.h>

/* Subroutines headers */
#include "proto_alloc.h"    /* proto_alloc.c    */
#include "rwsacs.h"      /* rwsacs.c      */
#include "rwtextfiles.h" /* rwtextfiles.c */
#include "sort_tree.h"   /* sort_tree.h   */
#include "travel_times.h"
#include "read_i_files.h"

#ifndef SAMPLEPERIOD
#define	SAMPLEPERIOD 1.
#endif /* not SAMPLEPERIOD */


/* External routines */
void distaz(double    cmt_lat,  double    cmt_lon,  float*    stlats,  float*    stlons, 
	    int       nstat,    float*    dists,    float*    azs,     float*    bazs,   
	    float*    xdegs,    long int* nerr);  


/* Internal routines */
void get_params(int argc, char **argv, int *un, char **i_sacs, char **o_sacs, 
		str_quake_params *eq) ;
void delta_t2(int *y1, int *mo1, int *j1, int *h1, int *mi1, int *s1, int *ms1,
	      int *y0, int *j0, int *h0, int *m0, int *s0, int *ms0, double *tdiff) ;


	  
int 
main(int argc, char *argv[])
{
  int    nc, nh=NDEPTHS, nd=NDISTAS ; 
  int    tmp, ierror = 1, tterr = 0 ; 
  int    sampstart, flag            ;
  long   int nerr                   ;
  float  dist, az, baz, xdeg        ;
  double xdegd, P_tt, tdiff, otime  ;
  double *x_in, *x_out, *dv, *tv    ;
  char   *i_sacs, *fil,*o_sacs      ;
  FILE   *i_sacf, *o_sacf   ;
  str_quake_params   eq     ;
  sachdr hdr                ;
  struct tree   *root, *mod ;
  
  /* OPTIONS */
  int un ;

  /* Input params */
  get_params(argc, argv, &un, &i_sacs, &o_sacs, &eq)     ;

  /* Allocate memory */
  fil    = char_alloc(FSIZE) ;
  x_in   = double_alloc((int)__LEN_SIG__);
  dv     = double_alloc(nd);
  tv     = double_alloc(nd);
  hdr_alloc(&hdr)    ;  
  root = alloctree();
  mod  = alloctree();

  /* Set travel time table */
  trav_time_init(&nh, &nd, &eq.pde_evdp, dv, tv, &tterr) ;

  /* Read sac file list */
  i_sacf = openfile_rt(i_sacs,&nc) ;
  /* ofil   = char_alloc2(nc, FSIZE)  ; */
  nc = 0   ;
  flag = 0 ;
  while( (tmp=fscanf (i_sacf, "%s", fil)) != EOF )
    {
      check_scan(1, tmp, i_sacs, i_sacf);

      /* Read sac file */
      rhdrsac(fil, &hdr, &ierror)       ;
      if (hdr.npts > (int)__LEN_SIG__)
	hdr.npts = (int)__LEN_SIG__ ;
      rdatsac(fil, &hdr, x_in, &ierror) ;

      /* Set epicentral distances, azimuth, backazimuth */
      dist = 0. ;
      az   = 0. ;
      baz  = 0. ;
      xdeg = 0. ;
      distaz(eq.pde_evla, eq.pde_evlo, &hdr.stla, &hdr.stlo, tmp, &dist, &az, &baz, &xdeg, &nerr) ;

      /* Set travel time */
      xdegd = (double) xdeg                         ;
      trav_time(&xdegd, tv, dv, &nd, &P_tt, &tterr) ;

      /* Set event origin time in sac header variable 'o' (relative to the reference time) */
      delta_t2(&eq.ot_ye, &eq.ot_mo, &eq.ot_dm, &eq.ot_ho, &eq.ot_mi, &eq.ot_se, &eq.ot_ms,
	       &hdr.nzyear, &hdr.nzjday, &hdr.nzhour, &hdr.nzmin, &hdr.nzsec, &hdr.nzmsec, &tdiff) ;
      otime = tdiff ;
      hdr.t[0]  = (float)(P_tt + tdiff) ; /* P arrival         */

      /* Windowing -- Screening by distance */
      tdiff    += P_tt - (double)hdr.b - eq.preevent - (double)SAFETY_DELAY ;      /* time for the 1st sample */
      sampstart = (int) (tdiff/((double)hdr.delta) + 0.5)                   ;
      if ((int)(hdr.delta * 1000.+0.5) != (int)((float)SAMPLEPERIOD * 1000.+0.5)) /* *************** MAY BE MODIFIED ***************** */
	{
	  fprintf(stderr, "WARNING: non uniform samp. period between sac files\n") ;
	  fprintf(stderr, "     ...file : %s with dt = %e is rejected\n", fil, hdr.delta)  ;
	  flag++   ;
	  continue ;
	}
      else if (sampstart >= 0 && sampstart < hdr.npts)
	{
	  /* Set new sac header variables */
	  hdr.delta = (float) SAMPLEPERIOD  ;
	  hdr.o     = (float) otime         ;   /* event origin time */
	  //hdr.npts  = hdr.npts + 1 - sampstart;              /* Error **** nb of samples             */
	  hdr.npts  = hdr.npts - sampstart ;                 /* nb of samples (corrected) */
	  hdr.b     = hdr.b + ((float)sampstart)*hdr.delta ; /* shift of the first sample (corrected) */
	  hdr.e     = hdr.b + (hdr.npts-1) * hdr.delta    ; 
	  hdr.dist  = dist ; /* epicentral distance (km)                      */
	  hdr.gcarc = xdeg ; /* station to event great circle arc length(deg) */
	  hdr.az    = az   ; /* event to station azimuth (deg)                */
	  hdr.baz   = baz  ; /* station to event azimuth (deg)                */
	  hdr.evla  = (float) eq.pde_evla ;
	  hdr.evlo  = (float) eq.pde_evlo ;
	  hdr.evdp  = (float) eq.pde_evdp ;
	  
	  /* Write trimmed sac file */
	  if (hdr.npts > (int)__LEN_SIG__)
	    fprintf(stderr,"Warning : traces cut to %d samples.\n", (int)__LEN_SIG__);
	  whdrsac(fil, &hdr);
	  //x_out = &x_in[sampstart-1]; /* Error **** shift of the first sample */
	  x_out = &x_in[sampstart];   /* shift of the first sample (corrected) */	  
	  wdatsac(fil, &hdr, x_out);

	  /* Sort sac files in a tree */
	  if (nc==0)
	    splithdr(fil, &hdr, root);
	  else
	    {
	      splithdr(fil, &hdr, mod) ;
	      build (root, mod, un)     ;
	    }
	  nc++ ;
	  /* strcpy(ofil[nc++],fil) ; */
	}
      else 
	{
	  fprintf(stderr,"Warning: (trim_sac_files) rejected incomplete file : %s (%d -- %d)\n", 
		  fil, sampstart, hdr.npts);
	}
    }
  fclose(i_sacf);
  if (flag > 1)
    {
      fprintf(stderr, "WARNING: %d files have been rejected because of non uniform sampling period\n",
	      flag);
      fprintf(stderr, "    ...: if to much files are rejected, please screen or decimate data files manually\n");
    }
  /* i_sacf = openfile_wt(i_sacs)     ; */
  /* for (i=0; i<nc; i++)               */
  /*   fprintf(i_sacf,"%s\n",ofil[i]) ; */
  /* fclose(i_sacf)                   ; */

  o_sacf = openfile_wt(o_sacs);
  savetree(root, o_sacf, &eq.dmin, &eq.dmax);
  fclose(o_sacf);

  /* Memory Freeing */
  free((void *)i_sacs) ;
  free((void *)o_sacs) ;
  free((void *)x_in)   ;
  free((void *) dv)    ;
  free((void *) tv)    ;
  freetree(root);
  freetree(mod);
  return 0;
}


void 
get_params(int argc, char **argv, int *un, char **i_sacs, char **o_sacs, 
	   str_quake_params *eq)
{
  int  i              ;
  char **keys, *i_tmp ;
  
  if ( argc < 4 ) {
    fprintf (stderr, "*** ERROR (minimum of 3 params needed (%d given)) ***",argc-1)        ;
    fprintf (stderr, "Syntax : %s i_master(in) i_sac_list(in) o_sac_list(out) [-u(unique network)  -a(all channels)]\n", argv[0]) ;
    exit(1) ; }

  *un = 0 ; 
  if (!strncmp(argv[argc-1],"-u",2))
    *un = 1  ;
  else if (!strncmp(argv[argc-1],"-a",2))
    *un = -1 ;

  i_tmp     = char_alloc(FSIZE) ;  
  (*i_sacs) = char_alloc(FSIZE) ;
  (*o_sacs) = char_alloc(FSIZE) ;
  strcpy(   i_tmp , argv[1]) ;
  strcpy((*i_sacs), argv[2]) ;
  strcpy((*o_sacs), argv[3]) ;
  eq->cmtfile[0] = '\0';

  i = 0;
  keys = char_alloc2(4, 16)   ;
  strcpy(keys[i++],"CMTFILE") ;
  strcpy(keys[i++],"IDEC_2")  ;
  strcpy(keys[i++],"DMIN")    ;
  strcpy(keys[i++],"DMAX")    ;

  get_i_master(i_tmp, keys, 4, eq) ;  
  get_cmtf( eq, 0) ;
  
  /* Memory Freeing */
  for(i=0 ; i<4 ; i++)
    free((void*)keys[i]) ;
  free((void**) keys )   ;
  free((void*) i_tmp )   ;
}


void 
date2epoch1(int *year, int *mont, int *jday, int *hour, int *min, int *sec, int *msec, double *epoch)
{
  struct tm date ;
  time_t    tmp  ;
  extern long timezone   ;
  date.tm_sec    = *sec  ;  
  date.tm_min    = *min  ;
  date.tm_hour   = *hour ;
  date.tm_mday   = *jday ;
  date.tm_mon    = *mont-1      ;
  date.tm_year   = *year - 1900 ;
  date.tm_isdst  =  0 ;
  tzset() ;
  tmp            = mktime(&date) ;
  *epoch         = tmp + (double)(*msec)/1000. - timezone ;
/*   
	  printf("in date2epoch: %4d %1d %3d %2d %2d %2d %2d %2d %1d %1ld %s = %lf\n",
	  date.tm_year, date.tm_wday, date.tm_yday, date.tm_hour, date.tm_min, 
	  date.tm_sec,  date.tm_mday, date.tm_mon,  date.tm_isdst, date.tm_gmtoff, 
	  date.tm_zone, *epoch); 
*/ 
}



void 
date2epoch2(int *year, int *jday, int *hour, int *min, int *sec, int *msec, double *epoch)
{
  struct tm date;
  time_t    tmp;
  extern long timezone;
  date.tm_sec    = *sec  ;  
  date.tm_min    = *min  ;
  date.tm_hour   = *hour ;
  date.tm_mday   = *jday ;
  date.tm_mon    =  0    ;
  date.tm_year   = *year - 1900 ;
  date.tm_wday   =  0    ;
  date.tm_yday   = *jday - 1 ;
  date.tm_isdst  =  0    ;
  //date.tm_gmtoff =  0;
  tzset();
  tmp            = mktime(&date) ;
  *epoch         = tmp + (double)(*msec)/1000. - timezone ;
/*   
	  printf("in date2epoch: %4d %1d %3d %2d %2d %2d %2d %2d %1d %1ld %s = %lf\n",
	  date.tm_year, date.tm_wday, date.tm_yday, date.tm_hour, date.tm_min, 
	  date.tm_sec,  date.tm_mday, date.tm_mon,  date.tm_isdst, date.tm_gmtoff, 
	  date.tm_zone, *epoch); 
*/ 
}



void 
delta_t2(int *y1, int *mo1, int *j1, int *h1, int *mi1, int *s1, int *ms1,
	 int *y0, int *j0, int *h0, int *m0, int *s0, int *ms0, double *tdiff)
{
  double t1, t0 ;
  int    one    ;
  one = 1       ;
  date2epoch1(y1,mo1,j1,h1,mi1,s1,ms1,&t1) ;
  date2epoch2(y0,&one,h0,m0,s0,ms0,&t0)    ;
  *tdiff = t1 - (t0 + (*j0-1)*86400)       ;
}

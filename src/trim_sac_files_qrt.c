#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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

#ifndef SAFETY_DIST
#define	SAFETY_DIST __SAFETY_DIST__
#endif /* not SAFETY_DIST */


/* External routines */
void distaz(double    cmt_lat,  double    cmt_lon,  float*    stlats,  float*    stlons, 
	    int       nstat,    float*    dists,    float*    azs,     float*    bazs,   
	    float*    xdegs,    long int* nerr);  
int  rmseed(str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1, double *o_x, sachdr *hdr);


/* Internal routines */
void get_params(int argc, char **argv, int *un, char **i_locs, char **o_sacdir, 
		char **o_sacs, str_quake_params *eq) ;
void delta_t2(int *y1, int *mo1, int *j1, int *h1, int *mi1, int *s1, int *ms1,
	      int *y0, int *j0, int *h0, int *m0, int *s0, int *ms0, double *tdiff) ;
void wp_end(double gcarc, double *wp_win4, double *twp_end);
void complete_with_blank(char *field, int N);
void date2epoch1(int *year, int *mont, int *jday, int *hour, int *min, int *sec, int *msec, double *epoch);
	  
int 
main(int argc, char *argv[])
{
  int    nc, nh=NDEPTHS, nd=NDISTAS ; 
  int    tmp, tterr = 0 ; 
  int    sampstart, flag        ;
  long   int nerr                   ;
  float  dist, az, baz, xdeg        ;
  double xdegd, P_tt, tdiff, otime  ;
  double *x_in, *x_out, *dv, *tv    ;
  char   *i_locs, *msfil, *scfil, *o_sacdir, *o_sacs ;
  FILE   *i_locf, *o_sacf ;
  str_quake_params   eq     ;
  sachdr hdr                ;
  time_t      t0, t1 ;
  struct tm   *tm0   ;
  struct tree *root, *mod ;
  
  /* OPTIONS */
  int un ;
  
  /* Input params */
  get_params(argc, argv, &un, &i_locs, &o_sacdir, &o_sacs, &eq) ;

  /* Allocate memory */
  msfil  = char_alloc(FSIZE) ;
  scfil  = char_alloc(FSIZE) ;
  x_in   = double_alloc((int)__LEN_SIG__);
  dv     = double_alloc(nd);
  tv     = double_alloc(nd);
  hdr_alloc(&hdr)    ;  
  root = alloctree();
  mod  = alloctree();

  /* Set travel time table */
  trav_time_init(&nh, &nd, &eq.pde_evdp, dv, tv, &tterr) ;
  /* Open loc list */
  i_locf = openfile_rt(i_locs,&nc) ;
  /* ofil   = char_alloc2(nc, FSIZE)  ; */
  nc = 0   ;
  flag = 0 ;
  while( (tmp=fscanf (i_locf, "%s %s %s %s %f %f %f", hdr.knetwk,hdr.kstnm,hdr.kcmpnm,hdr.khole,\
		      &hdr.stla,&hdr.stlo,&hdr.stel)) != EOF )
    {
      check_scan(7, tmp, i_locs, i_locf) ; 

      /* Set epicentral distances, azimuth, backazimuth */
      dist = 0. ;
      az   = 0. ;
      baz  = 0. ;
      xdeg = 0. ;
      distaz(eq.pde_evla, eq.pde_evlo, &hdr.stla, &hdr.stlo, 1, &dist, &az, &baz, &xdeg, &nerr) ;
      xdegd = (double) xdeg ;
      
      /* Rough distance screening */
      if (xdegd < eq.dmin - SAFETY_DIST || xdegd >eq.dmax + SAFETY_DIST)
	{
	  fprintf(stderr,"Warning: (trim_sac_files_qrt) sta %s rejected (gcarc= %f deg)\n",  
		  hdr.kstnm, xdegd); 
	  continue;
	}
      
      /* Set wp time window */
      trav_time(&xdegd, tv, dv, &nd, &P_tt, &tterr) ; /* P travel time */
      date2epoch1(&eq.ot_ye, &eq.ot_mo, &eq.ot_dm, &eq.ot_ho, &eq.ot_mi, &eq.ot_se, &eq.ot_ms, &otime) ;
      t0  = (time_t) (otime + P_tt - eq.preevent - (double)SAFETY_DELAY) ;
      tm0 = gmtime(&t0) ;
      wp_end(xdegd, eq.wp_win4, &tdiff);
      t1  = (time_t) (otime + P_tt + tdiff + eq.preevent + (double)SAFETY_DELAY + 1.) ;

      /* Read mseed file */
      add_slash(eq.seed);
      sprintf(scfil,"%s%4d.%03d.%02d.%02d.%02d.%04d.%s.%s.%s.%s.SAC",o_sacdir,tm0->tm_year+1900,
	      tm0->tm_yday+1,tm0->tm_hour,tm0->tm_min,tm0->tm_sec,0,hdr.knetwk,hdr.kstnm,hdr.khole,hdr.kcmpnm);
      if (rmseed(&eq, tm0, (int)t0, (int)t1, x_in, &hdr))
	  continue;

      /* Write rough sac file */
      whdrsac(scfil, &hdr);
      wdatsac(scfil, &hdr, x_in);

      /* Set event origin time in sac header variable 'o' (relative to the reference time)  */
      delta_t2(&eq.ot_ye, &eq.ot_mo, &eq.ot_dm, &eq.ot_ho, &eq.ot_mi, &eq.ot_se, &eq.ot_ms, 
	       &hdr.nzyear, &hdr.nzjday, &hdr.nzhour, &hdr.nzmin, &hdr.nzsec, &hdr.nzmsec, &tdiff) ; 
      otime = tdiff ; 
      hdr.t[0]  = (float)(P_tt + tdiff) ; /* P arrival */

      /* Windowing -- Screening by distance */
      tdiff   += P_tt - (double)hdr.b - eq.preevent - (double)SAFETY_DELAY ;  /* time for the 1st sample */
      sampstart = (int) (tdiff/((double)hdr.delta) + 0.5) ;
      if ((int)(hdr.delta * 1000.+0.5) != (int)((float)SAMPLEPERIOD * 1000.+0.5)) /* *************** MAY BE MODIFIED ***************** */
	{ 
	  fprintf(stderr, "WARNING: non uniform samp. period between sac files\n") ; 
	  fprintf(stderr, "     ...file : %s with dt = %e is rejected\n", scfil, hdr.delta)  ; 
	  flag++   ; 
	  continue ; 
	} 
      else if (sampstart >= 0 && sampstart < hdr.npts) 
	{ 
	  /* Set new sac header variables */ 
	  hdr.delta = (float) SAMPLEPERIOD  ; 
	  hdr.o     = (float) otime         ; /* event origin time */
	  //hdr.npts  = hdr.npts + 1 - sampstart;              /* Error **** nb of samples     */
	  //hdr.b     = hdr.b + ((float)sampstart)*hdr.delta ; /* Error **** shift of the first sample */
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
	  strcat(scfil,".scr.sac");
	  whdrsac(scfil, &hdr);
	  //x_out = &x_in[sampstart-1]; /* Error **** shift of the first sample */
	  x_out = &x_in[sampstart];   /* shift of the first sample (corrected) */	  
	  wdatsac(scfil, &hdr, x_out);


	  /* Sort sac files in a tree */
	  if (nc==0)
	    splithdr(scfil, &hdr, root);
	  else
	    {
	      splithdr(scfil, &hdr, mod) ;
	      build (root, mod, un)     ;
	    }

	  nc++ ;
	  /* strcpy(ofil[nc++],fil) ; */
	}
      else
	{
	  fprintf(stderr,"Warning: (trim_sac_files_qrt) rejected incomplete file : %s (%d -- %d)\n",
		  scfil, sampstart, hdr.npts);
	}
    }
  fclose(i_locf);
  if (flag > 1)
    {
      fprintf(stderr, "Warning: %d files have been rejected because of non uniform sampling period\n",
	      flag);
      fprintf(stderr, "    ...: if to much files are rejected, please screen or decimate data files manually\n");
    }

  o_sacf = openfile_wt(o_sacs);
  savetree(root, o_sacf, &eq.dmin, &eq.dmax);
  fclose(o_sacf);

  /* Memory Freeing */
  free((void *)i_locs) ;
  free((void *)msfil)  ;
  free((void *)scfil)  ;
  free((void *)o_sacs) ;
  free((void *)o_sacdir) ;
  free((void *)x_in)   ;
  free((void *) dv)    ;
  free((void *) tv)    ;  
  free((void *)eq.wp_win4);
  freetree(root);
  freetree(mod);
  return 0;
}


/********************************************************/
/*  W-phase time window used in the data fit.           */ 
/*  The epicentral distance is given in degrees and     */
/*  twp_end is given in seconds with respect to the     */
/*  P arrival time.                                     */
/*  The window is defined with 1, 2, 3 or 4 parameters  */
void 
wp_end(double gcarc, double *wp_win4, double *twp_end)
{
  double feakdist;

  if (wp_win4[2] > gcarc) 
    feakdist = wp_win4[2] ;
  else
    feakdist = gcarc ;
  if (wp_win4[3] < feakdist) 
    feakdist = wp_win4[3] ;

  *twp_end = wp_win4[1] * feakdist ;
}


void 
complete_with_blank(char *c, int N)
{
  int i,n ;
  n = strlen(c);
  for(i=n ; i<(N-1) ; i++)
    c[i]=' ' ;
  c[N-1]='\0';
}


void 
get_params(int argc, char **argv, int *un, char **i_locs, char **o_sacdir, char **o_sacs, str_quake_params *eq)
{
  int  i              ;
  char **keys, *i_tmp ;  

  /* Check syntax    */
  if ( argc < 5 ) {
    fprintf (stderr, "*** ERROR (minimum of 4 params needed (%d given)) ***\n",argc-1)        ;
    fprintf (stderr, "Syntax : %s i_master(in) i_loc_lst(in) o_sac_dir o_sac_lst(out) [-u]\n", argv[0]) ;
    exit(1) ;     }

  *un = 0 ; 
  if (!strncmp(argv[argc-1],"-u",2))
    *un = 1  ;
  else if (!strncmp(argv[argc-1],"-a",2))
    *un = -1 ;
  
  /* Allocate memory */
  i_tmp       = char_alloc(FSIZE) ;
  (*i_locs)   = char_alloc(FSIZE) ; 
  (*o_sacdir) = char_alloc(FSIZE) ; 
  (*o_sacs)   = char_alloc(FSIZE) ;
  eq->wp_win4 = double_alloc(4)   ;
  eq->cmtfile[0] = '\0';

  /* Get args        */
  strcpy(     i_tmp , argv[1])    ;
  strcpy((  *i_locs), argv[2])    ;
  strcpy((*o_sacdir), argv[3])    ;
  strcpy((  *o_sacs), argv[4])    ;

  add_slash(*o_sacdir);

  /* Read i_masterfile and cmtfile */
  i = 0;
  keys = char_alloc2(6, 16)   ;
  strcpy(keys[i++],"CMTFILE") ;
  strcpy(keys[i++],"SEED") ;
  strcpy(keys[i++],"IDEC_2")  ;
  strcpy(keys[i++],"WP_WIN")  ;
  strcpy(keys[i++],"DMIN")    ;
  strcpy(keys[i++],"DMAX")    ;
  get_i_master(i_tmp, keys, 6, eq) ;  
  get_cmtf( eq, 0) ;

  /* Memory Freeing */
  for(i=0 ; i<6 ; i++)
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

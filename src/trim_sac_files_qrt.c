/****************************************************************
*	W phase package - Trim from miniseed records
*                                           
*       History
*             2010  Original Coding          Zacharie Duputel
*       License
*             Distributed under the same License as the libmseed 
*             source and binaries, used here under permission by 
*             Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <time.h>
#include <locale.h>

/* Subroutines headers */
#include "proto_alloc.h" /* proto_alloc.c */
#include "rwsacs.h"      /* rwsacs.c      */
#include "rwtextfiles.h" /* rwtextfiles.c */
#include "sort_tree.h"   /* sort_tree.h   */
#include "travel_times.h"
#include "read_i_files.h"
#include "decimate.h"

#ifndef SAFETY_DIST
#define	SAFETY_DIST __SAFETY_DIST__
#endif /* not SAFETY_DIST */


/* External routines */
void distaz(double    cmt_lat,  double    cmt_lon,  float*    stlats,  float*    stlons, 
	    int       nstat,    float*    dists,    float*    azs,     float*    bazs,   
	    float*    xdegs,    long int* nerr);  
int rmseed(str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1, double *o_x, sachdr *hdr);
int rmseed_file(char *mseedfil, str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1,
	      double *o_x, sachdr *o_hdr);
int get_channel_name(char *msfil,char *netwk,char *stanm,char *loc, char *cmpnm);
int init_FIR(double *coeffs, int Ncoeffs, FIR_filter *FIR);
int decimate(FIR_filter *FIR, int dec_fac, double *yi, int ni, double *yo, int *no);

/* Internal routines */
void makeid(sachdr *hdr,char *id) ;
int  findid(char *id,char **ids,int n) ;
void get_params(int argc, char **argv, int *un, char **i_locs, char **i_mseeds, 
		char **o_sacdir, char **o_sacs, str_quake_params *eq);
void date2epoch(int year, int mm, int dd, int hour, int min, int sec, int msec, 
		double *epoch) ;
void delta_t(int y1, int j1, int h1, int m1, int s1, int ms1,int y0, int j0, 
	     int h0, int m0, int s0, int ms0, double *tdiff) ;
void ymd2jul(int yyyy,int mm,int dd, int *jul) ;
void wp_end(double gcarc, double *wp_win4, double *twp_end);
void complete_with_blank(char *field, int N);
	  
int 
main(int argc, char *argv[])
{
  int    i,nids,nc,nh=NDEPTHS,nd=NDISTAS ; 
  int    tmp, tterr = 0 ; 
  int    ot_jday, sampstart, flag   ;
  long   int nerr                   ;
  float  dist, az, baz, xdeg        ;
  float  *stla, *stlo, *stel        ;
  double xdegd, P_tt, tdiff, otime  ;
  double *x_in, *x_out, *dv, *tv    ;
  char *i_locs,*i_mseeds,msfil[FSIZE],scfil[FSIZE] ;
  char *o_sacdir,*o_sacs,**ids,id[IDSIZE] ;
  FILE   *fd, *o_sacf ;
  str_quake_params   eq ;
  sachdr hdr ;
  FIR_filter FIR2,FIR4,FIR5;
  time_t      t0, t1 ;
  struct tm   *tm0   ;
  struct tree *root, *mod ;
  
  /* OPTIONS */
  int un ;
  
  /* Input params */
  get_params(argc, argv, &un, &i_locs, &i_mseeds, &o_sacdir, &o_sacs, &eq) ;

  /* Allocate memory */
  x_in   = double_alloc(__LEN_SIG__*100); 
  dv     = double_alloc(nd);
  tv     = double_alloc(nd);
  hdr_alloc(&hdr)    ;  
  root = alloctree();
  mod  = alloctree();
  init_FIR(dec2,Ndec2,&FIR2) ;
  init_FIR(dec4,Ndec4,&FIR4) ;
  init_FIR(dec5,Ndec5,&FIR5) ;
  /* Set travel time table */
  trav_time_init(&nh, &nd, &eq.pde_evdp, dv, tv, &tterr) ;
  /* Open loc list */
  fd = openfile_rt(i_locs,&nids) ;
  ids  = char_alloc2(nids,IDSIZE)  ;
  stla = float_alloc(nc) ;
  stlo = float_alloc(nc) ;
  stel = float_alloc(nc) ;
  for(i=0;i<nids;i++)
    {
      tmp=fscanf (fd,"%s %s %s %s %f %f %f",hdr.knetwk,hdr.kstnm,hdr.kcmpnm,
		                                hdr.khole,stla+i,stlo+i,stel+i);
      check_scan(7,tmp,i_locs,fd) ; 
      makeid(&hdr,ids[i]);
    }
  fclose(fd);
  
  /* ofil   = char_alloc2(nc, FSIZE) ; */
  fd = openfile_rt(i_mseeds,&nc) ;
  nc = 0   ;
  flag = 0 ;
  while( (tmp=fscanf (fd,"%s",msfil)) != EOF )
    { 
      check_scan(1,tmp,i_mseeds,fd) ;       
      /* Channel location */
      if(get_channel_name(msfil,hdr.knetwk,hdr.kstnm,hdr.khole,hdr.kcmpnm)) /* get mseed channel name */
	continue; 
      makeid(&hdr,id);
      i = findid(id,ids,nids);  
      if (i<0)
	{
	  fprintf(stderr,"Warning: No location for %s\n",id);
	  continue;
	}
      hdr.stla = stla[i] ;
      hdr.stlo = stlo[i] ;
      hdr.stel = stel[i] ;
      dist = 0. ;
      az   = 0. ;
      baz  = 0. ;
      xdeg = 0. ;
      distaz(eq.pde_evla,eq.pde_evlo,&hdr.stla,&hdr.stlo,1,&dist,&az,&baz,&xdeg,&nerr) ;
      xdegd = (double) xdeg ;
      if (xdegd < eq.dmin - SAFETY_DIST || xdegd >eq.dmax + SAFETY_DIST) /* Rough distance screening */
	{
	  fprintf(stderr,"Warning: (trim) channel %s rejected (gcarc= %f deg)\n", id, xdegd); 
	  continue;
	}
      /* Set wp time window */
      trav_time(&xdegd, tv, dv, &nd, &P_tt, &tterr)   ; /* P travel time */
      ymd2jul(eq.ot_ye, eq.ot_mo, eq.ot_dm, &ot_jday) ;
      date2epoch(eq.ot_ye, eq.ot_mo, eq.ot_dm, eq.ot_ho, eq.ot_mi, eq.ot_se, eq.ot_ms, &otime) ;
      t0  = (time_t) (otime + P_tt - eq.preevent - (double)SAFETY_DELAY) ; 
      tm0 = gmtime(&t0) ;                                                  
      wp_end(xdegd, eq.wp_win4, &tdiff);
      t1  = (time_t) (otime + P_tt + tdiff + eq.preevent + (double)SAFETY_DELAY + 1.) ;
      sprintf(scfil,"%s%4d.%03d.%02d.%02d.%02d.%04d.%s.%s.%s.%s.SAC",o_sacdir,tm0->tm_year+1900,
	      tm0->tm_yday+1,tm0->tm_hour,tm0->tm_min,tm0->tm_sec,0,hdr.knetwk,hdr.kstnm,hdr.khole,
	      hdr.kcmpnm); 
      /* READ MINISEED FILE */
      if (rmseed_file(msfil,&eq, tm0, (int)t0, (int)t1, x_in, &hdr)) 
	  continue;      
      delta_t(eq.ot_ye, ot_jday, eq.ot_ho, eq.ot_mi, eq.ot_se, eq.ot_ms,
	      hdr.nzyear, hdr.nzjday, hdr.nzhour, hdr.nzmin, hdr.nzsec, hdr.nzmsec, &tdiff) ;
      otime = tdiff ;       /* Set the sac header variable 'o' */
      hdr.o     = (float) otime ;
      hdr.t[0]  = (float)(P_tt + tdiff) ; /* P arrival time    */
      whdrsac(scfil, &hdr); /* Write rough sac file            */
      wdatsac(scfil, &hdr, x_in);
      /* Decimation */
      if ((int)(hdr.delta * 10000.+0.5) == 100)           /* 100 sps */
	{
	  decimate(&FIR5,5,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR5,5,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR4,4,x_in,hdr.npts,x_in,&hdr.npts);
	  hdr.delta = 1.0;
	}
      else if ((int)(hdr.delta * 10000.+0.5) == 125)      /*  80 sps */
	{
	  decimate(&FIR5,5,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR4,4,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR4,4,x_in,hdr.npts,x_in,&hdr.npts);
	  hdr.delta = 1.0;
	}
      else if ((int)(hdr.delta * 10000.+0.5) == 250)      /*  40 sps */
	{
	  decimate(&FIR5,5,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR4,4,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR2,2,x_in,hdr.npts,x_in,&hdr.npts);
	  hdr.delta = 1.0;
	}
      else if ((int)(hdr.delta * 10000.+0.5) == 400)      /*  25 sps */
	{
	  decimate(&FIR5,5,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR5,5,x_in,hdr.npts,x_in ,&hdr.npts);
	  hdr.delta = 1.0;
	}
      else if ((int)(hdr.delta * 10000.+0.5) == 500)      /*  20 sps */
	{
	  decimate(&FIR5,5,x_in,hdr.npts,x_in,&hdr.npts);
	  decimate(&FIR4,4,x_in,hdr.npts,x_in ,&hdr.npts);
	  hdr.delta = 1.0;
	}
      else if ((int)(hdr.delta * 1000.+0.5) != (int)(1000.+0.5))
	{ 
	  fprintf(stderr, "WARNING: non uniform samp. period between sac files\n") ; 
	  fprintf(stderr, "     ...file : %s with dt = %e is rejected\n", scfil, hdr.delta)  ; 
	  flag++   ; 
	  continue ; 
	} 
      /* Windowing */
      tdiff   += P_tt - (double)hdr.b - eq.preevent - (double)SAFETY_DELAY ;  /* time for the 1st sample */
      sampstart = (int) (tdiff/((double)hdr.delta) + 0.5) ;
      if (sampstart >= 0 && sampstart < hdr.npts) /* Windowing -- Screening by distance */
	{ 
	  /* Set new sac header variables */ 
	  hdr.delta = (float) 1.    ; 
	  hdr.o     = (float) otime ; 
	  hdr.npts  = hdr.npts - sampstart ;                 /* nb of samples             */
	  hdr.b     = hdr.b + ((float)sampstart)*hdr.delta ; /* shift of the first sample */
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
	  strcat(scfil,".decim.scr.sac");
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
	}
      else
	fprintf(stderr,"Warning: (trim) rejected incomplete file : %s (%d -- %d)\n",scfil,sampstart,hdr.npts);
    }
  fclose(fd);
  if (flag > 1)
    {
      fprintf(stderr, "Warning: %d files have been rejected because of non uniform sampling period\n",flag);
      fprintf(stderr, "    ...: if to much files are rejected, please screen or decimate data files manually\n");
    }

  o_sacf = openfile_wt(o_sacs);
  savetree(root, o_sacf, &eq.dmin, &eq.dmax);
  fclose(o_sacf);

  /* Memory Freeing */
  free((void*)i_locs)   ;
  free((void*)i_mseeds) ;
  free((void*)o_sacs)   ;
  free((void*)o_sacdir) ;
  free((void*)x_in) ;
  free((void*) dv)  ;
  free((void*) tv)  ;  
  free((void*)eq.wp_win4);
  free((void*)stla);
  free((void*)stlo);
  free((void*)stel);
  for(i=0;i<nids;i++)
    free((void*)ids[i]);
  free((void**)ids);
  freetree(root);
  freetree(mod);
  return 0;
}


void 
makeid(sachdr *hdr, char *id)
{
  int tmp;
  strcpy(id, hdr->knetwk) ;
  id[nbchar(hdr->knetwk)] = '\0';
  strcat(id,"_") ;
  strncat(id, hdr->kstnm,nbchar(hdr->kstnm)) ;
  strcat(id,"_") ;
  
  tmp = nbchar(hdr->khole) ;
  if ( tmp == 0 )
    strcat(id, "--") ;
  else
    strncat(id,hdr->khole,tmp);
  
  strcat(id,"_") ;
  strncat(id, hdr->kcmpnm, nbchar(hdr->kcmpnm)) ;
  strcat(id,"");
}

int 
findid(char *id, char **ids, int n)
{
  int i;
  /* Find id */
  for (i=0; i<n; i++)
    {
      if (strcmp(id,ids[i])==0)
	  return i;
    }
  fprintf (stderr, "WARNING (rec_dec_filt): %s not found in id list (rejected)\n",id) ;
  return -1;
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
get_params(int argc, char **argv, int *un, char **i_locs, char **i_mseeds, char **o_sacdir, 
	   char **o_sacs, str_quake_params *eq)
{
  int  i              ;
  char **keys, *i_tmp ;  

  /* Check syntax    */
  if ( argc < 6 ) {
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
  (*i_mseeds) = char_alloc(FSIZE) ; 
  (*o_sacdir) = char_alloc(FSIZE) ; 
  (*o_sacs)   = char_alloc(FSIZE) ;
  eq->wp_win4 = double_alloc(4)   ;
  eq->cmtfile[0] = '\0';

  /* Get args        */
  strcpy(     i_tmp , argv[1]) ;
  strcpy((  *i_locs), argv[2]) ;
  strcpy((*i_mseeds), argv[3]) ;
  strcpy((*o_sacdir), argv[4]) ;
  strcpy((  *o_sacs), argv[5]) ;

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
ymd2jul(int yyyy,int mm,int dd, int *jul)
{
  int k ;
  int ndays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  if ( (yyyy%4)   == 0 ) ndays[1] ++ ;
  if ( (yyyy%100) == 0 ) ndays[1] -- ;
  if ( (yyyy%400) == 0 ) ndays[1] ++ ;
  (*jul) = dd ;
  for (k=0; k<mm-1; k++)
    (*jul) += ndays[k] ;
}

void 
jul2ymd(int yyyy,int jul,int *mm,int *dd)
{
  int k=0 ;
  int ndays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  if ( (yyyy%4)   == 0 ) ndays[1] ++ ;
  if ( (yyyy%100) == 0 ) ndays[1] -- ;
  if ( (yyyy%400) == 0 ) ndays[1] ++ ;
  while(jul>ndays[k])
    jul -= ndays[k++];
  (*mm) = k+1 ;
  (*dd) = jul ;
}

void
date2epoch(int year, int mm, int dd, int hour, int min, int sec, int msec, double *epoch)
{
  struct tm date;
  time_t     tmp;
  extern long timezone;
  date.tm_sec   = sec  ;
  date.tm_min   = min  ;
  date.tm_hour  = hour ;
  date.tm_mday  = dd   ;
  date.tm_mon   = mm-1 ;
  date.tm_year  = year - 1900 ;
  date.tm_isdst = 0    ;
  tzset();
  tmp           = mktime(&date) ;
  *epoch        = (double)tmp + (double)(msec)/1000. - timezone ;
}

void
delta_t(int y1, int j1, int h1, int m1, int s1, int ms1,
         int y0, int j0, int h0, int m0, int s0, int ms0, double *tdiff)
{
  double t1, t0 ;
  date2epoch(y1,1,1,h1,m1,s1,ms1,&t1) ;
  date2epoch(y0,1,1,h0,m0,s0,ms0,&t0) ;
  *tdiff = t1 - t0 + (j1 - j0)*86400    ;
}



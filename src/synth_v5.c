#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c  */
#include "rwsacs.h"       /* rwsacs.c       */
#include "rwtextfiles.h"  /* rwtextfiles.c  */
#include "read_i_files.h" /* read_i_files.c */
#include "butterworth.h"  /* butterworth_sin.c */
#include "syn_conv_sub.h" /* syn_conv_sub.c */

void get_params(char *file, str_quake_params *eq)      ;
void minimax(double *y, int np, double *ymin, double *ymax) ;


int 
main(int argc, char *argv[])
{
  int i, j, k, ns, flag, ierror, opt ;
  double *tmparray, *data, *x_conv, depmin, depmax, dt;  
  double *b1, *b2, *a1, *a2, gain; 
  char *i_master, *i_wpfilname, *datafile, *buf, *GF  ;
  char *gf_file, *o_dir, *o_file  ;
  char  gfdirs[6][7] = {"gf_rr/","gf_tt/", "gf_pp/",
			"gf_rt/", "gf_rp/","gf_tp/"} ;
  str_quake_params eq   ;
  sachdr hd_data, hd_GF ; 
  FILE *i_wp            ;

  /* Input params     */
  if (argc < 5)
    {
      fprintf(stderr,"Error input params \n");
      fprintf(stderr,"Syntax : %s i_master cmtfile i_wpinversion o_direct\n",
	      argv[0]);
      exit(1);
    }
  i_master     = char_alloc(FSIZE) ;
  i_wpfilname  = char_alloc(FSIZE) ;
  o_dir        = char_alloc(FSIZE) ;
  strcpy(   i_master, argv[1])     ;
  strcpy(i_wpfilname, argv[3])     ;
  strcpy(      o_dir, argv[4])     ;
  get_params(i_master, &eq)   ;
  strcpy( eq.cmtfile, argv[2])     ;

  /* Allocates memory */
  datafile = char_alloc(FSIZE) ;
  buf      = char_alloc(200) ;
  gf_file  = char_alloc(FSIZE) ;  

  eq.wp_win4  = double_alloc(4);
  eq.vm    = double_alloc2p(2) ;
  eq.vm[1] = double_calloc(6)  ;

  tmparray = double_alloc((int)__LEN_SIG__) ;  
  data     = double_alloc((int)__LEN_SIG__) ;  
  x_conv   = double_alloc((int)__LEN_SIG__) ;  

  hdr_alloc(&hd_data)     ;
  hdr_alloc(&hd_GF)       ;  

  get_cmtf(&eq,2) ;
  /* Reads list of data files */
  ierror = 1 ;
  opt = 0; /* *** TO BE MODIFIED *** */
  i_wp   = openfile_rt(i_wpfilname, &ns);
  for(i=0; i<ns; i++)
    {
      flag = fscanf (i_wp, "%s", datafile) ;
      fgets(buf,200,i_wp); /* end of line */
      check_scan(1, flag, i_wpfilname, i_wp)  ;
      rhdrsac(datafile,  &hd_data, &ierror)   ;

      /* Set GF filenames */
      for(j=0; j<6; j++)      /* GF  */
	{
	  strcpy(gf_file,eq.gf_dir) ;
	  GF = get_gf_filename(gfdirs[j], hd_data.kstnm, hd_data.knetwk, hd_data.kcmpnm, ".SAC") ;
	  strcat(gf_file,GF)         ;
	  free((void*) GF)           ;
	  rhdrsac( gf_file, &hd_GF, &ierror)          ;
	  if (hd_GF.npts > (int)__LEN_SIG__)
	    {
	      fprintf(stderr, "WARNING : too much samples in file %s (truncated)", gf_file) ;
	      hd_GF.npts = (int)__LEN_SIG__ ;
	    }
	  rdatsac( gf_file, &hd_GF, tmparray, &ierror) ;
	  if (j==0)
	    for(k=0 ; k<hd_GF.npts ; k++)
	      data[k]  = tmparray[k] * eq.vm[1][j] ;
	  else
	    for(k=0 ; k<hd_GF.npts ; k++)
	      data[k] += tmparray[k] * eq.vm[1][j] ;
	}
      /* Perform convolution */
      conv_by_stf(&eq.ts, &eq.hd, "l", &hd_GF, data, x_conv) ;
      minimax(x_conv, hd_GF.npts, &depmin, &depmax) ;
      hd_GF.depmin = (float)depmin ;
      hd_GF.depmax = (float)depmax ;
      
      /* Write output file 1 */
      o_file = get_gf_filename(o_dir, hd_data.kstnm, hd_data.knetwk, hd_data.kcmpnm, ".complete_synth.sac") ;
      printf("Writing sac file : %s\n",o_file) ;
      whdrsac(o_file, &hd_GF);
      wdatsac(o_file, &hd_GF, x_conv);
      free((void*)o_file) ;
      if (i == 0)
	{
	  dt = (double)hd_GF.delta;
	  butter_sos(&eq.flow, &eq.fhigh, &eq.filtorder, &dt, &b1, &b2, &a1, &a2, &gain) ;
	}

      /* Apply sos */ 
      filter_data(x_conv, &hd_GF.npts, &eq.filtorder, b1, b2, a1, a2, &gain, &opt) ;
      minimax(x_conv, hd_GF.npts, &depmin, &depmax) ;
      hd_GF.depmin = (float)depmin ;
      hd_GF.depmax = (float)depmax ;

      /* Write output file 2 */
      o_file = get_gf_filename(o_dir, hd_data.kstnm, hd_data.knetwk, hd_data.kcmpnm, ".complete_synth.bp.sac") ;
      printf("Writing sac file : %s\n",o_file) ;
      whdrsac(o_file, &hd_GF);
      wdatsac(o_file, &hd_GF, x_conv);
      free((void*)o_file) ;
    }
  free((void*)tmparray);
  free((void*)data); 
  free((void*)x_conv);
  free((void*)b1);
  free((void*)b2);
  free((void*)a1);
  free((void*)a2);
  free((void*)i_master);
  free((void*)i_wpfilname);
  free((void*)datafile);
  free((void*)buf);
  free((void*)gf_file);
  free((void*)o_dir);
  return 0;
}



/************************************************/
/*           get_params(file, eq)          */
/************************************************/
/*  > Read input file for recursive filtering   */
void 
get_params(char *file, str_quake_params *eq)
{
  int  i     ;
  char **keys ;

  keys = char_alloc2(5, 16) ;
  i = 0 ;
  strcpy(keys[i++],"filt_order")  ;
  strcpy(keys[i++],"filt_cf1")    ;
  strcpy(keys[i++],"filt_cf2")    ;
  strcpy(keys[i++],"IDEC_2")      ;
  strcpy(keys[i++],"GFDIR")      ;
  
  eq->cmtfile[0] = '\0'           ;
  get_i_master(file, keys, 5, eq) ;

  for(i=0 ; i<5 ; i++)
    free((void*)keys[i]) ;
  free((void**) keys )   ;
}




void 
minimax(double *y, int np, double *ymin, double *ymax) 
{
  int i ;
  *ymin = y[0] ;
  *ymax = y[0] ;
  for (i=1 ; i<np ; i++)
    {
      if (y[i] > *ymax)
	*ymax = y[i] ;
      if (y[i] < *ymin)
	*ymin = y[i] ;
    }
}

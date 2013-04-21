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

/*     Centroid grid-search for LTZ components     */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <locale.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c  */
#include "rwsacs.h"       /* rwsacs.c       */
#include "rwtextfiles.h"  /* rwtextfiles.c  */
#include "travel_times.h" /* travel_times.c */
#include "read_i_files.h" /* read_i_files.c */
#include "butterworth.h" 
#include "syn_conv_sub.h" 
#include "wpinversion_LTZ.h"

void get_opt(int numarg1, int numarg2, char **argv, structopt *opt, str_quake_params *eq) ;
void get_param1(int argc, char **argv, int *M, structopt *opt, str_quake_params *eq);
void get_param2(char *file, structopt *opt, str_quake_params *eq) ;

int 
main(int argc, char *argv[])
{
  int i, j, nsac, M, ierror ;
  int    nh = NDEPTHS, nd = NDISTAS ;
  double s1, d1, r1, s2, d2, r2, gap, Cond ;
  double latopt,lonopt,depopt,tsopt,rmsopt,M0,Mw,diplow,M0_12 ;
  double *eval3, *global_rms, *data_norm, **TM, sdrM0[4] ;
  double **data,  **rms, ***G = NULL, ***dcalc ;
  double *dv, *tv;
  char **sacfiles ;
  FILE *o_log, *o_tmp ;
  structopt opt    ;
  sachdr *hd_synt  ;
  str_quake_params eq ;  
  /* Allocate memory for input parameters */
  eq.wp_win4  = double_alloc(4) ;
  eq.vm    = double_alloc2p(2)  ;
  eq.vm[0] = double_calloc(NM)  ; 
  eq.vm[1] = double_calloc(NM)  ;
  TM       = double_alloc2(3,3) ;
  eval3    = double_alloc(3)    ;
  dv       = double_alloc(nd)   ;
  tv       = double_alloc(nd)   ;
  /* Initialize some variables */
  data       = NULL ; G       = NULL ; opt.wgt = NULL ;
  opt.rms_in = NULL ; opt.p2p = NULL ; opt.avg = NULL ;
  sacfiles   = NULL ; o_tmp = NULL ; 
  ierror     = 1;
  /* Get input parameters */
  get_param1(argc,argv,&M,&opt,&eq) ;
  get_param2(opt.i_master,&opt,&eq) ;  
  trav_time_init(nh,nd,eq.pde_evdp,dv,tv,&ierror) ;  
  /* Write log header     */
  o_log = openfile_wt(opt.log)                     ;
  if (opt.dts_step <= 0. && opt.xy_dx <= 0.)
    o_tmp = o_log ;
  w_log_header(argv,&opt,&eq,eq.wp_win4,o_log) ;
  /* Set the data vector  */
  set_data_vector(nd,dv,tv,&nsac,&data,&sacfiles,&hd_synt,&eq,&opt,o_log);
  /* RMS per channel*/
  data_norm = double_alloc(nsac);
  calc_data_norm(data,hd_synt,nsac,data_norm);
  /* Compute G */
  G    = double_alloc3p(nsac) ; /* Memory allocation for G */
  for(i=0;i<nsac;i++)
    G[i] = double_alloc2(NM,hd_synt[i].npts) ; 
  calc_kernel(&eq,&opt,hd_synt,nsac,"l",nd,dv,tv,G,o_log);
  /* Inversion       */
  global_rms = double_calloc(2*(opt.ref_flag+1)) ;
  if (opt.dc_flag) /* Double Couple inversion                 */
    {              /* Warning: This has not been fully tested */
      inversion(M,nsac,hd_synt,G,data,eq.vm[0],&Cond,&opt,NULL) ;
      get_planes(eq.vm[0],TM,eval3,&sdrM0[0],&sdrM0[1],&sdrM0[2],&s2,&d2,&r2) ;
      sdrM0[3] = (fabs(eval3[0]) + fabs(eval3[2])) / 2.  ;
      for(i=0;i<opt.ip;i++)
	sdrM0[opt.ib[i]-1] = opt.priorsdrM0[opt.ib[i]-1] ;
      for(i=0;i<4;i++)
	opt.priorsdrM0[i]=sdrM0[i];
      fprintf(stderr,"WARNING: **** Double couple inversion have not been fully tested yet ****\n");
      inversion_dc(nsac,hd_synt,G,data,sdrM0,global_rms,&opt,NULL) ;
      sdr2mt(eq.vm[0],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2])         ;
    }
  else
    inversion(M,nsac,hd_synt,G,data,eq.vm[0],&Cond,&opt,o_log) ;
  /* Predicted data  */
  dcalc = double_alloc3p(nsac) ;
  calc_data(nsac,hd_synt,G,eq.vm,data,dcalc,&opt,o_tmp);
  /* Get RMS and Gap */  
  rms = double_calloc2(nsac, 2*(opt.ref_flag+1)) ;
  calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms,&opt) ;
  get_gap(hd_synt,nsac,&gap) ;
  if (opt.dts_step <= 0. && opt.xy_dx <= 0.)
    {
      /* Set stike/dip/rake */      
      get_planes(eq.vm[0],TM,eval3,&s1,&d1,&r1,&s2,&d2,&r2) ;
      write_cmtf(opt.o_cmtf, &eq, eq.vm[0]) ;
      /* Set Moment and Magnitude (Harvard Definition) */
      M0     = ((fabs(eval3[0]) + fabs(eval3[2])) * (double)POW) / 2. ;
      Mw     = (log10(M0) - 16.1) / 1.5 ;
      diplow = d2 ;
      if (d1 < d2) diplow = d1 ;
      M0_12  = M0 * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ;
      /* Output */
      w_o_saclst(nsac,sacfiles,hd_synt,rms,data_norm,&opt) ; 
      output_products(&opt,&eq,s1,d1,r1,s2,d2,r2,TM,eval3,M0,M0_12,Mw,
		      global_rms,gap,Cond,nsac,hd_synt,o_log) ;
    }
  /* Time-shift grid search */
  if (opt.dts_step > 0.)
    {      
      strcpy(opt.o_cmtf,"ts_WCMTSOLUTION")  ;
      strcpy(opt.psfile,"ts_p_wpinversion") ;
      strcpy(opt.p_data,"ts_fort.15")       ;
      strcpy(opt.o_saclst,"ts_o_wpinversion") ;
      strcpy(opt.o_covf,"ts_o_covariance")  ;
      if (opt.hdsafe)
	{
	  ts_gridsearch(nsac,M,nd,dv,tv,hd_synt,data,global_rms,
			&opt,&eq,&tsopt,&rmsopt,o_log);
	  eq.ts = tsopt  ;      
	  eq.hd = eq.ts ;
	}
      else
	{
	  fast_ts_gridsearch(nsac,M,nd,dv,tv,hd_synt,data,G,dcalc,rms,global_rms, 
			     &opt,&eq,&tsopt,&rmsopt,o_log);
	  opt.dts_val = 0.; /* Optimum solution */
	  eq.ts += tsopt  ;      
	  eq.hd = eq.ts ;
	}
      realloc_gridsearch(nsac,rms,global_rms,dcalc,opt.ref_flag+1) ;      
      calc_kernel(&eq,&opt,hd_synt,nsac,"l",nd,dv,tv,G,o_log) ;
      if (opt.dc_flag) /* Double Couple inversion                 */
	{              /* Warning: This has not been fully tested */
	  for(i=0;i<4;i++)
	    sdrM0[i] = opt.priorsdrM0[i];	  
	  inversion_dc(nsac,hd_synt,G,data,sdrM0,global_rms,&opt,o_log) ;
	  sdr2mt(eq.vm[0],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2])          ;
	}
      else
	inversion(M,nsac,hd_synt,G,data,eq.vm[0],&Cond,&opt,o_log) ;
      calc_data(nsac,hd_synt,G,eq.vm,data,dcalc,&opt,o_log) ;
      calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms,&opt) ;
      get_gap(hd_synt,nsac,&gap) ;
      /* Set stike/dip/rake */
      get_planes(eq.vm[0],TM,eval3,&s1,&d1,&r1,&s2,&d2,&r2) ;
      write_cmtf(opt.o_cmtf, &eq, eq.vm[0]) ;
      /* Set Moment and Magnitude (Harvard Definition) */
      M0     = ((fabs(eval3[0]) + fabs(eval3[2])) * (double)POW) / 2. ;
      Mw     = (log10(M0) - 16.1) / 1.5 ;
      diplow = d2 ;
      if (d1 < d2) diplow = d1 ;
      M0_12  = M0 * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ;
      /* Output */
      w_o_saclst(nsac,sacfiles,hd_synt,rms,data_norm,&opt) ;
      output_products(&opt,&eq,s1,d1,r1,s2,d2,r2,TM,eval3,M0,M0_12,Mw,
		      global_rms,gap,Cond,nsac,hd_synt,o_log) ;
    }
  /* Centroid position Grid-search */
  if (opt.xy_dx > 0. || opt.dz > 0.)
    {
      strcpy(opt.o_cmtf,"xy_WCMTSOLUTION")  ;
      strcpy(opt.psfile,"xy_p_wpinversion") ;
      strcpy(opt.p_data,"xy_fort.15")       ;
      strcpy(opt.o_covf,"xy_o_covariance")  ;
      strcpy(opt.o_saclst,"xy_o_wpinversion") ;
      xy_gridsearch(nsac,M,nd,dv,tv,hd_synt,data,G,dcalc,rms,global_rms,
		    &opt,&eq,&rmsopt,&latopt,&lonopt,&depopt,o_log) ;
      eq.evla = latopt ;
      eq.evlo = lonopt ;
      eq.evdp = depopt ;
      realloc_gridsearch(nsac, rms, global_rms, dcalc, opt.ref_flag+1) ;
      calc_kernel(&eq,&opt,hd_synt,nsac,"l",nd,dv,tv,G,o_log);
      if (opt.dc_flag) /* Double Couple inversion                 */
	{              /* Warning: This has not been fully tested */
	  for(i=0;i<4;i++)
	    sdrM0[i] = opt.priorsdrM0[i];	  
	  inversion_dc(nsac,hd_synt,G,data,sdrM0,global_rms,&opt,o_log) ;
	  sdr2mt(eq.vm[0],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2])         ;
	}
      else
	inversion(M,nsac,hd_synt,G,data,eq.vm[0],&Cond,&opt,o_log) ;
      calc_data(nsac,hd_synt,G,eq.vm,data,dcalc,&opt,o_log) ;
      calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms,&opt) ;
      get_gap(hd_synt,nsac,&gap) ;
      /* Set stike/dip/rake */
      get_planes(eq.vm[0],TM,eval3,&s1,&d1,&r1,&s2,&d2,&r2) ;
      write_cmtf(opt.o_cmtf, &eq, eq.vm[0]) ;
      /* Set Moment and Magnitude (Harvard Definition) */
      M0     = ((fabs(eval3[0]) + fabs(eval3[2])) * (double)POW) / 2. ;
      Mw     = (log10(M0) - 16.1) / 1.5 ;
      diplow = d2 ;
      if (d1 < d2) diplow = d1 ;
      M0_12  = M0 * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ;
      /* Output */
      w_o_saclst(nsac,sacfiles,hd_synt,rms,data_norm,&opt) ;
      output_products(&opt,&eq,s1,d1,r1,s2,d2,r2,TM,eval3,M0,M0_12,Mw,
		      global_rms,gap,Cond,nsac,hd_synt,o_log) ;
    }
  fclose(o_log);

  /* Memory Freeing */
  free((void*)eq.wp_win4) ;
  free((void*)eq.vm[0])   ;
  free((void*)eq.vm[1])   ;
  free((void**)eq.vm)     ;
  free((void*)eval3)      ;
  free((void*)global_rms) ;
  free((void*)data_norm)  ;
  free((void*)dv) ;
  free((void*)tv) ;
  for(i=0 ; i<nsac ; i++)
    {
      free((void*)data[i])     ;
      free((void*)rms[i] )     ;
      free((void*)sacfiles[i]) ;
      free_G(&G[i]);
      for(j=0 ; j<opt.ref_flag+1 ; j++)
	free((void*)dcalc[i][j]) ;
      free((void**)dcalc[i]);
    }
  free((void**)data)     ;
  free((void**)rms )     ;
  free((void**)sacfiles) ;
  free((void***)G)       ;
  free((void***)dcalc)   ;
  for(i=0 ; i<3 ; i++)
    free((void*)TM[i]) ;
  free((void**)TM)     ;
  free((void*)opt.rms_in) ;
  free((void*)opt.rms_r)  ;
  free((void*)opt.p2p)    ;
  free((void*)opt.avg)    ;
  free((void*)opt.wgt)    ;
  free((void*)hd_synt)    ;
  return 0; 
}


void get_param2(char *file, structopt *opt, str_quake_params *eq)
{
  int  i=0, nimas ;
  char **keys     ;
  if (strlen(eq->gf_dir) == 0)
    {
      nimas = 11;
      keys = char_alloc2(nimas, 16) ;
      strcpy(keys[i++],"GFDIR")   ;
    }
  else
    {
      nimas = 10;
      keys = char_alloc2(nimas, 16) ;
    }
  strcpy(keys[i++],"EVNAME")     ;
  strcpy(keys[i++],"CMTFILE")    ;
  strcpy(keys[i++],"WP_WIN")     ;
  strcpy(keys[i++],"DMIN")       ;
  strcpy(keys[i++],"DMAX")       ;
  strcpy(keys[i++],"filt_order") ;
  strcpy(keys[i++],"filt_cf1")   ;
  strcpy(keys[i++],"filt_cf2")   ;
  strcpy(keys[i++],"filt_pass")  ;
  strcpy(keys[i++],"IDEC_2")     ;

  get_i_master(file,keys,nimas,eq) ;  
  opt->ref_flag = get_cmtf(eq,opt->ref_flag+1)-1 ;
  for(i=0 ; i<nimas ; i++)
    free((void*)keys[i]) ;
  free((void**) keys )   ;
}  

void 
dispsynt(char **argv)
{
  fprintf(stderr,"Syntax: %s [-imas imaster(in)] [-ifil stalist(in)] [-ofil stalist(out)]\n",argv[0]);
  fprintf(stderr,"              [-log logfil(out)] [-icmtf cmtfil(in)] [-ocmtf cmtfil(out)] [-osyndir out_synt_dir] \n")   ;
  fprintf(stderr,"              [-pdata calc_dat_txtfil(out)] [-wpbm wp_bitmap(out)] [-ocovf o_covariance(out)] \n")        ;
  fprintf(stderr,"              [-ps ps_filename(out)] [-refbm ref_bitmap(out)] [-th rms_threshold(in)] [-cth cond_thre(in)] \n") ;
  fprintf(stderr,"              [-wz wgt] [-wl wgt] [-wt wgt] [-df damp_fac(in)] [-med] [-old] [-nont] [-noref] [-h (help)]\n")                                                 ;
}


void 
disphelp(char **argv,structopt *opt)
{
  fprintf(stderr,"WPHASE CENTROID GRID-SEARCH \n\n") ;
  dispsynt(argv) ;  
  fprintf(stderr,"\nAll parameters are optional :\n");
  
  fprintf(stderr,"\nTime-shift grid-search :\n");
  fprintf(stderr,"  -dts real_value         apply a time-shift to Green's functions (no time-shift)\n");
  fprintf(stderr,"  -ts tsmin dts tsmax     fast time-shift grid-search (no grid-search)\n") ;
  fprintf(stderr,"  -ts_Nit                 nb. of iteration for centroid location grid-search (%d)\n",opt->ts_Nit) ; 

  fprintf(stderr,"\nCentroid grid-search parameters:\n");
  fprintf(stderr,"  -xy_dx                  initial spacial sampling period in degree (%f)\n",opt->xy_dx) ;
  fprintf(stderr,"  -xy_Nx                  half_width = xy_Nx*xy_dx (%d)\n",opt->xy_Nx) ;
  fprintf(stderr,"  -xy_Nopt                nb. of optimal points (%d)\n",opt->xy_Nopt) ;
  fprintf(stderr,"  -Nit                    nb. of iterations (%d)\n",opt->xy_Nit) ;

  fprintf(stderr,"\nInput files: \n");
  fprintf(stderr,"  -imas imasterfile       imaster file (%s)\n",opt->i_master);
  fprintf(stderr,"  -ifil stalistfile       input sac file list (%s)\n",opt->i_saclst);
  fprintf(stderr,"  -icmtf cmtfile          input CMTSOLUTION file (specified in i_master)\n");
  fprintf(stderr,"  -old                    using output sta. list file as input sta. list file, and read\n") ;
  fprintf(stderr,"                                additional parameters from this file (don't read additional params)\n") ;   
  
  fprintf(stderr,"\nInput paths: \n");
  fprintf(stderr,"  -gfdir path             Green's function directory (./GF/)\n");
  
  fprintf(stderr,"\nOutput files: \n");  
  fprintf(stderr,"  -noref                  do not read the reference solution in the input CMTSOLUTION file (ref. sol. used)\n") ;
  fprintf(stderr,"  -ocovf o_covariance     covariance file (%s)\n",opt->o_covf) ;
  fprintf(stderr,"  -ofil stalistfile       output sac file list (%s)\n",opt->o_saclst) ;
  fprintf(stderr,"  -log  logfile           output log file (%s)\n",opt->log) ;
  fprintf(stderr,"  -wpbm wp_pgm            wphase solution pgm file (no pgmfile)\n") ;
  fprintf(stderr,"  -refbm ref_pgm          reference solution pgm file (no pgmfile)\n") ; 
  
  fprintf(stderr,"\nOutput paths: \n");
  fprintf(stderr,"  -osyndir out_synt_dir   output synthetic directory (%s)\n",opt->osacdir) ;
  
  fprintf(stderr,"\nInversion :\n");  
  fprintf(stderr,"  -nont                   no constraints on the moment tensor trace (zero trace)\n") ;
  fprintf(stderr,"  -cth real_value         set conditioning threshold (no conditioning)\n") ;
  fprintf(stderr,"  -df real_value          set the damping factor to use when the conditioning threshold\n");
  fprintf(stderr,"                                specified by -cth is reached (no conditioning)\n") ;
  
  fprintf(stderr,"\n data weighting and screening: \n");  
  fprintf(stderr,"  -wz real_value          weight for LHZ channels (%f)\n",opt->wZ);
  fprintf(stderr,"  -wl real_value          weight for LHL channels (%f)\n",opt->wL);
  fprintf(stderr,"  -wt real_value          weight for LHT channels (%f)\n",opt->wT);

  fprintf(stderr,"\n  -h, --help              display this help and exit\n\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>\n") ;
  fflush(stdout);
  exit(0);
}

void 
error_syntax(argv, comp)
     char **argv, *comp ;
{
  if (!strlen(comp))
    fprintf(stderr,"*** ERROR ***\n");
  else
    fprintf(stderr,"*** ERROR %s***\n",comp);
  dispsynt(argv) ;
  fflush(stderr);
  exit(1);
}

void 
get_char_arg(argv, j, i, numarg2, str)
     int   i, j, numarg2 ;
     char  **argv, *str  ;
{
  char *comp ;
  comp = char_alloc(32);
  sprintf(comp,": Missing argument (%s) ",argv[j]);
  if (i==numarg2)        error_syntax(argv,comp) ;
  if (argv[j+1][0]=='-') error_syntax(argv,comp) ;
  strcpy(str,argv[j+1]) ;
  free((void*)comp)     ;
}

void 
get_num_arg(char **argv, int j,int i,int numarg2,char *type,
	    void *value)
{
  char comp[64] ;
  
  sprintf(comp,": Missing argument (%s) ",argv[j]);
  if (i==numarg2) error_syntax(argv,comp) ;
  sprintf(comp,": Incompatible argument (%s) ",argv[j]);
  if (argv[j+1][0]=='-') 
    {
      if (strlen(argv[j+1]) == 1 ) error_syntax(argv,comp) ;
      if (!isdigit(argv[j+1][1]) && argv[j+1][0] != '.')  error_syntax(argv,comp) ; 
    }
  else
    if (!isdigit(argv[j+1][0]) && argv[j+1][0] != '.')  error_syntax(argv,comp) ;
  sscanf(argv[j+1],type,value) ;
}

void 
get_num_arg3(char **argv, int j, int i, int numarg2, double *value1, 
	     double *value2, double *value3)
{
  char comp[32] ;

  sprintf(comp,": Missing argument (%s) ",argv[j]);
  if (i==numarg2 || i+1==numarg2 || i+2==numarg2)
    error_syntax(argv,comp) ;
  sprintf(comp,": Incompatible argument (%s) ",argv[j]) ;

  if (argv[j+1][0]=='-') 
    {
      if (strlen(argv[j+1]) == 1 ) error_syntax(argv,comp) ;
      if (!isdigit(argv[j+1][1]) && argv[j+1][0] != '.')  error_syntax(argv,comp) ; 
    }
  else
    if (!isdigit(argv[j+1][0]) && argv[j+1][0] != '.')  error_syntax(argv,comp) ;
  sscanf(argv[j+1],"%lf",value1) ;

  if (argv[j+2][0]=='-') 
    {
      if (strlen(argv[j+2]) == 1 ) error_syntax(argv,comp) ;
      if (!isdigit(argv[j+2][1]))  error_syntax(argv,comp) ;
    }
  else
    if (!isdigit(argv[j+2][0]))  error_syntax(argv,comp) ;
  sscanf(argv[j+2],"%lf",value2) ;

  if (argv[j+3][0]=='-') 
    {
      if (strlen(argv[j+3]) == 1 ) error_syntax(argv,comp) ;
      if (!isdigit(argv[j+3][1]))  error_syntax(argv,comp) ;
    }
  else
    if (!isdigit(argv[j+3][0]))  error_syntax(argv,comp) ;
  sscanf(argv[j+3],"%lf",value3) ;
}

void 
get_opt(numarg1, numarg2, argv, opt, eq)
     int   numarg1, numarg2 ; 
     char  **argv   ;
     str_quake_params *eq      ;
     structopt *opt ; 
{
  int i, j, k;

  strcpy(opt->i_master,"i_master")       ;
  strcpy(   eq->gf_dir,"")               ;
  strcpy(opt->i_saclst,"o_wpinversion")  ;
  strcpy(opt->o_saclst,"gs_o_wpinversion") ;
  strcpy(opt->log,"gs_wpinversion.log")  ;
  strcpy(opt->o_covf,"gs_o_covariance")  ;
  strcpy(opt->o_cmtf,"gs_WCMTSOLUTION")  ;
  strcpy(opt->p_data,"gs_fort.15")       ;
  strcpy(opt->psfile,"gs_p_wpinversion") ;
  strcpy(opt->wpbmfile,"")  ;
  strcpy(opt->refbmfile,"") ;
  strcpy(opt->osacdir,"./") ;
  strcpy(opt->tsgsfile,"grid_search_ts_out") ;
  strcpy(opt->xygsfile,"grid_search_xy_out") ;
  eq->cmtfile[0] = '\0' ;
  opt->th_val    = 0.  ;
  opt->cth_val   = 0.  ;
  opt->df_val    = 0.  ;
  opt->med_val   = 0.  ;
  opt->op_pa     = 0.  ;
  opt->ref_flag  = 1   ;
  opt->ntr_val   = 1.  ;
  opt->wZ        = 1.  ;
  opt->wL        = 1.  ;
  opt->wT        = 1.  ;
  opt->azp       = 0.  ;
  opt->ps        = 1   ;
  opt->dts_val   = 0.  ;
  opt->ts_Nit    = 4   ;
  opt->dts_min   = 1.  ;
  opt->dts_step  = 4.  ;
  opt->dts_max   = 168.;
  opt->xy_Nit    = 2   ;
  opt->xy_dx     = 0.4 ;
  opt->xy_Nx     = 3   ;
  opt->xy_Nopt   = 5   ;
  opt->hdsafe    = 0   ;
  opt->mindep    = 3.5 ;
  opt->dz        = 0.  ;
  opt->dc_flag   = 0   ;
  opt->ip        = 0   ;
  opt->ib[0]     = 0   ;
  k = 0 ;
  for( i = 0; i<numarg2; i++ )
    {
      j = i + numarg1 + 1 ;
      if (!strncmp(argv[j],"-imas",5))
	{
	  get_char_arg(argv, j, i, numarg2, opt->i_master) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-ifil",5))
	{
	  get_char_arg(argv, j, i, numarg2, opt->i_saclst) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-ofil",5))
	{
	  get_char_arg(argv, j, i, numarg2, opt->o_saclst) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-log",4))
	{
	  get_char_arg(argv, j, i, numarg2, opt->log) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-ocovf",6))
	{
	  get_char_arg(argv, j, i, numarg2, opt->o_covf) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-ocmtf",6))
	{
	  get_char_arg(argv, j, i, numarg2, opt->o_cmtf) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-icmtf",6))
	{
	  get_char_arg(argv, j, i, numarg2, eq->cmtfile) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-gfdir",6))
	{
	  get_char_arg(argv, j, i, numarg2, eq->gf_dir) ;
	  add_slash(eq->gf_dir);
	  k+=2 ;
	}   
      else if (!strncmp(argv[j],"-osyndir",8))
	{
	  get_char_arg(argv, j, i, numarg2, opt->osacdir) ;
	  k+=2 ;
	} 
      else if (!strncmp(argv[j],"-wpbm",6))
	{
	  get_char_arg(argv, j, i, numarg2, opt->wpbmfile) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-refbm",6))
	{
	  get_char_arg(argv, j, i, numarg2, opt->refbmfile) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-cth",4))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf",(double*)&opt->cth_val) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-df",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf",(double*)&opt->df_val) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-ts_Nit",7)) 
	{
	  get_num_arg(argv, j, i, numarg2, "%d", (int*)&(opt->ts_Nit));
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-dts",4))
	{ 
	  get_num_arg(argv, j, i, numarg2,"%lf",(double*)&opt->dts_val); 
	  k+=2 ;
	} 
      else if (!strncmp(argv[j],"-ts",3)) 
	{ 
      	  get_num_arg3(argv, j, i, numarg2, &opt->dts_min, &opt->dts_step, &opt->dts_max) ;
      	  if (opt->dts_min > opt->dts_max)
      	    error_syntax(argv,"Incompatible arguments -- tsmin < tsmax is required (-ts)") ;
      	  else if (opt->dts_step <= 0)
      	    error_syntax(argv,"Incompatible arguments --  dts > 0 is required (-ts)") ;
      	  else if (opt->dts_min < 0)
      	    error_syntax(argv,"Negative time-shift -- tsmax > tsmin > 0 is required(-ts)") ;
      	  k+=4 ;
      	  i++  ;
      	}
      else if (!strncmp(argv[j],"-hdsafe",7))
	{ 
	  opt->hdsafe = 1;
	  k+=1 ;}       
      else if (!strncmp(argv[j],"-xy_Nit",7)) 
	{
	  get_num_arg(argv, j, i, numarg2,"%d", (int*)&opt->xy_Nit);
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-dx",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf",(double*)&opt->xy_dx) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-Nx",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%d",(int*)&opt->xy_Nx) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-dz",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf",(double*)&opt->dz) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-minz",5))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf",(double*)&opt->mindep) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-noxy",5)) 
	{
	  opt->xy_dx = 0. ;
	  k+=1 ;
	}
      else if (!strncmp(argv[j],"-nots",5)) 
	{
	  opt->dts_step = 0. ;
	  k+=1 ;
	}
      else if (!strncmp(argv[j],"-Nopt",5))
	{
	  get_num_arg(argv, j, i, numarg2,"%d",(int*)&opt->xy_Nopt) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-otsgsf",7))
	{
	  get_char_arg(argv, j, i, numarg2, opt->tsgsfile) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-oxygsf",7))
	{
	  get_char_arg(argv, j, i, numarg2, opt->xygsfile) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-wz",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf", (double*)&opt->wZ) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-wl",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf", (double*)&opt->wL) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-wt",3))
	{
	  get_num_arg(argv, j, i, numarg2,"%lf", (double*)&opt->wT) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-strike",7))
	{
	  opt->dc_flag = 1;
	  get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+0) ;	
	  opt->ib[opt->ip++] = 1;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-dip",4))
	{	
	  opt->dc_flag = 1;
	  get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+1) ;	
	  opt->ib[opt->ip++] = 2;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-rake",5))
	{
	  opt->dc_flag = 1;
	  get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+2) ;	
	  opt->ib[opt->ip++] = 3;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-mom",4))
	{
	  opt->dc_flag = 1;
	  get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+3) ;       
	  opt->priorsdrM0[3] /= POW;
	  opt->ib[opt->ip++] = 4;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-dc",3))
	{
	  opt->dc_flag = 1 ;
	  k++ ;
	}
      else if (!strncmp(argv[j],"-old",4))
	{
	  opt->op_pa   = 1. ;	
	  k++ ;
	}
      else if (!strncmp(argv[j],"-nont",5))
	{
	  opt->ntr_val = 0. ;
	  k++ ;
	}
      else if (!strncmp(argv[j],"-noref",6))
	{
	  opt->ref_flag = 0 ;
	  k++ ;
	}
      else if (!strncmp(argv[j],"-nops",5))
	{
	  opt->ps = 0 ;
	  k++ ;
	}
      else if (!strncmp(argv[j],"-h",2))
	disphelp(argv,opt) ;
      else if (!strncmp(argv[j],"--help",6))
	disphelp(argv,opt) ;
    }
  if (k != numarg2)
    error_syntax(argv,"") ;
  add_slash(opt->osacdir) ;
}

void 
get_param1(int argc, char **argv, int *M, structopt *opt, str_quake_params *eq)
{
  int numarg1, numarg2 ;
  int max = 55 ;

  numarg1 = 0              ;
  numarg2 = argc-numarg1-1 ;

  if ( (argc < numarg1+1 ) || (numarg2 > max))
    error_syntax(argv," (nb of arguments) ");
  get_opt(numarg1, numarg2, argv, opt, eq);
  *M    = NM  ;
  if (opt->ntr_val > 0.)
    *M = NM-1 ;
}

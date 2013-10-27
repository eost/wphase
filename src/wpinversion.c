/***************************************************************************
*
*	              W phase source inversion package 	            
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Universite de Strasbourg / CNRS 
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

/*     W phase inversion   */

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
#include "wpinversion.h"


/* *** EXTERNAL FUNCTIONS *** */
void jacobi(double **a,int n, int np, double *d, double **v, int *nrot) ;
void eigsrt(double *d, double **v, int n) ;

void get_opt(int numarg1, int numarg2, char **argv, structopt *opt, str_quake_params *eq) ;
void get_param1(int argc, char **argv, int *M, structopt *opt, str_quake_params *eq);
void get_param2(char *file, structopt *opt, str_quake_params *eq) ;


int 
main(int argc, char *argv[])
{
  int i,j,nsac,M,nsini ;
  double s1,d1,r1,s2,d2,r2,gap,Cond ;
  double M0,Mw,diplow,M0_12,**TM ;
  double *eval3,*global_rms,*data_norm ;
  double **data,**rms,***G = NULL,***dcalc ;  
  double sdrM0[4];
  char **sacfiles ;
  FILE *o_log  ;
  structopt opt       ;
  sachdr    *hd_synt  ;
  str_quake_params eq ;    
  /* Allocate memory for input parameters */
  eq.wp_win4  = double_alloc(4) ;
  eq.vm    = double_alloc2p(2)  ;
  eq.vm[0] = double_calloc(NM)  ; 
  eq.vm[1] = double_calloc(NM)  ;
  TM       = double_alloc2(3,3) ;
  eval3    = double_alloc(3)    ;
  /* Initialize some pointers */
  data       = NULL ; G       = NULL ; opt.wgt = NULL ;
  opt.rms_in = NULL ; opt.p2p = NULL ; opt.avg = NULL ;
  sacfiles   = NULL ; opt.rms_r = NULL ;
  /* Get input parameters */
  get_param1(argc, argv, &M, &opt, &eq) ;
  get_param2(opt.i_master, &opt, &eq)    ;
  fflush(stdout) ;  
  /* Write log header     */
  o_log = openfile_wt(opt.log)                     ;
  w_log_header(argv, &opt, &eq, eq.wp_win4, o_log) ;
  /* Set G and data       */
  set_matrices (&nsac,&nsini,&sacfiles,&hd_synt,&data,&G,&opt,&eq,o_log) ; 
  fflush(stdout);
  /* Screening            */
  if (opt.med_val > 0.) 
    {
      median(nsac, &opt)   ;
      screen_med(&nsac, sacfiles, data, G, hd_synt, &opt, o_log) ; 
    }
  if (opt.th_val > 0.)
    screen_rms(&nsac, sacfiles, data, G, hd_synt, &opt, o_log) ;
  if (opt.rms_r_th > 0.)
    screen_ratio(&nsac,sacfiles,data,G,hd_synt,&opt,o_log) ;
  if (nsac < 1) 
    {
      fprintf(stderr,"Error : Too few accepted channels\n")     ;
      fprintf(stderr,"....... Clean manually the input file\n") ;
      fprintf( o_log,"Aborded: too few accepted channels\n")     ;
      fclose(o_log) ;
      exit(1)       ;
    }  
  printf("%4d accepted_channels (%d rejected)\n",nsac,nsini - nsac)          ;
  fprintf(o_log,"accepted_channels: %4d (%d rejected)\n",nsac,nsini - nsac)  ;
  fflush(o_log) ;
  /* RMS per channel*/
  data_norm = double_alloc(nsac) ;
  calc_data_norm(data,hd_synt,nsac,data_norm) ;
  /* Inversion       */
  global_rms = double_calloc(2*(opt.ref_flag+1)) ;  
  if (opt.dc_flag) /* Double Couple inversion                 */
    {              /* Warning: This has not been fully tested */
      inversion(M,nsac,hd_synt,G,data,eq.vm[0],&Cond,&opt,NULL)    ;
      get_planes(eq.vm[0],TM,eval3,&sdrM0[0],&sdrM0[1],&sdrM0[2],&s2,&d2,&r2) ;
      sdrM0[3] = (fabs(eval3[0]) + fabs(eval3[2])) / 2.            ;
      for(i=0;i<opt.ip;i++)
	sdrM0[opt.ib[i]-1] = opt.priorsdrM0[opt.ib[i]-1];
      fprintf(stderr,"WARNING: **** Double couple inversion have not been fully tested yet ****\n");
      inversion_dc(nsac,hd_synt,G,data,sdrM0,global_rms,&opt,o_log);
      sdr2mt(eq.vm[0],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2])         ;
    }
  else
    inversion(M,nsac,hd_synt,G,data,eq.vm[0],&Cond,&opt,o_log) ;
  
  /* Predicted data  */
  dcalc = double_alloc3p(nsac) ;
  calc_data(nsac,hd_synt,G,eq.vm,data,dcalc,&opt,o_log) ;
  /* Get RMS and Gap */
  rms = double_calloc2(nsac, 2*(opt.ref_flag+1))        ;
  calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms,&opt) ;
  get_gap(hd_synt, nsac, &gap) ;
  w_o_saclst(nsac,sacfiles,hd_synt,rms,data_norm,&opt)  ; 
  /* Set stike/dip/rake */
  get_planes(eq.vm[0], TM, eval3, &s1,&d1,&r1, &s2,&d2,&r2) ;
  write_cmtf(opt.o_cmtf, &eq, eq.vm[0]) ;  
  /* Set Moment and Magnitude (Harvard Definition) */
  M0     = ((fabs(eval3[0]) + fabs(eval3[2])) * (double)POW) / 2. ; 
  Mw     = (log10(M0) - 16.1) / 1.5 ;
  diplow = d2 ;
  if (d1 < d2) diplow = d1 ;
  M0_12  = M0 * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ;
  /* Output */
  output_products(&opt,&eq,s1,d1,r1,s2,d2,r2,TM,eval3,M0,M0_12,Mw,
		  global_rms,gap,Cond,nsac,hd_synt,o_log) ;
  fclose(o_log);
  /* Memory Freeing */
  free((void*)eq.wp_win4) ;
  free((void*)eq.vm[0])   ;
  free((void*)eq.vm[1])   ;
  free((void**)eq.vm)     ;
  free((void*)eval3)      ;
  free((void*)global_rms) ;
  free((void*)data_norm)  ;
  for(i=0 ; i<nsac ; i++)
    {
      free((void*)data[i])     ;
      free((void*)rms[i] )     ;
      free((void*)sacfiles[i]) ;
      free_G(&G[i]);
      for(j=0 ; j<(opt.ref_flag+1) ; j++)
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
      nimas = 10;
      keys = char_alloc2(nimas, 16) ;
      strcpy(keys[i++],"GFDIR")   ;
    }
  else
    {
      nimas = 9;
      keys = char_alloc2(nimas, 16) ;
    }
  strcpy(keys[i++],"EVNAME")    ;
  strcpy(keys[i++],"CMTFILE")   ;
  strcpy(keys[i++],"WP_WIN")    ;
  strcpy(keys[i++],"DMIN")      ;
  strcpy(keys[i++],"DMAX")      ;
  strcpy(keys[i++],"filt_order");
  strcpy(keys[i++],"filt_cf1")  ;
  strcpy(keys[i++],"filt_cf2")  ;
  strcpy(keys[i++],"filt_pass") ;

  get_i_master(file,keys,nimas,eq) ;
  opt->ref_flag = get_cmtf(eq,opt->ref_flag+1) - 1 ;
  
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
  fprintf(stderr,"WPHASE INVERSION \n\n") ;
  dispsynt(argv) ;  
  fprintf(stderr,"\nAll parameters are optional :\n");
  
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
  fprintf(stderr,"  -ocmtf cmtfil           output CMTSOLUTION file (%s)\n",opt->o_cmtf) ;
  fprintf(stderr,"  -ocovf o_covariance     covariance file (%s)\n",opt->o_covf) ;
  fprintf(stderr,"  -pdata calc_dat_txtfil  predicted data filename (%s)\n",opt->p_data) ;
  fprintf(stderr,"  -ofil stalistfile       output sac file list (%s)\n",opt->o_saclst) ;
  fprintf(stderr,"  -log  logfile           output log file (%s)\n",opt->log) ;
  fprintf(stderr,"  -ps ps_filename         ps filename (%s)\n",opt->psfile) ;
  fprintf(stderr,"  -nops                   no output ps file\n") ;
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
  fprintf(stderr,"  -wn real_value          weight for LHN channels (%f)\n",opt->wN);
  fprintf(stderr,"  -we real_value          weight for LHE channels (%f)\n",opt->wE);
  fprintf(stderr,"  -azp                    azimuth weighting (no az. ponderation)\n");
  fprintf(stderr,"  -med                    screening data before inversion (no pre-screening)\n") ;
  fprintf(stderr,"  -th real_value          reject data using a rms threshold (no rms threshold)\n") ; 
  
  fprintf(stderr,"\nTime-shift :\n");
  fprintf(stderr,"  -dts real_value         apply a time-shift to Green's functions (no time-shift)\n");
  fprintf(stderr,"\n  -h, --help              display this help and exit\n\nReport bugs to: <zacharie.duputel@unistra.fr>\n") ;
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
get_num_arg(argv,j,i,numarg2,type,value)
     int    i, j, numarg2 ;
     char   **argv, *type;
     double *value ;
{
  char *comp ;
  comp = char_alloc(32);

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
  free((void*)comp)             ;
}

void 
get_num_arg3(argv, j, i, numarg2, value1, value2, value3)
     int    i, j, numarg2 ;
     char   **argv ;
     double *value1, *value2, *value3;
{
  char *comp ;

  comp = char_alloc(32);

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
  free((void*)comp)              ;
}


void 
get_opt(numarg1, numarg2, argv, opt, eq)
     int   numarg1, numarg2 ; 
     char  **argv   ;
     str_quake_params *eq      ;
     structopt *opt ; 
{
  int i, j, k;

  strcpy(opt->i_master,"i_master")   ;
  strcpy(   eq->gf_dir,"")   ;
  strcpy(opt->i_saclst,"i_wpinversion")  ;
  strcpy(opt->o_saclst,"o_wpinversion")  ;
  strcpy(opt->log,"wpinversion.log")     ;
  strcpy(opt->o_covf,"o_covariance")     ;
  strcpy(opt->o_cmtf,"WCMTSOLUTION")     ;
  strcpy(opt->p_data,"fort.15")  ;
  strcpy(opt->psfile,"p_wpinversion.ps") ;
  strcpy(opt->wpbmfile,"")  ;
  strcpy(opt->refbmfile,"") ;
  strcpy(opt->osacdir,"./") ;
  strcpy(opt->tsgsfile,"grid_search_ts_out") ;
  eq->cmtfile[0] = '\0' ;
  opt->th_val    = 0. ;
  opt->rms_r_th  = 0. ;
  opt->cth_val   = 0. ;
  opt->df_val    = 0. ;
  opt->med_val   = 0. ;
  opt->op_pa     = 0. ;
  opt->ref_flag  = 1  ;
  opt->ntr_val   = 1. ;
  opt->wZ        = 1. ;
  opt->wN        = 1. ;
  opt->wE        = 1. ;
  opt->azp       = 0. ;
  opt->ps        = 1  ;
  opt->dts_val   = 0. ;
  opt->dts_min   = 0. ;
  opt->dts_max   = 0. ;
  opt->dts_step  = 0. ;
  opt->ts_Nit    = 1  ;
  opt->xy_dx     = 0. ;
  opt->xy_Nx     = 0  ;
  opt->xy_Nit    = 0  ;
  opt->xy_Nopt   = 0  ;
  opt->dc_flag   = 0  ;
  opt->ip        = 0  ;
  opt->ib[0]     = 0  ;
  opt->ncom      = 0  ;
  k = 0 ;
  for( i = 0; i<numarg2; i++ )
    {
      j = i + numarg1 + 1 ;
      fflush(stdout);
      if (!strncmp(argv[j],"-imas",5)){
		get_char_arg(argv,j,i,numarg2,opt->i_master) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-ifil",5)){
		get_char_arg(argv,j,i,numarg2,opt->i_saclst) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-ofil",5)){
		get_char_arg(argv,j,i,numarg2,opt->o_saclst) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-log",4)){
		get_char_arg(argv,j,i,numarg2,opt->log) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-ocovf",6)){
		get_char_arg(argv,j,i,numarg2,opt->o_covf) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-ocmtf",6)){
		get_char_arg(argv,j,i,numarg2,opt->o_cmtf) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-pdata",6)){
		get_char_arg(argv,j,i,numarg2,opt->p_data) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-icmtf",6)){
		get_char_arg(argv,j,i,numarg2,eq->cmtfile) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-gfdir",6)){
		get_char_arg(argv,j,i,numarg2,eq->gf_dir) ;
		add_slash(eq->gf_dir);
		k+=2 ;}   
      if (!strncmp(argv[j],"-comment",8)){
		if (opt->ncom >= NCOM)
		  fprintf(stderr,"WARNING: Too many comments: Comment %s ignored.",argv[j+1]);
		else
		  get_char_arg(argv,j,i,numarg2,opt->comments[opt->ncom++]) ;
		k+=2 ;}   
      if (!strncmp(argv[j],"-osyndir",8)){
		get_char_arg(argv,j,i,numarg2,opt->osacdir) ;
		k+=2 ;}      
      if (!strncmp(argv[j],"-ps",6)){
		get_char_arg(argv,j,i,numarg2,opt->psfile) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-wpbm",6)){
		get_char_arg(argv,j,i,numarg2,opt->wpbmfile) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-refbm",6)){
		get_char_arg(argv,j,i,numarg2,opt->refbmfile) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-th",3)){
		get_num_arg(argv,j,i,numarg2,"%lf",&opt->th_val) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-nr",3)){
		get_num_arg(argv,j,i,numarg2,"%lf",(double*)&opt->rms_r_th) ;
		if (opt->rms_r_th < 1.)
		  opt->rms_r_th = 1./opt->rms_r_th ;
		k+=2 ;}
      if (!strncmp(argv[j],"-cth",4)){
		get_num_arg(argv,j,i,numarg2,"%lf",&opt->cth_val) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-df",3)){
		get_num_arg(argv,j,i,numarg2,"%lf",&opt->df_val) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-dts",4)){
		get_num_arg(argv,j,i,numarg2,"%lf",&opt->dts_val);
		k+=2 ;}
      if (!strncmp(argv[j],"-ogsf",6)){
		get_char_arg(argv,j,i,numarg2,opt->tsgsfile) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-wz",3)){
		get_num_arg(argv,j,i,numarg2,"%lf", &opt->wZ) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-wn",3)){
		get_num_arg(argv,j,i,numarg2,"%lf", &opt->wN) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-we",3)){
		get_num_arg(argv,j,i,numarg2,"%lf", &opt->wE) ;
		k+=2 ;}
      if (!strncmp(argv[j],"-strike",7)){
		opt->dc_flag = 1;
		get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+0) ;	
		opt->ib[opt->ip++] = 1;
		k+=2 ;}      
      if (!strncmp(argv[j],"-dip",4)){	
		opt->dc_flag = 1;
		get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+1) ;	
		opt->ib[opt->ip++] = 2;
		k+=2 ;}      
      if (!strncmp(argv[j],"-rake",5)){
		opt->dc_flag = 1;
		get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+2) ;	
		opt->ib[opt->ip++] = 3;
		k+=2 ;}      
      if (!strncmp(argv[j],"-mom",4)){
		opt->dc_flag = 1;
		get_num_arg(argv,j,i,numarg2,"%lf",(opt->priorsdrM0)+3) ;
		opt->priorsdrM0[3] /= POW;
		opt->ib[opt->ip++] = 4;
		k+=2 ;}      
      if (!strncmp(argv[j],"-dc",3)){
		opt->dc_flag = 1 ;
		k++ ;}      
      if (!strncmp(argv[j],"-azp",4)){
		opt->azp = 1 ;
		k+=1 ;}
      if (!strncmp(argv[j],"-med",4)){
		opt->med_val = 1 ;	
		k++ ;}
      if (!strncmp(argv[j],"-old",4)){
		opt->op_pa   = 1. ;	
		k++ ;}
      if (!strncmp(argv[j],"-nont",5)){
		opt->ntr_val = 0. ;
		k++ ;}
      if (!strncmp(argv[j],"-noref",6)){
		opt->ref_flag = 0 ;
		k++ ;}
      if (!strncmp(argv[j],"-nops",5)){
		opt->ps = 0 ;
		k++ ;}
      if (!strncmp(argv[j],"-h",2))
		disphelp(argv,opt) ;
      if (!strncmp(argv[j],"--help",6))
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
  int max = 128 ;

  numarg1 = 0              ;
  numarg2 = argc-numarg1-1 ;
  if ( (argc < numarg1+1 ) || (numarg2 > max))
    error_syntax(argv," (nb of arguments) ");
  get_opt(numarg1, numarg2,argv, opt, eq);
  *M    = NM  ;
  if (opt->ntr_val > 0.)
    *M = NM-1 ;
}

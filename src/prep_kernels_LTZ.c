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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/* Subroutines headers */
#include "rwsacs.h"
#include "proto_alloc.h"
#include "read_i_files.h"
#include "travel_times.h"
#include "butterworth.h" 
#include "syn_conv_sub.h" 

#define  YES  1
#define  NON  0

/* External routines */
int r_scr_dat_fil_list(char* file, char*** stats, char*** nets, float** stlats, 
		       float** stlons);
void taper_syn(double *R, double *T, double *P, sachdr *hdr);
void distaz(double cmt_lat, double cmt_lon, float* stlats, float* stlons, 
	    int nstat, float* dists, float* azs, float* bazs, float* xdegs,
	    long int* nerr);  
void rotate_traces(double *T, double *P, float baz, int npts, double *N, double *E);
void save_gf_sac(char *sac_filename, char *stnm, char *netwk, char *chan, 
		 float *lat, float *lon, sachdr *hdr,  double*   depval);
void fast_synth_sub(double az, double baz, double xdeg, double *tv, double *dv, int nd, 
		    str_quake_params *eq, sachdr *hdr, double **GFs, double *Z, double *TH, double *PH);
void crea_dir(const char *path);

/* Internal routines */
void get_params(int argc, char **argv, char *stat_file, char *itype, 
		char *i_master, str_quake_params *eq, int *tapering);

int 
main(int argc, char **argv)
{
  int i,j,k,flag,jstat,nstat,ngfcomp=6,nsects,ierror=1 ;
  int tapering=NON,nh=NDEPTHS,nd=NDISTAS ;
  long int nerr ; 
  char stat_file[FSIZE],i_master[FSIZE],path[FSIZE];
  char sacfile[FSIZE],itype[2],*chans[]={"LHZ","LHL","LHT"};
  char *gfcomp[]={"rr","tt","pp","rt","rp","tp"}; 
  char **stats, **nets ; 
  float *stlats,*stlons,*dists,*azs,*bazs,*xdegs,b ;
  double **GFs,*Z,*TH,*PH, *x_conv,**WAV,*tv,*dv   ;
  double *b1,*b2,*a1,*a2,gain,dt=1.; 
  str_quake_params eq; 
  sachdr   hdr; 

  /* Set input parameters */
  get_params(argc, argv, stat_file, itype, i_master, &eq, &tapering) ;
  get_cmtf(&eq, 1) ;
  /* Allocations */
  GFs = double_alloc2(10,__LEN_SIG__);/* GFs: Rrr, Rtt, Rpp, Rrt */
  Z   = double_alloc(__LEN_SIG__) ;   /* Vertical components     */
  TH  = double_alloc(__LEN_SIG__) ;   /* Radial components       */
  PH  = double_alloc(__LEN_SIG__) ;   /* Transverse components   */
  x_conv   = double_alloc((int)__LEN_SIG__) ;
  eq.vm    = double_alloc2p(2) ;
  eq.vm[1] = double_alloc(6)   ;
  eq.vm[0] = double_alloc(6)   ;   
  WAV      = double_alloc2p(3) ;
  *WAV     = Z  ;
  *(WAV+1) = TH ;
  *(WAV+2) = PH ;
  hdr_alloc(&hdr) ; /* SAC header allocation */
  nsects = (eq.flow > 0.)? eq.filtorder : eq.filtorder/2 ;
  b1 = double_alloc(nsects) ; 
  b2 = double_alloc(nsects) ;
  a1 = double_alloc(nsects) ; 
  a2 = double_alloc(nsects) ;
  tv = double_alloc(nd) ; /* travel times */
  dv = double_alloc(nd) ; /* distances    */
  /* Reading CMTSOLUTION and STAT_FILE */
  nstat = r_scr_dat_fil_list(stat_file, &stats, &nets, &stlats, &stlons) ;
  dists = float_calloc(nstat) ; 
  azs   = float_calloc(nstat) ;
  bazs  = float_calloc(nstat) ;
  xdegs = float_calloc(nstat) ;
  /* Distance, azimuth, back-azimuth, etc  */
  distaz(eq.evla, eq.evlo, stlats, stlons, nstat, dists, azs, bazs, xdegs, &nerr);
  /* Set travel time table for depth = dep */
  trav_time_init(nh,nd,eq.evdp,dv,tv,&ierror);
  /* Exitation kernels calculation */
  crea_dir(eq.gf_dir);
  flag = 0;
  for(i=0;i<ngfcomp;i++)
    {
      printf("**************************************\n"); 
      printf("Computing synthetics for M_%s...\n",gfcomp[i]);    
      strcpy(path,eq.gf_dir) ;
      strcat(path,"gf_")     ;
      strcat(path,gfcomp[i]) ;
      crea_dir(path)         ;
      strcat(path,"/")       ; 
      for(j=0;j<ngfcomp;j++)/* Inititializing the MT components */
	eq.vm[1][j] = 0. ;
      eq.vm[1][i]   = 1. ;
      for (jstat=0;jstat<nstat;jstat++) /* Computing exitation kernels for MT component #i at each station */
	{ 
	  /* Computing Z, TH, PH  */
	  printf("%-5s", stats[jstat]) ;
	  fast_synth_sub(azs[jstat],bazs[jstat],xdegs[jstat],tv,dv,nd,&eq,&hdr,GFs,Z,TH,PH);
	  if (tapering == YES) taper_syn(Z,TH,PH,&hdr);
	  b = hdr.b;
	  for(k=0;k<3;k++)
	    {
	      hdr.b = b                    ;
	      strcpy(sacfile,path)         ; /* Save Raw GF SAC */
	      strcat(sacfile,stats[jstat]) ;
	      strcat(sacfile,".")          ;
	      strcat(sacfile,nets[jstat])  ;
	      strcat(sacfile,".")          ;
	      strcat(sacfile,chans[k])     ;
	      strcat(sacfile,".SAC")       ;
	      save_gf_sac(sacfile,stats[jstat],nets[jstat],chans[k],&stlats[jstat],&stlons[jstat],&hdr,WAV[k]) ; 
	      conv_by_stf(eq.ts,eq.hd,itype,&hdr,WAV[k],x_conv) ;/* Perform convolution  */
	      if (flag == 0) /* Set the butterworth sos (dt must be the same for all stations)  */
		{
		  flag = 1 ; 
		  dt = (double)hdr.delta;
		  if (eq.flow>0.)
		    bpbu2sos(eq.flow,eq.fhigh,dt,eq.filtorder,&gain,b1,b2,a1,a2);
		  else
		    lpbu2sos(eq.fhigh,dt,eq.filtorder,&gain,b1,b2,a1,a2);		  
		}
	      else if (dt != (double)hdr.delta)
		{
		  fprintf(stderr, "ERROR: non uniform samp. period between sac files, file : %s\n",sacfile);
		  exit(1);
		}	  
	      strcat(sacfile,".sac") ; /* Save SAC after STF convolution   */
	      save_gf_sac(sacfile,stats[jstat],nets[jstat],chans[k],&stlats[jstat],&stlons[jstat],&hdr,x_conv) ; 
	      filter_with_sos(gain,b1,b2,a1,a2,nsects,x_conv,hdr.npts) ; /* Apply sos */
	      strcat(sacfile,".bp") ;  /* Save SAC after bandpass filtering */
	      save_gf_sac(sacfile,stats[jstat],nets[jstat],chans[k],&stlats[jstat],&stlons[jstat],&hdr,x_conv) ;
	    }
	}
      printf("\n");
    }
  /* Freeing memory */
  free((void*)b1) ;
  free((void*)b2) ;
  free((void*)a1) ;
  free((void*)a2) ;
  free((void*)Z) ; 
  free((void*)PH); 
  free((void*)TH); 
  free((void*)x_conv);
  for(i=0; i<10; i++)
    free((void *)GFs[i]) ;
  free((void**)GFs)      ;  
  free((void**)WAV)      ;  
  free((void*)eq.vm[0])  ;
  free((void*)eq.vm[1])  ;
  free((void**)eq.vm)    ;
  free((void*)dists)     ;
  free((void*)azs)       ;
  free((void*)bazs)      ;
  free((void*)xdegs)     ;
  for (i=0;i< nstat;i++)
    {
      free((void*)stats[i]) ;
      free((void*)nets[i])  ;
    }
  free((void**)stats) ;
  free((void**)nets)  ;
  free((void*)stlats) ;
  free((void*)stlons) ;
  if(tapering == YES) 
    {
      free((void*)tv);
      free((void*)dv);
    }
  printf("\n");
  return 0;
}


void 
save_gf_sac(char *sac_filename, char *stnm, char *netwk, char *chan, 
	    float *lat, float *lon, sachdr *hdr,  double*   depval)
{
  int  i,nbc;
  nbc = strlen(stnm);
  strncpy(hdr->kstnm, stnm, nbc);
  for (i=nbc; i<8; i++) 
    hdr->kstnm[i] = ' ';
  hdr->kstnm[8] = '\0';
  nbc = strlen(netwk);
  strncpy(hdr->knetwk,netwk,nbc);
  for (i=nbc; i<8; i++)
    hdr->knetwk[i] = ' ';
  hdr->knetwk[8] = '\0';  
  nbc = strlen(chan);
  strncpy(hdr->kcmpnm,chan, nbc);
  for (i=nbc; i<8; i++)
    hdr->kcmpnm[i] = ' ';
  hdr->kcmpnm[8] = '\0';
  hdr->stla = *lat;
  hdr->stlo = *lon;
  wsac(sac_filename, hdr, depval);
}


void
dispsynt(char **argv)
{
  fprintf (stderr, "Syntax: %s i_sacs stftype [-imas] [-gfdir] [-t] [-h]\n",argv[0]) ;
}

void
disphelp(char **argv)
{
  fprintf(stderr,"Convolution of Green's functions with source time function\n") ;
  dispsynt(argv) ;
  fprintf(stderr,"\n");
  fprintf(stderr,"  i_sacs            input sac files list\n");
  fprintf(stderr,"  stftype           STF type : 'g':gaussian\n")  ;
  fprintf(stderr,"                               'q':parabolic\n") ;
  fprintf(stderr,"                               'l':triangle\n")  ;
  fprintf(stderr,"                               'b':box car\n")   ;
  fprintf(stderr,"                               'c':cosine\n")    ;
  fprintf(stderr,"Optional parameters :\n");
  fprintf(stderr,"  -imas filename   i_master filename (i_master)\n");
  fprintf(stderr,"  -gfdir path      Green's function directory (if not used,\n");
  fprintf(stderr,"                     'GFDIR' must be specified in the i_master file)\n") ;
  fprintf(stderr,"  -t               Tapering between Origin and P_arrival times\n");
  fprintf(stderr,"  -h, --help       display this help and exit\n\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>\n") ;
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
get_num_arg(argv, j, i, numarg2, value)
     int    i, j, numarg2 ;
     char   **argv ;
     double *value ;
{
  char *comp ;
  comp = char_alloc(32);

  sprintf(comp,": Missing argument (%s) ",argv[j]);
  if (i==numarg2) error_syntax(argv,comp) ;
  sprintf(comp,": Incompatible argument (%s) ",argv[j]);
  if (argv[j+1][0]=='-') {
    if (strlen(argv[j+1]) == 1 ) error_syntax(argv,comp) ;
    if (!isdigit(argv[j+1][1]))  error_syntax(argv,comp) ; }
  else
    if (!isdigit(argv[j+1][0]))  error_syntax(argv,comp) ;
  sscanf(argv[j+1],"%lf",value) ;
  free((void*)comp)             ;
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
get_opt(int numarg1, int numarg2, char **argv, str_quake_params *eq, 
	char *i_master,int *tapering)
{
  int i,j,k;
  /* Default values */
  strcpy(i_master,"i_master")  ; 
  strcpy(eq->gf_dir,"")  ; 
  k = 0;
  for (i=0; i<numarg2; i++)
    {
      j = i+numarg1+1;
      if (!strncmp(argv[j],"-imas",5))
	{
	  get_char_arg(argv, j, i, numarg2, i_master) ;
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-gfdir",6))
	{
	  get_char_arg(argv, j, i, numarg2, eq->gf_dir) ;
	  add_slash(eq->gf_dir);
	  k+=2 ;
	}
      else if (!strncmp(argv[j],"-t",2))
	{
	  (*tapering) = YES ;
	  k++;
	}
      else if (!strncmp(argv[j],"-h",2))
	{
	  disphelp(argv)    ;
	  k++;
	}
      else if (!strncmp(argv[j],"--help",6))
	{
	  disphelp(argv)    ;
	  k++;
	}
    }
  if (k != numarg2)
    error_syntax(argv,"") ;
}

void
get_params(int argc, char **argv, char *stat_file, char *itype, 
	   char *i_master, str_quake_params *eq, int *tapering)
{
  int  numarg1, numarg2, i, nimas;
  char **keys ;

  numarg1 = 2              ;
  numarg2 = argc-numarg1-1 ;

  if(argc<numarg1+1)
    {
      if (!strncmp(argv[argc-1],"-h",2) || !strncmp(argv[argc-1],"-help",6))
	disphelp(argv);
      else
	error_syntax(argv," (nb of arguments) ");
    }
  
  /* Allocates filenames */
  strcpy(stat_file, argv[1]) ;
  strcpy(itype    , argv[2]) ;
  get_opt(numarg1,numarg2,argv,eq,i_master,tapering) ;

  /* Set keys to read in i_masterfile */
  i = 0 ;
  if (strlen(eq->gf_dir) == 0)
    {
      nimas = 6;
      keys = char_alloc2(nimas, 16) ;
      strcpy(keys[i++],"GFDIR")   ;
    }
  else
    {
      nimas = 5;
      keys = char_alloc2(nimas, 16) ;
    }
  strcpy(keys[i++],"CMTFILE")     ;  
  strcpy(keys[i++],"filt_order")  ;
  strcpy(keys[i++],"filt_cf1")    ;
  strcpy(keys[i++],"filt_cf2")    ;
  strcpy(keys[i++],"IDEC_2")      ;
  
  eq->cmtfile[0] = '\0'           ;
  get_i_master(i_master, keys, nimas, eq) ;
  printf("GFDIR: %s\n",eq->gf_dir); 
  for(i=0 ; i<nimas ; i++)
    free((void*)keys[i]) ;
  free((void**) keys )   ;
}


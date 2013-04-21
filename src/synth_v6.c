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

/*	Synthetic seismograms computation       */

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
#include "travel_times.h" /* travel_times.c */

void fast_synth_sub(double az, double baz, double xdeg, double *tv, double *dv, int nd, 
		    str_quake_params *eq, sachdr *hdr, double **GFs, double *Z, double *TH, double *PH);
void rotate_traces(double *T, double *P, float baz, int npts, double *N, double *E);
void distaz(double cmt_lat, double cmt_lon, float* stlats, float* stlons, 
	    int nstat, float* dists, float* azs, float* bazs, float* xdegs,
	    long int* nerr);  
void get_params(char *file, str_quake_params *eq)      ;

int 
main(int argc, char *argv[])
{
  int i, j, ns, flag, flagr, ierror, nsects, nh=NDEPTHS, nd=NDISTAS ;
  long int nerr ;
  double **GFs,*Z,*TH,*PH,*E,*N,*x_conv,**WAV ;
  double *b1, *b2, *a1, *a2, gain, dt=1.0     ; 
  double *tv, *dv ;
  float dist,az,baz,xdeg;
  char i_master[FSIZE], i_wpfilname[FSIZE], datafile[FSIZE], buf[200] ;
  char o_dir[FSIZE], *o_file,stnm[9],netwk[9],cmpnm[9];
  char stacmp[]={'Z','N','E','L','T'}  ;
  char itype[2]="l";
  str_quake_params eq ;
  sachdr hd_data, hd_synt ; 
  FILE *i_wp ;

  /* Input params     */
  if (argc < 5)
    {
      fprintf(stderr,"Error input params \n");
      fprintf(stderr,"Syntax : %s i_master cmtfile i_wpinversion o_direct [stftype]\n", argv[0]);
      fprintf(stderr,"stftype (optionnal) can be either:\n g (gaussian),\n q (parabolic),\n l (triangle,\n default),\n b(boxcar) or\n c (cosine)\n");	  
      exit(1);
    }
  strcpy(   i_master, argv[1]) ;
  strcpy(i_wpfilname, argv[3]) ;
  strcpy(      o_dir, argv[4]) ;
  get_params(i_master, &eq)    ;
  strcpy( eq.cmtfile, argv[2]) ;
  if (argc==6)
    {
      if (strlen(argv[5])==1)
	strcpy(itype,argv[5]);
      else
	{
	  fprintf(stderr,"Error input params \n");
	  fprintf(stderr,"Syntax : %s i_master cmtfile i_wpinversion o_direct [stftype]\n", argv[0]);
	  fprintf(stderr,"stftype (optionnal) can be either:\n g (gaussian),\n q (parabolic),\n l (triangle,\n default),\n b(boxcar) or\n c (cosine)\n");	  
	  exit(1);	  
	}
    }
  /* Allocates memory */
  eq.vm    = double_alloc2p(2) ;
  eq.vm[0] = double_calloc(6)  ;  
  eq.vm[1] = double_calloc(6)  ;
  GFs    = double_alloc2(10,__LEN_SIG__) ;/* GFs: Rrr, Rtt, Rpp, Rrt  */
  Z      = double_alloc(__LEN_SIG__) ;/*    Vertical components   */
  TH     = double_alloc(__LEN_SIG__) ;/*    Radial components     */
  PH     = double_alloc(__LEN_SIG__) ;/*    Transverse components */
  E      = double_alloc(__LEN_SIG__) ;/*    East components       */
  N      = double_alloc(__LEN_SIG__) ;/*    North components      */
  x_conv = double_alloc(__LEN_SIG__) ;
  WAV      = double_alloc2p(5) ;
  *WAV     = Z  ;
  *(WAV+1) = N  ;
  *(WAV+2) = E  ;  
  *(WAV+3) = TH ;  
  *(WAV+4) = PH ;  
  hdr_alloc(&hd_data) ;
  hdr_alloc(&hd_synt) ;  
  nsects = (eq.flow > 0.)? eq.filtorder : eq.filtorder/2 ;
  b1 = double_alloc(nsects) ; 
  b2 = double_alloc(nsects) ;
  a1 = double_alloc(nsects) ; 
  a2 = double_alloc(nsects) ;
  tv = double_alloc(nd); /* travel times */
  dv = double_alloc(nd); /* distances    */  
  /* Read CMTFILE */
  get_cmtf(&eq,2) ;    
  /* Set travel time table for depth = dep */
  ierror = 1 ;
  trav_time_init(nh,nd,eq.evdp,dv,tv,&ierror) ;
  /* Read list of data files */
  flag = 0   ;
  i_wp   = openfile_rt(i_wpfilname, &ns);
  for(i=0; i<ns; i++)
    {
      flagr = fscanf (i_wp, "%s", datafile) ;
      fgets(buf,200,i_wp); /* end of line */
      check_scan(1, flagr, i_wpfilname, i_wp)  ;
      rhdrsac(datafile,  &hd_data, &ierror)   ;
      /* Calculate azimuths, back-azimuths */
      dist = 0. ;
      az   = 0. ;
      baz  = 0. ;
      xdeg = 0. ;
      distaz(eq.evla,eq.evlo,&hd_data.stla,&hd_data.stlo,1,&dist,&az,&baz,&xdeg,&nerr) ;
      /* Computing Z, TH, PH  */	  
      fast_synth_sub(az,baz,xdeg,tv,dv,nd,&eq,&hd_synt,GFs,Z,TH,PH);
      rotate_traces(TH,PH,baz,hd_synt.npts,N,E) ; /* Rotating TH, PH to N, E */
      sscanf(hd_data.kstnm,"%s",stnm);
      sscanf(hd_data.knetwk,"%s",netwk);
      sscanf(hd_data.kcmpnm,"%s",cmpnm);
      for(j=0;j<5;j++)
	{
	  if (cmpnm[2] == stacmp[j])
	    break;
	}
      if (j==5)
	{
	  fprintf(stderr,"*** ERROR: Unknownk component %s for sta %s\n",cmpnm,stnm) ;
	  fprintf(stderr,"    -> Exiting\n") ;
	  fflush(stderr);
	  exit(1);
	}
      conv_by_stf(eq.ts,eq.hd,itype,&hd_synt,WAV[j],x_conv) ;/* Perform convolution */
      strcpy(hd_synt.kstnm,hd_data.kstnm)   ;
      strcpy(hd_synt.kcmpnm,hd_data.kcmpnm) ;
      strcpy(hd_synt.knetwk,hd_data.knetwk) ;
      hd_synt.stla = hd_data.stla ;
      hd_synt.stlo = hd_data.stlo ;
      hd_synt.evla = eq.pde_evla;
      hd_synt.evlo = eq.pde_evlo;
      hd_synt.evdp = eq.pde_evdp;
      /* Write output file 1 */
      o_file = get_gf_filename(o_dir,stnm,netwk,cmpnm[2],".complete_synth.sac") ;
      wsac(o_file,&hd_synt,x_conv);
      free((void*)o_file) ;
      if (flag == 0) /* Set the butterworth sos (dt must be the same for all stations)   */
	{
	  flag = 1 ; 
	  dt = (double)hd_data.delta;
	  if (eq.flow>0.)
	    bpbu2sos(eq.flow,eq.fhigh,dt,eq.filtorder,&gain,b1,b2,a1,a2);
	  else
	    lpbu2sos(eq.fhigh,dt,eq.filtorder,&gain,b1,b2,a1,a2);		  
	}
      else if ((int)(dt*1000+0.5) != (int)((double)hd_data.delta*1000+0.5))
	{
	  fprintf(stderr, "ERROR: non uniform samp. period between sac files, file : %s\n",datafile);
	  exit(1);
	}	  
      filter_with_sos(gain,b1,b2,a1,a2,nsects,x_conv,hd_synt.npts) ; /* Apply sos */
      /* Write output file 2 */
      o_file = get_gf_filename(o_dir,stnm,netwk,cmpnm[2],".complete_synth.bp.sac") ;
      printf("Writing sac file : %s\n",o_file) ;
      wsac(o_file,&hd_synt,x_conv);
      free((void*)o_file) ;
    }
  fclose(i_wp);
  free((void*)Z);
  free((void*)N);
  free((void*)E);
  free((void*)TH);
  free((void*)PH);
  free((void**)WAV);
  for(j=0;j<10;j++)
    free((void*)GFs[j]);
  free((void**)GFs);
  free((void*)x_conv);
  free((void*)b1);
  free((void*)b2);
  free((void*)a1);
  free((void*)a2);
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



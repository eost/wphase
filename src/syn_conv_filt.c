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

/*	STF convolution and butterworth filtering        */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c          */
#include "rwsacs.h"       /* rwsacs.c               */
#include "rwtextfiles.h"  /* rwtextfiles.c          */
#include "butterworth.h"  /* butterworth_sin.c */
#include "read_i_files.h" /* read_i_files.c */
#include "syn_conv_sub.h" /* syn_conv_sub.c */

/* internal functions */
void get_sacs(char *file, int *n, char ***sacs) ;
void get_params(int argc, char **argv, char **i_sacs, \
		char **itype, char **i_master, str_quake_params *eq);

int 
main(int argc, char *argv[])
{
  int    i, j, nbsta, ierror, count, nsects ;
  char   *i_sac_bp, *i_sacs, *path, *sacfile ;
  char   **sfiles, *itype ;
  double *b1, *b2, *a1, *a2, gain, dt=1.; 
  double *x_in, *x_conv;
  sachdr hdr ;
  str_quake_params eq ;
  int    ngfcomp  = 6 ;
  char   *gfcomp[]={"rr","tt","pp","rt","rp","tp"};  
  /* SET INPUT PARAMETERS */
  get_params(argc, argv, &i_sacs, &itype, &i_sac_bp, &eq) ;
  get_sacs(i_sacs, &nbsta, &sfiles) ;
  get_cmtf(&eq, 1) ;
  printf("       delay,  half_width\n");
  printf("%12.4f %12.4f\n", eq.ts, eq.hd);
  /* ALLOCATION */
  path = char_alloc(FSIZE) ;
  x_in   = double_alloc((int)__LEN_SIG__);
  x_conv = double_alloc((int)__LEN_SIG__);
  hdr_init(&hdr);  
  sacfile  = char_alloc(FSIZE) ;
  nsects = (eq.flow > 0.)? eq.filtorder : eq.filtorder/2 ;
  b1 = double_alloc(nsects) ; 
  b2 = double_alloc(nsects) ;
  a1 = double_alloc(nsects) ; 
  a2 = double_alloc(nsects) ;
  /* TREATING SYNTHETICS */
  count = 0;
  for ( i = 0; i<ngfcomp; i++)
    {
      printf("**************************************\n"); 
      printf("Treating synthetics for M_%s...\n",gfcomp[i]);
      strcpy(path,eq.gf_dir);
      strcat(path,"gf_") ;
      strcat(path,gfcomp[i]) ;
      strcat(path,"/") ;
      for (j = 0; j<nbsta; j++)
	{
	  strcpy( sacfile,      path);
	  strcat( sacfile, sfiles[j]);

	  /* Read data */
	  ierror = 0;
	  rhdrsac(sacfile, &hdr ,&ierror) ;
	  if (ierror)
	    continue;
	  rdatsac(sacfile, &hdr, x_in ,&ierror) ;
	  /* Perform convolution */
	  conv_by_stf(eq.ts,eq.hd,itype,&hdr,x_in,x_conv) ;
	  /* Write output-file 1 (to be removed???)*/
	  strcat( sacfile, ".sac");
	  wsac(sacfile, &hdr, x_conv);
	  /* Set the butterworth sos (samp. rate must be the same for all stations)*/
	  if (count == 0)
	    {
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
	  /* Apply sos */
	  filter_with_sos(gain,b1,b2,a1,a2,nsects,x_conv,hdr.npts) ; /* Apply sos */	  
	  /* Write output file 2 */
	  strcat( sacfile, ".bp");
	  printf("output: %s\n",sacfile);
	  wsac(sacfile, &hdr, x_conv);
	  count++;
	}
    }

  /* Memory freeing */
  free((void *)i_sac_bp);
  free((void *)i_sacs);
  free((void *)path);
  free((void *)sacfile);
  for(i=0;i<nbsta;i++)
    free((void *)sfiles[i]);
  free((void**)sfiles);
  free((void *)itype);
  free((void *)b1);
  free((void *)b2);
  free((void *)a1);
  free((void *)a2);
  free((void *)x_in);
  free((void *)x_conv);

  return 0; 
}

void 
get_sacs(char *file, int *n, char ***sacs)
{
  int   i, tmp ;
  FILE *saclist ;  

  /* Open file */
  saclist = openfile_rt(file,n);  
  
  /* Allocates memory */
  (*sacs) = char_alloc2((*n),FSIZE);
  
  for(i=0 ; i<(*n) ; i++)
    {
      tmp = fscanf (saclist, "%s", (*sacs)[i]);
      check_scan(1,tmp,file,saclist);
    }
  fclose(saclist);
}

void
dispsynt(char **argv)
{
  fprintf (stderr, "Syntax: %s i_sacs stftype [-imas] [-gfdir] [-h]\n",argv[0]) ;
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
get_opt(int numarg1, int numarg2, char **argv, str_quake_params *eq, char **i_master)
{
  int i,j,k;

  /* Default values */
  *i_master = char_alloc(FSIZE) ;
  strcpy(*i_master,"i_master")  ; 
  strcpy(eq->gf_dir,"")  ; 
  k = 0;
  for (i=0; i<numarg2; i++)
    {
      j = i+numarg1+1;
      if (!strncmp(argv[j],"-imas",5))
	{
	  get_char_arg(argv, j, i, numarg2, *i_master) ;
	  k+=2 ;
	}
      if (!strncmp(argv[j],"-gfdir",6))
	{
	  get_char_arg(argv, j, i, numarg2, eq->gf_dir) ;
	  add_slash(eq->gf_dir);
	  k+=2 ;
	}
      if (!strncmp(argv[j],"-h",2))
	disphelp(argv) ;
      if (!strncmp(argv[j],"--help",6))
	disphelp(argv) ;
    }
  if (k != numarg2)
    error_syntax(argv,"") ;
}

void
get_params(int argc, char **argv, char **i_sacs, char **itype, char **i_master, str_quake_params *eq)
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
  *i_sacs = char_alloc(FSIZE) ;
  *itype  = char_alloc(2) ;  
  
  strcpy(*i_sacs, argv[1]) ;
  strcpy( *itype, argv[2]) ;
  
  get_opt(numarg1,numarg2,argv,eq,i_master) ;

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
  get_i_master(*i_master, keys, nimas, eq) ;
  printf("GFDIR: %s\n",eq->gf_dir); 
  for(i=0 ; i<nimas ; i++)
    free((void*)keys[i]) ;
  free((void**) keys )   ;
}


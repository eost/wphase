#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c          */
#include "rwsacs.h"       /* rwsacs.c               */
#include "rwtextfiles.h"  /* rwtextfiles.c          */
#include "read_i_files.h"


typedef struct
{
  int fad;
} structopt ;


/* Internal routines */
void get_params(int argc, char **argv, char **scr_lst, char **o_lst, char **o_dir,
		str_quake_params *eq, structopt *opt);

/* External routines */
void distaz(double    cmt_lat,  double    cmt_lon,  float*    stlats,  float*    stlons, 
	    int       nstat,    float*    dists,    float*    azs,     float*    bazs,   
	    float*    xdegs,    long int* nerr);  


int 
main(int argc, char *argv[])
{
  int    MAX, flag,ierror, nl ;
  long   int nerr             ; 
  float  az, xdeg, dum        ;
  double stla, stlo, stel     ;
  double *x_in1 ;
  char   *i_fil, **h_fil, *o_fil, *scr_lst, *o_lst ;
  char   *o_dir, *sta, *net, *cmp, *prev_sta       ;
  FILE   *istaf, *ostaf       ;
  sachdr hdr1;
  str_quake_params eq ;
  structopt opt ;

  /* Input params */
  get_params(argc, argv, &scr_lst, &o_lst, &o_dir, &eq,&opt) ;


  /* Allocate data tabs     */
  MAX   = (int) __LEN_SIG__ ;
  x_in1 = double_alloc(MAX) ;

  /* Allocate sac header    */
  hdr_alloc(&hdr1) ; 
  
  /* String Allocation      */
  i_fil = char_alloc(FSIZE) ;
  h_fil = char_alloc2(2,FSIZE) ;
  sta   = char_alloc(9) ;
  net   = char_alloc(9) ;
  cmp   = char_alloc(9) ;
  prev_sta  = char_alloc(80) ;

  /* Open Data File List    */
  istaf = openfile_rt( scr_lst, &nl);
  ostaf = openfile_wt( o_lst);

  flag   = 0 ;
  ierror = 1 ;
  while( (nl=fscanf (istaf, "%s %s %s %s %lf %lf %lf %f %f", i_fil, sta, net, cmp, &stla, &stlo, &stel, &az, &xdeg)) != EOF )
    {
      if (nl == 0)
	break;
      else
	check_scan(9, nl, scr_lst, istaf);
      
      o_fil = get_gf_filename(o_dir, sta, net, cmp, ".data.sac") ;
      rhdrsac(i_fil, &hdr1, &ierror)        ;
      rdatsac(i_fil, &hdr1, x_in1, &ierror) ;
      /* az, baz, xdeg are not writen in sac files since we use PDE for the W phase time window */
      distaz(eq.evla, eq.evlo, &hdr1.stla, &hdr1.stlo, 1, &dum, &az, &dum, &xdeg, &nerr) ; 
      fprintf(ostaf,"%-50s %-9s %-9s %-9s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
	      o_fil,sta,net,cmp,stla,stlo,stel,az,xdeg) ;
      whdrsac(  o_fil, &hdr1) ;
      wdatsac(  o_fil, &hdr1, x_in1) ;
      free((void*)o_fil) ;
    }
  fclose(istaf);
  fclose(ostaf);
  return 0;
}


void 
dispsynt(char **argv)
{
  fprintf (stderr, "Syntax : %s i_sac_list(input) o_sac_list(output) o_dir(output) [-icmtf] [-h]\n", argv[0]) ;
}

void 
disphelp(char **argv)
{
  fprintf(stderr,"Write Centroid Lat, Lon, Dep and ep. distances in sac headers\n") ;
  dispsynt(argv) ;
  fprintf(stderr,"\n  i_sac_list      input sac files list\n");
  fprintf(stderr,"  o_sac_list      output sac files list\n");
  fprintf(stderr,"  o_dir           output directory for sac files\n");
  fprintf(stderr,"optional parameters :\n");
  fprintf(stderr,"  -icmtf   cmtfile   cmtfile (CMTSOLUTION)\n") ;
  fprintf(stderr,"  -h, --help         display this help and exit\n\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>\n") ;
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
get_opt(int numarg1, int numarg2, char **argv, str_quake_params *eq, structopt *opt)
{
  int i,j,k;

  /* Default values */
  k = 0        ;
  strcpy(eq->cmtfile,"CMTSOLUTION");
  for (i=0; i<numarg2; i++)
    {
      j = i+numarg1+1;
      if (!strncmp(argv[j],"-icmtf",6))
	{
	  get_char_arg(argv, j, i, numarg2, eq->cmtfile) ;
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
get_params(int argc, char **argv, char **scr_lst, char **o_lst, char **o_dir, 
	   str_quake_params *eq, structopt *opt)
{
  int numarg1, numarg2 ;

  numarg1 = 3              ;
  numarg2 = argc-numarg1-1 ;

  if(argc<numarg1+1)
    {
      if (!strncmp(argv[argc-1],"-h",2) || !strncmp(argv[argc-1],"-help",6))
	disphelp(argv);
      else
	error_syntax(argv," (nb of arguments) ");
    }
      

  /* Allocates filenames */
  *scr_lst = char_alloc(FSIZE) ;
  *o_lst   = char_alloc(FSIZE) ;
  *o_dir   = char_alloc(FSIZE) ;
 
  strcpy(*scr_lst, argv[1])    ;
  strcpy(*o_lst, argv[2])      ;
  strcpy(*o_dir, argv[3])      ;

  get_opt(numarg1,numarg2,argv,eq,opt) ;
  get_cmtf(eq,1) ;
}

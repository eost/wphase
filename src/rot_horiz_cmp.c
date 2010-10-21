/****************************************************************
*	W phase package - Hozizontal components rotation
*                                           
*       History
*             2010  Original Coding
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

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


/* Internal routines */
void get_params(int argc, char **argv, char **scr_lst, char **o_lst, char **o_dir,
		str_quake_params *eq);

/* External routines */
void distaz(double    cmt_lat,  double    cmt_lon,  float*    stlats,  float*    stlons, 
	    int       nstat,    float*    dists,    float*    azs,     float*    bazs,   
	    float*    xdegs,    long int* nerr);  


int 
main(int argc, char *argv[])
{
  int    i,MAX, flag,ind, ierror, nl;
  long   int nerr                   ; 
  float  dum, az, baz=0, xdeg ;
  double stla, stlo, stel   ;
  double *x_in1, *x_in2, tmp, co, si ;
  char   *i_fil, **h_fil, *o_fil, *scr_lst, *o_lst ;
  char   *o_dir, *sta, *net, *cmp, *prev_sta       ;
  FILE   *istaf, *ostaf       ;
  sachdr hdr1, hdr2, *hdproto ;
  str_quake_params eq ;

  /* Input params */
  get_params(argc, argv, &scr_lst, &o_lst, &o_dir, &eq) ;


  /* Allocate data tabs     */
  MAX   = (int) __LEN_SIG__ ;
  x_in1 = double_alloc(MAX) ;
  x_in2 = double_alloc(MAX) ;

  /* Allocate sac header    */
  hdr_alloc(&hdr1) ; 
  hdr_alloc(&hdr2) ;
  
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
  dum   = -1 ;
  flag   = 0 ;
  ierror = 1 ;
  while( (nl=fscanf (istaf, "%s %s %s %s %lf %lf %lf %f %f", i_fil, sta, net, cmp, &stla, &stlo, &stel, &az, &xdeg)) != EOF )
    {
      if (nl == 0)
	break;
      else
	check_scan(9, nl, scr_lst, istaf);
      
      if (strcmp(cmp, "LHZ") == 0)  /* Vertical component */
	{
	  o_fil = get_gf_filename(o_dir, sta, net, "LHZ", ".data.sac") ;
	  rhdrsac(i_fil, &hdr1, &ierror)        ;
	  rdatsac(i_fil, &hdr1, x_in1, &ierror) ;
	  /* az, baz, xdeg are not writen in sac files since we use PDE for the W phase time window */
	  az   = 0. ;
	  baz  = 0. ;
	  xdeg = 0. ;
	  distaz(eq.evla, eq.evlo, &hdr1.stla, &hdr1.stlo, 1, &dum, &az, &baz, &xdeg, &nerr) ; 
	  fprintf(ostaf,"%-50s %-9s %-9s %-9s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
		  o_fil,sta,net,cmp,stla,stlo,stel,az,xdeg) ;
	  whdrsac(  o_fil, &hdr1) ;
	  wdatsac(  o_fil, &hdr1, x_in1) ;
	  free((void*)o_fil) ;
	  continue ;
	}
      else if ( !flag || strcmp(prev_sta, sta) != 0 )  /* 1st horizontal component */
	{
	  if (strcmp(cmp, "LHE") == 0)
	    ind = 0;
	  else if (strcmp(cmp, "LHN") == 0)
	    ind = 1;
	  else {
	    fprintf(stderr,"Warning (rot_horiz_cmp) : unknown component %s in file %s\n", cmp, i_fil);
	    continue ; }
	  if (flag != 0) 
	    fprintf(stderr,"Warning (rot_horiz_cmp): incomplete channels for sta %s to perform rotation (horizontal data rejected)\n", prev_sta) ;
	  flag = 1 ;
	  strcpy(h_fil[ind],i_fil) ;
	  strcpy(  prev_sta,  sta) ;
	  continue ;
	}
      else if (strcmp(prev_sta, sta) == 0)             /* 2nd horizontal component */
	{
	  if (strcmp(cmp, "LHE") == 0) 
	    {
	      if (ind == 0) {
		fprintf(stderr,"Warning (rot_horiz_cmp) : %s redundancy for sta %s\n", cmp, prev_sta);
		continue ;  }
	      ind = 0;
	    }
	  else if (strcmp(cmp, "LHN") == 0)
	    {
	      if (ind == 1) {
		fprintf(stderr,"Warning (rot_horiz_cmp) : %s redundancy for sta %s\n", cmp, prev_sta);
		continue ;  }
	      ind = 1;	      
	    }
	  else
	    {
	      fprintf(stderr,"Warning (rot_horiz_cmp) : unknown component %s in file %s\n", 
		      cmp, i_fil);
	      continue ;
	    }
	  flag = 0 ;
	  strcpy( h_fil[ind], i_fil) ;
	  strcpy(   prev_sta,   sta) ;
	}
      else
	fprintf(stderr, "WOOPS... unpredicted case\n");
      
      /* Read input sac file */
      rhdrsac(h_fil[0], &hdr1, &ierror) ;
      rhdrsac(h_fil[1], &hdr2, &ierror) ;
      /* az, baz, xdeg are not writen in sac files since we use PDE for the W phase time window */
      distaz(eq.evla, eq.evlo, &hdr1.stla, &hdr1.stlo, 1, &dum, &az, 
	     &baz, &xdeg, &nerr) ; 
      hdr2.dist  = hdr1.dist     ;	  
      hdr2.az    = hdr1.az       ;
      hdr2.baz   = hdr1.baz      ;
      hdr2.gcarc = hdr1.gcarc    ;
      hdproto = &hdr1 ;
      if (hdr1.npts > hdr2.npts) 
	hdproto = &hdr2 ;

      rdatsac(h_fil[0], &hdr1, x_in1, &ierror) ; 
      rdatsac(h_fil[1], &hdr2, x_in2, &ierror) ;
      /* Rotate E/N components to L/T */
      co = cos(M_PI*((double)baz)/180.);
      si = sin(M_PI*((double)baz)/180.);
      for(i=0; i<hdr1.npts; i++)
	{
	  tmp  = -co*x_in2[i] -si*x_in1[i]     ;
	  x_in2[i] = -si*x_in2[i] + co*x_in1[i] ;
	  x_in1[i] = tmp ;
	}      
      /* Writing Longitudinal cmp */
      strcpy(hdproto->kcmpnm, "LHL     ") ;
      o_fil = get_gf_filename(o_dir, sta, net, "LHL", ".data.sac") ;
      printf("%-9s L : %-f", hdproto->kstnm, hdproto->cmpaz) ;
      whdrsac( o_fil, hdproto) ;
      wdatsac( o_fil, hdproto, x_in1) ;
      fprintf(ostaf,"%-50s %-9s %-9s %-9s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
	      o_fil,sta,net,"LHL",stla,stlo,stel,az,xdeg);      
      free((void*)o_fil) ;

      /* Writing Transverse cmp */
      strcpy(hdproto->kcmpnm, "LHT     ") ;
      o_fil = get_gf_filename(o_dir, sta, net, "LHT", ".data.sac") ;
      printf(" T : %-f\n", hdproto->cmpaz) ;
      whdrsac( o_fil, hdproto) ;
      wdatsac( o_fil, hdproto, x_in2) ;
      fprintf(ostaf,"%-50s %-9s %-9s %-9s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
	      o_fil,sta,net,"LHT",stla,stlo,stel,az,xdeg);
      free((void*)o_fil) ;
    }
  fclose(istaf);
  fclose(ostaf);

  free((void*)x_in1) ;
  free((void*)x_in2) ;
  free((void*)i_fil) ;
  free((void*)h_fil) ;
  free((void*)  sta) ;
  free((void*)  net) ;
  free((void*)  cmp) ;
  free((void*)prev_sta);

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
  fprintf(stderr,"Rotation of horizontal components East, North => radial,transverse components\n") ;
  dispsynt(argv) ;
  fprintf(stderr,"\n  i_sac_list      input sac files list\n");
  fprintf(stderr,"  o_sac_list      output sac files list\n");
  fprintf(stderr,"  o_dir           output directory for sac files\n");
  fprintf(stderr,"optional parameters :\n");
  fprintf(stderr,"  -icmtf   cmtfile   use centroid location in cmtfile\n") ;
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
get_opt(int numarg1, int numarg2, char **argv, str_quake_params *eq)
{
  int i,j,k;

  /* Default values */
  strcpy(eq->cmtfile,"CMTSOLUTION") ;
  k = 0        ;
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
	   str_quake_params *eq)
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

  get_opt(numarg1,numarg2,argv,eq) ;
  printf("Using centroid location in %s : \n",eq->cmtfile);
  printf("      corresponding azim. and ep. distances \n");
  printf("      will be used and writen in sac headers \n");
  get_cmtf(eq,1) ;
}

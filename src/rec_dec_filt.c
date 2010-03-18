#define _GNU_SOURCE
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Subroutines headers */
#include "proto_alloc.h"  /* proto_alloc.c          */
#include "rwsacs.h"       /* rwsacs.c               */
#include "rwtextfiles.h"  /* rwtextfiles.c          */
#include "butterworth.h"  /* butterworth_sin.c */
#include "read_i_files.h" /* read_i_files.c */
#include "sort_tree.h"

/* internal functions */
void makeid(char *id, sachdr *hdr) ;
int  findid(char *id, char **ids, int n) ;
void trapz(double *x_in, int npts, double dt, double *x_int) ;
void get_id_dec(char *file, int *n, char ***id, double **a1, double **a2, double **a3) ;
void set_sos(double fl, double fh, double dt, int *order, double **b1, double **b2, double **a1, double **a2, double *gain);
void get_filt_params(char *file, str_quake_params *eq) ;


int 
main(int argc, char *argv[])
{
  int     i, n ,tmp, count, opt, ierror = 1;
  int     MAX, nsac;
  double  scale, gain;
  double  *c1, *c2, *c3, *b1, *b2, *a1, *a2 ;
  double  *x_in, *x_int, dt ;
  char    **ids, *x_int_fil ;
  char    *decfile, *i_master, *i_sacfile ;
  char    *id, *o_sacfile ;
  FILE    *i_sacf, *o_sacf;
  sachdr  hdr;
  str_quake_params eq ;
  struct tree *datfil ;
  /* Set input files */

  if(argc<5){
    fprintf (stderr, "*** ERROR (input parameters : 4 params needed (%d given))\n", argc-1) ;
    fprintf (stderr, "Syntax : %s i_coeffs(input) i_master i_scr_data_list(input) o_scr_data_list(output)\n", argv[0]) ;
    exit (1) ; }

  decfile   = char_alloc(FSIZE) ;
  i_sacfile = char_alloc(FSIZE) ;
  o_sacfile = char_alloc(FSIZE) ;
  i_master  = char_alloc(FSIZE) ;
  strcpy(  decfile, argv[1])    ;
  strcpy( i_master, argv[2])    ;
  strcpy(i_sacfile, argv[3])    ; 
  strcpy(o_sacfile, argv[4])    ; 

  /* Allocate data tabs */
  MAX   = (int) __LEN_SIG__ ;
  x_in  = double_alloc(MAX) ;
  x_int = double_alloc(MAX) ;

  /* Allocate tree and sac header */
  hdr_alloc(&hdr)      ;
  datfil = alloctree();
  datfil->d    = NULL ;
  datfil->g    = NULL ;

  /* Allocate sac filenames */
  x_int_fil        = char_alloc(FSIZE);
  
  /* Allocate id */
  id = char_alloc(IDSIZE);

  /* Set input parameters */
  /*  > if opt == 1 the two first samples of the
		    deconvolved signal are zero.
      > if opt == 0 the convolution is computed
		    as the remaining sos.         */
  opt = 1 ;
  get_filt_params(i_master, &eq)               ;
  get_id_dec(decfile, &n, &ids, &c1, &c2, &c3) ;

  /* Open Data File Lists */
  i_sacf = openfile_rt(i_sacfile,&nsac) ;
  o_sacf = openfile_wt(o_sacfile)       ;

  i     = 0;
  count = 0;
  datfil->occur = 1      ;
  datfil->npts  = -12345 ;
  while( ( tmp = fscanf(i_sacf,"%s %s %s %s %lf %lf %lf %lf %lf\n",
			datfil->file,datfil->sta,datfil->net,datfil->cmp,
			&datfil->stla,&datfil->stlo,&datfil->stel,
			&datfil->az,&datfil->xdeg) ) != EOF )
    {
      if (tmp == 0)
	break;
      else
	check_scan(9, tmp, i_sacfile, i_sacf);
     
      /* Read input sac file */
      rhdrsac(datfil->file, &hdr, &ierror) ;
      if (hdr.npts > MAX) hdr.npts = MAX ;	
      rdatsac(datfil->file, &hdr, x_in, &ierror) ;
      if (hdr.npts < 2)
	{
	  fprintf(stderr,"Warning (rec_dec_filt) : %s (rejected) npts < 2\n",datfil->file);
	  continue;
	}

      /* Set the butterworth sos (samp. rate must be the same for all stations)*/
      if (count == 0 && i == 0)	
	{
	  dt = (double)hdr.delta;
	  set_sos(eq.flow, eq.fhigh, dt, &eq.filtorder, &b1, &b2, &a1, &a2, &gain) ;
	}
      else if (dt != (double)hdr.delta)
	{
	  fprintf(stderr, "ERROR (rec_dec_filt): non uniform samp. period between sac files, file : %s\n",datfil->file);
	  fclose(i_sacf);
	  exit(1);
	}

      /* Set filenames */
      strcpy(x_int_fil ,datfil->file);
      strcat(x_int_fil ,".id");
      strcat(datfil->file,".dec.bp.int");

      /* Base line operations */
      if ( eq.idtr == 1 )  /*   detrend and taper, eq.preevent and eq.fend are taper constants */
	{                  /*   usually 0.1 and 0.1                                            */
	  dtrd( x_in, hdr.npts, hdr.npts) ;               /* rm mean and trend */
	  taper(x_in, hdr.npts, eq.preevent, eq.fend) ;   /* tapering          */
	}
      else if (eq.idtr == 2)                  /*   shift the base line, eq.preevent=duration over which the  */
        dtrd (x_in, hdr.npts, eq.preevent) ;  /*   base line is defined, eq.fend = dummy                     */
                                              /*   if eq.idtr=0 no detrend and taper                         */

      /* Trapezoidal numerical integration */
      trapz(x_in, hdr.npts, (double)hdr.delta, x_int);
      
      /* Write header values and data in the sac input file */
      whdrsac(x_int_fil, &hdr)        ;
      wdatsac(x_int_fil, &hdr, x_int) ;

      makeid(id,&hdr);
      i = findid(id,ids,n);
      if (i == -1)
	continue ;

      /* Set the first line of the sos (Deconvolution) */
      b1[0] = c2[i]/c3[i] ;
      b2[0] = c1[i]/c3[i] ;
      a1[0] = -1.0 ;
      a2[0] =  0.0 ;
      scale = gain * c3[i] ;
      
      /* Applying sos */
      filter_data(x_in, &hdr.npts, &eq.filtorder, b1, b2, a1, a2, &scale, &opt) ;
      
      /* Trapezoidal numerical integration */
      trapz(x_in, hdr.npts, (double)hdr.delta, x_int);
      
      /* Integration again */
      trapz(x_int, hdr.npts, (double)hdr.delta, x_in);
      
      /* Write output sac  */
      whdrsac(datfil->file, &hdr);
      wdatsac(datfil->file, &hdr, x_in);
      savetree(datfil, o_sacf, &eq.dmin, &eq.dmax);
      count++;
    }
  fclose(i_sacf) ;
  fclose(o_sacf) ;
  printf("%d channels remain after deconvolution and filtering (%d rejected)\n", count, nsac-count) ;
  
  /*fprintf(o_sacf,"%-50s %-9s %-9s %-9s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
	  datfil->file, datfil->sta, datfil->net, datfil->cmp, datfil->stla,
	  datfil->stlo, datfil->stel, datfil->az, datfil->xdeg);*/
  

  /* Memory Freeing */
  free((void *)decfile);
  free((void *)i_master);
  free((void *)i_sacfile);
  free((void *)o_sacfile);
  free((void *)x_in);
  free((void *)x_int);
  free((void *)x_int_fil);
  free((void *)a1);
  free((void *)a2);
  free((void *)b1);
  free((void *)b2);
  free((void *)id);
  for(i=0; i<n; i++)
    free((void *)ids[i]);
  free((void *)c1);
  free((void *)c2);
  free((void *)c3);
  freetree(datfil);
  return 0;
}


/************************************************/
/*           get_filt_params(file, eq)          */
/************************************************/
/*  > Read input file for recursive filtering   */
void 
get_filt_params(char *file, str_quake_params *eq)
{
  int  i = 0  ;
  char **keys ;

  keys = char_alloc2(6, 16) ;
  strcpy(keys[i++],"filt_order") ;
  strcpy(keys[i++],"filt_cf1")   ;
  strcpy(keys[i++],"filt_cf2")   ;
  strcpy(keys[i++],"IDEC_2")     ;
  strcpy(keys[i++],"DMIN")     ;
  strcpy(keys[i++],"DMAX")     ;

  get_i_master(file, keys, 6, eq) ;

  for(i=0 ; i<6 ; i++)
    free((void*)keys[i]) ;
  free((void**) keys )   ;
}

void 
makeid(char *id, sachdr *hdr)
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


void 
trapz(double *x_in, int npts, double dt, double *x_int)
{
  int     i;
  
  //x_int = (double *) malloc (npts * sizeof(double));
  x_int[0] = 0.0 ;
  for (i=1 ; i < npts ; i++)
      x_int[i] = x_int[i-1] + (x_in[i-1] + x_in[i])*dt/2.0 ;
}


void 
get_id_dec(char *file, int *n, char ***id, double **a1, double **a2, double **a3)
{
  int   i, tmp ;
  FILE *decin ;  

  /* Open Data File List */
  decin = openfile_rt(file,n);  

  /* Allocates memory */
  (*id) = char_alloc2((*n),IDSIZE);
  (*a1) = double_alloc((*n));
  (*a2) = double_alloc((*n));
  (*a3) = double_alloc((*n));
  
  for(i=0 ; i<(*n) ; i++)
    {
      tmp = fscanf (decin, "%s %lf %lf %lf", (*id)[i],(*a1)+i,(*a2)+i,(*a3)+i);
      check_scan(4,tmp,file,decin);
    }
  fclose(decin);
}



/************************************************************************/
/*             set_sos(fl, fh, dt, order, b1, b2, a1, a2)               */
/************************************************************************/
/* Set the second order section for butterworth                         */
/*  the first line of the sos (corresponding to the deconvolution is    */
/*  proto_allocd but remains to be filled ...                              */
/* Input : fl : low cutoff freq.                                        */
/*         fh : high cutoff freq.                                       */
/*         dt : sampling period                                         */
/*         order : pointer to filter order                              */
/*                                                                      */
/* Output : order : order of the sos (int *)                            */
/*          b1,b2,a1,a2 : coef. of the sos (double **)                  */
void  
set_sos(double fl, double fh, double dt, int *order, double **b1, double **b2, double **a1, double **a2, double *gain) 
{
  int i;
  double  *bb1, *bb2, *ab1, *ab2;
  butter_sos(&fl, &fh, order, &dt, &bb1, &bb2, &ab1, &ab2, gain) ;
  /* add deconvolution to sos */
  (*b1) = double_alloc((*order)+1);
  (*b2) = double_alloc((*order)+1);
  (*a1) = double_alloc((*order)+1);
  (*a2) = double_alloc((*order)+1);
  for ( i=0; i< (*order); i++)
    {
      (*b1)[i+1]=bb1[i] ;
      (*b2)[i+1]=bb2[i] ;
      (*a1)[i+1]=ab1[i] ;
      (*a2)[i+1]=ab2[i] ;
    }
  (*order)++;
}

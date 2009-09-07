#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "complex.h"
#include "proto_alloc.h"    /* proto_alloc.c  */
#include "rwtextfiles.h"    /* rwtextfiles.c  */
#include "read_i_files.h"   /* read_i_files.c */

void plzr2resp(char *pzfilename , double tolerance, double fl, double fh, int nf , double fref, double *gain_factor, double *fc, double *h, double *s2, int *ierror) ;
void lsqenp_ (int *nf, int *mp, int *mv, float *yl, float *x, float *b, int *ip, int *ib, int *idvt, int *icon, int *iquit, int *iprnt) ;
void pzfile2id(char *pzfile, char *id, int *ierror) ;
void get_params(char *file, str_quake_params *eq)   ;

int 
main(int argc, char *argv[])
{
  int    tmp, ierror                                     ;
  double f0, fref, g, h, s2, w0, dt, a1, a2, a3          ; 
  char   *net, *stat, *locid, *chan, *id, *line          ;
  char   *pzfilename, *pzfilelist,*outputfile, *irecfile ; 
  FILE   *pzlist, *out ;
  str_quake_params eq  ;

  /* Allocates memory */
  line = char_alloc(FSIZE) ;
  pzfilename = char_alloc(FSIZE) ;
  pzfilelist = char_alloc(FSIZE) ;
  outputfile = char_alloc(FSIZE) ;  
  irecfile   = char_alloc(FSIZE) ;  
  net   = char_alloc(5) ;
  stat  = char_alloc(5) ;
  locid = char_alloc(5) ;
  chan  = char_alloc(5) ;
  id    = char_alloc(IDSIZE);

  /* INPUT PARAMETERS */
  
  dt = 1       ; /* *** TO BE MODIFIED *** */
  if(argc<4) {
    fprintf (stderr, "*** ERROR (input parameters\n 3 params needed (%d given)) ***",
	     argc-1)  ;
    fprintf (stderr, "Syntax: %s pzfilelist i_master coeff_list_file\n",
	     argv[0]) ;
    exit (1) ;}

  strcpy (pzfilelist, argv[1]) ;
  strcpy (irecfile, argv[2])   ;
  strcpy (outputfile, argv[3]) ;

  get_params(irecfile, &eq) ;
  fref = eq.fh ;  
 
  /* OPEN PZLIST FILE */
  pzlist = openfile_rt(pzfilelist, &tmp) ;

  /* OPEN OUTPUT FILE */
  out = openfile_wt(outputfile) ;

  while( (tmp=fscanf (pzlist, "%s", pzfilename)) != EOF ) 
    {
      if( tmp == 0 )
	break ;
      /* COMPUTE RESPONSE PARAMETERS */
      plzr2resp(pzfilename , eq.tol, eq.fl, eq.fh, eq.nf ,fref, &g, &f0, &h, &s2, &ierror) ;
      printf("%s %f %f %f %d\n",pzfilename,h, f0, g, ierror) ;
      if (ierror!=0)
	continue;
      pzfile2id(pzfilename, id, &ierror) ;
      if (!ierror)
	{
	  w0 = 2*M_PI*f0              ;
	  a1 = 1./dt                  ;
	  a2 = -2.*(a1+h*w0)/g        ;
	  a3 = (a1+2*h*w0+w0*w0*dt)/g ; 
	  a1 = a1/g                   ;
	  if (finite(a1) && finite(a2) && finite(a3))
	    fprintf (out, "%-14s%18.9e%18.9e%18.9e\n", id, a1, a2, a3) ; 
	  else {
	    fprintf(stderr,"WARNING (make_resp_lookup_table): Non-Finite coefficient detected\n")                 ;
	    fprintf(stderr,"          id : %s ; c1=%e ; c2=%e ; c3=%e\n", id, a1, a2, a3) ;}
	}
    }

  fclose (pzlist) ;
  fclose (out)    ;
  /* Memory Freeing */
  free((void *)line)       ;
  free((void *)pzfilename) ;
  free((void *)pzfilelist) ;
  free((void *)outputfile) ;  
  free((void *)irecfile)   ;  
  free((void *)net)   ;
  free((void *)stat)  ;
  free((void *)locid) ;
  free((void *)chan)  ;
  free((void *)id)    ;

  return 1 ;
}

void 
get_params(char *file, str_quake_params *eq) 
{
  int i=0     ;
  char **keys ;
  keys = char_alloc2(1, 16)  ;
  strcpy(keys[i++],"IDEC_3") ;

  get_i_master(file, keys, 1, eq) ;
  
  free((void*)keys[0]) ;
  free((void**) keys ) ;
}  


/**********************************************/
/*               readzeros(nz,zr,in)          */
/**********************************************/
/* Read zeros  from pzfile                    */
/* Input parameters  : in : pzfile descriptor */
/*                                            */
/* Output parameters : zr : zeros (nz values) */
/*                                            */
void 
readzeros(int *nz, complex *zr, FILE *in)
{
  int i, nzn0, nz0, tmp ;
  i=0 ; 
  while( (tmp=fscanf (in, "%lf%lf", &zr[i].real, &zr[i].imag)) != EOF ) 
    {
      if( tmp == 0 )
	break ;
      i++ ;
    }
  /* nzn0 is the number of zeros specified (often null)*/
  nzn0 = i ;
  /* nz0 is the number of zeros at 0 for the velocity response */
  nz0  = *nz - nzn0 ;
  if( nzn0 >= 1 )
    for( i=0; i<nzn0; i++ )
	zr[(*nz)-i-1] = zr[nzn0-i-1];
  for( i=0; i<nz0 ; i++ ) // on pourait ne travailler que sur 
    {                     // des pointeurs, mais il faudrais
      zr[i].real=0 ;      // déclarer un deuxième pointeur
      zr[i].imag=0 ;      // sur une struc complex
    }
}


/**********************************************/
/*               readpoles(np,pl,in)          */
/**********************************************/
/* Read poles from pzfile                     */
/* Input parameters  : in : pzfile descriptor */
/*                                            */
/* Output parameters : pl : poles (np values) */
/*                                            */
void 
readpoles(int *np, complex *pl, FILE *in, int *ierror) 
{
  int i;
  i=0 ;
  while(i < *np)
    {
      if(fscanf (in, "%lf%lf", &pl[i].real, &pl[i].imag) == 0 )
	{
	  fprintf (stderr, "ERROR reading poles of pzfiles\n") ;
	  fclose (in) ;
	  /*exit (1)  ; */
	  *ierror = 99;
	  return;
	}
      i++ ;
    }
}




/************************************************************/
/*           readpzfile(pzfilename,nz,zr,np,pl,sg)          */
/************************************************************/
/* Extract poles, zeros and the constant from pzfile        */
/* Input parameters  : pzfilename                           */
/*                                                          */
/* Output parameters : zr : zeros (nz complex values)       */
/*                     pl : poles (np complex values)       */
/*                     sg : constant value                  */
/*                                                          */
void 
readpzfile (char *pzfilename, int *nz, complex **zr, int *np, complex **pl, double *sg, int *ierror) 
{
  int  ns, zflag=0, pflag=0, cflag=0;
  char *line, *head ;
  FILE *in ;
  
  line = char_alloc(32) ;
  head = char_alloc(32) ;

  if ((in =fopen (pzfilename,"rt"))==NULL)
    {
      fprintf (stderr, "WARNING (make_resp_lookup_table): pzfile %s does not exist...\n     ======>this channel is skipped\n",pzfilename) ;
      /*exit (1) ;*/
      *ierror = 99;
      return;
    }

  
  while(fgets(line, 32, in) != NULL)
    {
      if(!strncmp(line,"ZEROS",5))/* We ONLY read the zeros following the first header "ZEROS" */
	{
	  if( !zflag )
	    {
	      ns = sscanf (line, "%s%d", head, nz) ;
	      *nz   = *nz-1; /* dispacement -> velocity (nzeros -> nzeros-1) */
	      printf("%d poles\n",(*nz)) ;
	      *zr   = complex_alloc(*nz) ;
	      zflag = 1;
	      readzeros(nz,*zr,in) ;
	    }
	  else
	    fprintf (stderr, "WARNING (make_resp_lookup_table): ZEROS redundancy in pzfile : %s\n",pzfilename) ;
	}
      if(!strncmp(line,"POLES",5))/* We ONLY read the poles following the first header "POLES" */
	{
	  if( !pflag )
	    {
	      ns    = sscanf (line, "%s%d", head, np) ;
	      *pl   = complex_alloc(*np) ;
	      pflag = 1 ;
	      readpoles(np,*pl,in,ierror) ;
	      if (*ierror==99)
		{
		  fclose (in);
		  return;
		}
	    }
	  else
	    fprintf (stderr, "WARNING (make_resp_lookup_table): POLES redundancy in pzfile : %s\n",pzfilename) ;
	}
      if(!strncmp(line,"CONSTANT",8))/* We ONLY read the constant following the first header "CONSTANT" */
	{
	  if( !cflag )
	    {
	      ns    = sscanf (line, "%s%lf", head, sg) ;
	      cflag = 1;	  
	    }
	  else
	    fprintf (stderr, "WARNING (make_resp_lookup_table): CONSTANT redundancy in pzfile : %s\n",pzfilename) ;
	}
    }


  if(!zflag || !pflag || !cflag)
    {
      fprintf (stderr, "WARNING (make_resp_lookup_table): Inclomplete pzfile : %s\n   ======>this channel is skipped \n", pzfilename) ;
      fclose (in);
      /*exit (1) ;*/
      *ierror = 99 ;
      return;
    }

  fclose (in);
  /* Memory Freeing */
  free((void *)line) ;
  free((void *)head) ;
}




/************************************************************/
/*               resp=pz2res (w,nz,zr,np,pl,sg)             */
/************************************************************/
/* Compute the response at w from poles, zeros and constant */
/* Input parameters  : w  : angular frequency               */
/*                     zr : zeros (nz complex values)       */
/*                     pl : poles (np complex values)       */
/*                     sg : constant value                  */
/* Output parameters : resp : (complex) response at w       */
/*                                                          */
complex 
pz2res (double w, int nz, complex *zr, int np, complex *pl, double sg, int *ierror)
{
  int i;
  complex num, den, resp;
  
  num = cmplx ( 1., 0. ) ;
  for ( i=0 ; i<nz ; i++ )
    {
      num = mul_c( num , sub_c(cmplx(0., w),zr[i]) ) ;
    }

  den = cmplx ( 1., 0. ) ;
  for ( i=0 ; i<np ; i++ )
    {
      den = mul_c( den , sub_c(cmplx(0., w),pl[i]) ) ;
    }
  
  resp = div_c( mul_c(num, cmplx(sg, 0)) , den , ierror) ;
  return (resp) ;
}





/*************************************************************/
/*              compresp(nf,f,nz,zr,np,pl,sg,y,yl)           */
/*************************************************************/
/* Compute the response for all f using the poles and zeros  */
/* Input parameters  : f  : frequency (nf values)            */
/*                     zr : zeros (nz values)                */
/*                     pl : poles (np values)                */
/*                     sg : constant value                   */
/*                                                           */
/* Output parameters : y  : velocity-response (nf doubles)   */
/*                     yl : log10(y) (nf floats)             */
/*                                                           */
void 
compresp (int nf, double *f, double nz, complex *zr, double np, complex *pl, double sg, double **y, float **yl, int *ierror)
{
  int     i ;
  double  w ;
  complex resp ;
  
  (*y)   = double_alloc(nf) ;
  (*yl)    = float_alloc(nf) ;
  for ( i=0; i<nf; i++ )
    {
      w     = 2*M_PI*f[i] ;
      resp  = pz2res (w, nz, zr, np, pl, sg, ierror) ;
      if (*ierror == 99)
	return;
      (*y)[i]  = abs_c(resp) ;
      (*yl)[i] = (float)log10((*y)[i]) ;
    }
}
  
  

/*************************************************************/
/*   firstparams(fref,fl,nz,zr,np,pl,sg,gain_factor,fc,h)    */
/*************************************************************/
/* Compute the first approximation of gain_factor, fc and h  */
/* Input parameters  : fref : reference frequency            */
/*                     fl   : lower frequency                */
/*		       zr   : zeros (nz value)               */
/*                     pl   : poles (np values)              */
/*                     sg   : constant value                 */
/*                                                           */
/* Output parameters : gain_factor : gain factor (1st approx)*/
/*                     fc : natural frequency (1st approx)   */
/*                     h  : damping factor (1st approx)      */
/*                                                           */
void 
firstparams(double fref, double fl, double nz, complex *zr, double np, complex *pl, double sg, double *gain_factor, double *fc, double *h, int *ierror)
{
  double  w0, w, amp01;
  complex resp ;
  /* Compute the (a priori) gain factor at the reference frequency */
  w0           = 2*M_PI*fref ;
  resp         = pz2res (w0, nz, zr, np, pl, sg, ierror) ;
  if (*ierror == 99)
    return;
  *gain_factor = abs_c(resp) ;

  /* Compute the first approximation of fc */
  w     = 2*M_PI*fl ;
  resp  = pz2res (w, nz, zr, np, pl, sg, ierror) ;
  if (*ierror == 99)
    return;
  amp01 = abs_c(resp) ;
  *fc   = fl*sqrt((*gain_factor)/amp01) ;

  /* Compute the first approximation of h */
  w     = 2*M_PI*(*fc) ;
  resp  = pz2res (w, nz, zr, np, pl, sg, ierror) ;
  if (*ierror == 99)
    return;
  amp01 = abs_c(resp) ;
  *h    = (*gain_factor)/(2.*amp01) ;
  
}


/*****************************************************************/
/*           s2=rmsestimates(gain_factor,fc,h,nf,f,y)            */
/*****************************************************************/
/* Error estimates on the evaluation of gain_factor, fc and h    */
/* Input parameters  : gain_factor : gain factor (inverted)      */
/*                     fc : natural frequency (inverted)         */
/*		       h  : damping factor                       */
/*                     f  : frequency (nf values)                */
/*                     y  : "true" response from zeros and poles */
/*                                                               */
/* Output parameter  : s2 : "rms" error                          */
/*                                                               */
double 
rmsestimates(double gain_factor, double fc, double h, int nf, double *f, double *y, int *ierror) 
{
  int     i ;
  double  s2, *yc, w, wc;
  complex num, den;
  
  wc = 2.*M_PI*fc ;
  yc = double_alloc(nf) ;
  s2 = 0. ;
  for (i=0; i<nf; i++)
    {
      w     = 2.*M_PI*f[i] ;
      num   = cmplx(-pow(w,2),0.) ;
      den   = add_c(add_c(cmplx(-pow(w,2),0.),cmplx(0.,2.*h*wc*w)),cmplx(pow(wc,2),0.)) ;
      yc[i] = abs_c(div_c(num,den,ierror)) * gain_factor ; /* inverted velocity-response */
      if (*ierror == 99)
	return 99;
      s2    = s2 + pow(y[i]/yc[i] - 1.,2) ;         
    }
  s2 = sqrt(s2/(double)(nf-1)) ; // pourquoi nf-1 et pas nf???
  free(yc);
  return (s2) ;
}



/*******************************************************************/
/*  plzr2resp(pzfilename,tol,fl,fh,nf,fref,&g,&fc,&h,&s2,&ierror)  */
/*******************************************************************/
/* computation of the gain factor, the natural freq., the damping  */
/* factor, and the rms error associated with this estimation       */
/* Input parameters  : pzfilename : name of the pzfile             */
/*                     tol        : rms tolerance                  */
/*                     fl         : lower frequency  =>|           */
/*                     fh         : higher frequency =>|           */
/*                     nf         : nb of samples    =>|=> used in */
/*                                     the estimates of g,fc and h */
/*                     fref       : reference frequency used for a */
/*                                     first estimates of the gain */
/*                                                                 */
/* Output parameter  : g  : gain factor (inverted)                 */
/*                     fc : natural frequency (inverted)           */
/*		       h  : damping factor                         */
/*                     s2 : "rms" error                            */
/* if s2 > tol ierror is set to 1 if everything is O.K. ierror=0,  */
/* if something is wrong i_error=99                                */
void 
plzr2resp(char *pzfilename , double tolerance, double fl, double fh, int nf , double fref, 
	       double *gain_factor, double *fc, double *h, double *s2, int *ierror)
{
  int     i, nz, np, mp, mv, ip, *ib, idvt, icon, iquit, iprnt;
  double  *f, *y, sg;
  float   *yl,*x,*b;
  complex *zr, *pl;

  *ierror = 0 ;

  /* Set the lsq parameters */
  mp    = 3 ;   /* no. of param */
  mv    = 1 ;   /* no. of indep. variables */
  ip    = 1 ;   /* no. of parameters which are held constant (0)*/
  
  ib    = int_alloc(3) ;
  ib[0] = 1;    /* indices of param which are held constant 
  		   uncomment if you want gain_factor to be constant  
		   (and set ip=1) */
  ib[1] = 0;
  ib[2] = 0;
  idvt  = 1 ;   /* we use "estimated" derivatives */ 
  icon  = 1 ;
  iquit = 0 ;   /* no force off */
  iprnt = 0 ;   /* abbreviated printout */


  /* Set the frequency vector (scaled linearly in log) */
  f = double_alloc(nf) ;
  x = float_alloc(nf) ;
  for ( i=0; i<nf; i++)
    {
      f[i] = pow( fh/fl , (double)i/(double)(nf-1) ) * fl ;
      x[i] = (float)f[i] ;
    }

  /* Read the pzfile */
  readpzfile (pzfilename, &nz, &zr, &np, &pl, &sg, ierror) ;
  if (*ierror == 99)
    return;

  /* Compute the velocity-response */
  compresp (nf, f, nz, zr, np, pl, sg, &y, &yl, ierror) ;
  if (*ierror == 99)
    return;

  /* first approximation of gain_factor, fc and h */
  firstparams(fref, fl, nz, zr, np, pl, sg, gain_factor, fc, h, ierror) ;
  if (*ierror == 99)
    return;
  printf(" 1st approx. GF, fc and h =    %.7e  %.7e  %.7f\n",*gain_factor, *fc, *h) ;

  /* LSQENP */
  b    = float_alloc(3) ;
  b[0] = (float)*gain_factor ; /*  lsq "parameters" */
  b[1] = (float)*fc ;
  b[2] = (float)*h ;
  lsqenp_ (&nf,&mp,&mv,yl,x,b,&ip,ib,&idvt,&icon,&iquit,&iprnt);
  *gain_factor = (double)b[0] ;/* lsq output */
  *fc          = (double)b[1] ;
  *h           = (double)b[2] ;
  
  /* RMS estimates */
  *s2 = rmsestimates(*gain_factor, *fc, *h, nf, f, y, ierror) ;
  if (*ierror == 99)
    return;
  printf(" lsq :  GF, fc, h and rms =    %.7e  %.7e  %.7f      %.7e\n",*gain_factor,*fc,*h,*s2);
  if (*s2 >= tolerance)
    {
      fprintf (stderr, "%s misfit too large (=%f)\n",pzfilename,*s2) ; 
      *ierror = 1 ;
      /*exit (1) ;*/
    }
  /* freeing memory */
  free((void *) b);
  free((void *) zr);
  free((void *) pl);
  free((void *) f);
  free((void *) y);
  free((void *) yl);
  free((void *) x);
}



void 
pzfile2id(char *pzfile, char *id, int *ierror) 
{
  int   i, n, k, ksep[7], shift;
  char *buf;
  
  buf    = char_alloc(6) ;
  id[0]  = '\0' ;
  n      = strlen(pzfile) ;
  k      = 0 ;
  for(i=0; i<n; i++)
    {
      if (pzfile[i] == '_') 
	{
	  ksep[k] = i ;
	  k++ ;
	}
    }

  if (k == 7)
    {
      i = 1 ; 
      pzfile+= ksep[0]+1 ;
      shift  = ksep[1]-ksep[0] ;
      buf[0] = 0 ;
      while ( i < 4 )
	{
	  strcat( id, buf) ;
	  pzfile+= shift ;
	  shift  = ksep[i+1] - ksep[i] ;
	  strncpy( buf, pzfile, shift ) ;
	  buf[shift] = '\0' ;
	  i++ ;
	}
      if ( (ksep[i+1] - ksep[i]) > 1 )
	{
	  pzfile+= shift ;
	  shift  = ksep[i+1] - ksep[i] ;
	  strncat( id, pzfile, shift ) ;
	}
      else
	{
	  strcat(id,"--_");
	}
      strncat( id, buf, strlen(buf)-1) ;
      *ierror = 0 ;
    }
  else
    *ierror = 1;

  /* freeing memory */
  free((void *) buf);
}

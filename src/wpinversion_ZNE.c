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

#ifndef POW
#define POW 1.e28
#endif

#ifndef DEG2RAD
#define DEG2RAD M_PI/180.
#endif

#ifndef RX
#define	RX 20
#endif

#ifndef RY
#define	RY 9
#endif

#ifndef NM
#define NM 6
#endif 

typedef struct
{
  int    Nit ;
  double th_val, cth_val, df_val, med_val  ;
  double op_pa, p2p_med, p2p_low, p2p_high ;
  double dts_min,dts_step,dts_max, dts_val ; 
  double  ntr_val, ref_val, dmin, dmax,azp ;
  double *rms_in, *p2p,*avg, *wgt,wZ,wN,wE ;
  char   i_master[FSIZE], i_saclst[FSIZE]  ;
  char   o_saclst[FSIZE], log[FSIZE]       ; 
  char   o_cmtf[FSIZE], p_data[FSIZE]      ;
  char   osacdir[FSIZE], psfile[FSIZE]     ; 
  char   wpbmfile[FSIZE], refbmfile[FSIZE] ;
  char   gsfile[FSIZE], o_covf[FSIZE]      ; 
} structopt ;


/* *** INTERNAL FUNCTIONS *** */


/* Parameter setting functions */
void get_opt(int numarg1, int numarg2, char **argv, structopt *opt, str_quake_params *eq) ;
void get_param1(int argc, char **argv, int *M, structopt *opt, str_quake_params *eq, int *flag);
void get_param2(char *file, structopt *opt, str_quake_params *eq, int *flag) ;

/* Screening routines */
void screen_med(int *nsac, char **data_name, double **data, 
		double ***G, sachdr *hd_synt, structopt *opt, FILE *o_log) ;
void screen_rms(int *nsac, char **data_name, double **data, 
		double ***G, sachdr *hd_synt, structopt *opt, FILE *o_log) ;

/* Stats routines */
void median(int *nsac, structopt *opt) ;
void sort(double *tab, int *ns) ;
void calc_stat(int npts, double *data, double *p2p, double *ave) ;
void calc_rms_sub(int npts, double *data, double **dcalc, double *rms, 
		   int flag) ;
void calc_rms(int *ns, sachdr *hd_synt, double **data, double ***dcalc, 
	       double **rms, double *global_rms, structopt *opt, int flag) ;
void azpond(sachdr *hd_synt, int ns, double *wgt) ;

/* Inversion routines and linear algebra */
void free_G(double ***G) ;
int  fill_G(char *gf_file, char *datafile, sachdr *hd_GF, sachdr *hd_data, int npts, double Ptt, 
	    double twp_beg, double twp_end, double *buffer, double *G, structopt *opt, FILE *o_log) ;
void comp_GtG( int *M, int *nsac, sachdr *hd_synt, double ***G, double **GtG, structopt *opt) ;
void jacobi(double **a,int n, int np, double *d, double **v, int *nrot) ;
void eigsrt(double *d, double **v, int n) ;
void inversion(int *M, int *nsac, sachdr *hd_synt, double ***G, double **d, 
		 double *vma, double *Cond, structopt *opt, FILE *o_log) ;
void set_matrices (char *i_saclst, double *evdp, double *wp_win4, int *nsac, 
		     int *nsini, char ***sacfiles, sachdr **hd_synt, double ***data, 
		     double ****G, structopt *opt, str_quake_params *eq, FILE *o_log) ;
void calc_data(int *nsac, sachdr *hd_synt, double ***G, double **vm, double **data,
		 double ***d, structopt *opt, int flag) ;
void residual_moment(double **vm, double *ma, double *mb, double *mc) ;
void mt2sm(double *vm,double *sm) ;
double **set_mt(double *vm) ;
void get_planes(double *vm, double ***TM, double **eval3, double *s1, 
		double *d1, double *r1, double *s2,double *d2,double *r2) ;
void realloc_gridsearch(int nsac, double ***rms, double **global_rms, double **vm, 
			double ****dcalc, int flag);
void fast_ts_gridsearch(int nsac, int M, char **sacfiles, sachdr *hd_synt, double **data, 
			double ***G, double ****dcalc, double ***rms, double **global_rms, 
			structopt *opt, str_quake_params *eq, double *rmsini, double *tsopt, 
			double *rmsopt, FILE *o_log);

/* Writing/Printing functions */
void w_log_header(char **argv, structopt *opt, str_quake_params *eq, double *wp_win4, 
		  FILE *o_log) ;
void write_cmtf(char *filename, str_quake_params *eq, double *vm, int flag) ;
void w_o_saclst(char *o_saclst, int *ns, char **sacfiles, sachdr *hd_synt,
		double **rms, structopt *opt, int flag) ;
void get_gap(sachdr *hd_synt, int *ns, double *gap) ;
int charplot(double *M, double s1, double d1, double s2, double d2, 
	     char D, char P, char W, char B, char sep, char pnod, 
	     int rx, int ry, FILE *stream) ;
void output_products(structopt *opt, str_quake_params *eq, double *s1, double *d1, 
		     double *r1, double *s2, double *d2, double *r2, double **TM, 
		     double *eval3, double *M0, double *M0_12, double *Mw, 
		     double *global_rms, double *gap, double *Cond, int *nsac, 
		     sachdr *hd_synt, FILE *o_log, int flag) ;
void prad_pat(double **TM, FILE *ps)               ;
void pnod_pat(double *s1, double *d1, FILE *ps)    ;


int 
main(int argc, char *argv[])
{
  int i, j, nsac, M, flag, nsini   ;

  double s1, d1, r1, s2, d2, r2, gap, Cond            ;
  double tsopt,rmsopt,rmsini,M0, Mw, diplow, M0_12    ;
  double *eval3, *global_rms, **TM  ;
  double **data,  **rms, ***G = NULL, ***dcalc ;
  

  char **sacfiles ;
  FILE *o_log  ;
  structopt opt       ;
  sachdr    *hd_synt  ;
  str_quake_params eq ;  
  
  /* Allocate memory for input parameters */
  eq.wp_win4  = double_alloc(4) ;
  eq.vm    = double_alloc2p(2)  ;
  eq.vm[0] = double_calloc(NM)   ; 
  eq.vm[1] = double_calloc(NM)   ;

  /* Initialize some pointers */
  data       = NULL ; G       = NULL ; opt.wgt = NULL ;
  opt.rms_in = NULL ; opt.p2p = NULL ; opt.avg = NULL ;
  sacfiles   = NULL ; 

  /* Get input parameters */
  get_param1(argc, argv, &M, &opt, &eq, &flag) ;
  get_param2(opt.i_master, &opt, &eq, &flag)    ;
  fflush(stdout) ;
  
  /* Write log header     */
  o_log = openfile_wt(opt.log)                     ;
  w_log_header(argv, &opt, &eq, eq.wp_win4, o_log) ;
  
  /* Set G and data       */
  set_matrices (opt.i_saclst, &eq.pde_evdp, eq.wp_win4, &nsac, &nsini,
		&sacfiles, &hd_synt, &data, &G, &opt, &eq, o_log) ; 
  
  /* Screening            */
  if (opt.med_val > 0.) 
    {
      median(&nsac, &opt)   ;
      screen_med(&nsac, sacfiles, data, G, hd_synt, &opt, o_log) ; 
    }

  if (opt.th_val > 0.)
    screen_rms(&nsac, sacfiles, data, G, hd_synt, &opt, o_log) ;
  
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
  
  /* Inversion */
  inversion(&M, &nsac, hd_synt, G, data, eq.vm[0], &Cond, &opt, o_log) ;
  
  /* Predicted data */
  dcalc = double_alloc3p(nsac) ;
  calc_data(&nsac, hd_synt, G, eq.vm, data, dcalc, &opt, flag);
  
  /* Get RMS and Gap */
  global_rms = double_calloc(2*flag) ;
  rms = double_calloc2(nsac, 2*flag) ;
  calc_rms(&nsac, hd_synt, data, dcalc, rms, global_rms, &opt, flag)  ;
  get_gap(hd_synt, &nsac, &gap) ;
  w_o_saclst(opt.o_saclst, &nsac, sacfiles, hd_synt, rms, &opt, flag) ; 

  /* Time-shift grid search */
  if (opt.dts_step > 0.)
    {
      fast_ts_gridsearch(nsac,M,sacfiles,hd_synt,data,G,&dcalc,&rms,&global_rms,&opt,&eq,&rmsini,&tsopt,&rmsopt,o_log);
      /* Optimum solution */
      opt.dts_val = tsopt ;
      realloc_gridsearch(nsac, &rms, &global_rms, &eq.vm[0], &dcalc, flag) ;
      set_matrices (opt.i_saclst, &eq.pde_evdp, eq.wp_win4, &nsac, &nsini,
		    &sacfiles, &hd_synt, &data, &G, &opt, &eq, o_log) ;    
      inversion(&M, &nsac, hd_synt, G, data, eq.vm[0], &Cond, &opt, o_log) ;   
      calc_data(&nsac, hd_synt, G, eq.vm, data, dcalc, &opt, flag)    ;
      calc_rms(&nsac, hd_synt, data, dcalc, rms, global_rms, &opt, flag)    ;
      eq.ts += opt.dts_val ;
    }
  
  /* Set stike/dip/rake */
  get_planes(eq.vm[0], &TM, &eval3, &s1,&d1,&r1, &s2,&d2,&r2) ;
  write_cmtf(opt.o_cmtf, &eq, eq.vm[0], flag) ;  
  
  /* Set Moment and Magnitude (Harvard Definition) */
  M0     = ((fabs(eval3[0]) + fabs(eval3[2])) * (double)POW) / 2. ; 
  Mw     = (log10(M0) - 16.1) / 1.5 ;
  diplow = d2 ;
  if (d1 < d2) diplow = d1 ;
  M0_12  = M0 * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ;

  /* Output */
  output_products(&opt, &eq, &s1,&d1,&r1, &s2,&d2,&r2, TM, eval3, &M0, &M0_12, &Mw, 
		  global_rms, &gap, &Cond, &nsac, hd_synt, o_log, flag) ;
  fclose(o_log);

  /* Memory Freeing */
  free((void*)eq.wp_win4) ;
  free((void*)eq.vm[0])        ;
  free((void*)eq.vm[1])   ;
  free((void**)eq.vm)     ;
  free((void*)eval3)      ;
  free((void*)global_rms) ;
  for(i=0 ; i<nsac ; i++)
    {
      free((void*)data[i])     ;
      free((void*)rms[i] )     ;
      free((void*)sacfiles[i]) ;
      free_G(&G[i]);
      for(j=0 ; j<flag ; j++)
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
  return 0;
  
  free((void*)opt.rms_in) ;
  free((void*)opt.p2p)    ;
  free((void*)opt.avg)    ;
  free((void*)opt.wgt)    ;	 
  free((void*)hd_synt)    ;
}


void 
output_products(opt, eq, s1a, d1a, r1a, s2a, d2a, r2a, TMa, eval3a, M0a, M0_12a, Mwa, 
		     global_rms, gap, Cond, nsac, hd_synt, o_log, flag) 
     int flag, *nsac ;
     double *s1a, *s2a, *d1a, *d2a, *r1a, *r2a ;
     double *M0a, *M0_12a, *Mwa, **TMa         ;
     double  *eval3a, *global_rms, *gap, *Cond ;
     structopt *opt       ;
     sachdr    *hd_synt   ;
     str_quake_params *eq ;  
     FILE       *o_log    ;
{
  int    i,j,k, nb,nb2, size, *nbcmp  ;
  double M0b, Mwb, M0_12b, ma, mb, mc ;
  double s1b, s2b, d1b, d2b, r1b, r2b ;
  double diplow, **TMb, *eval3b       ;
  char   *buf, *buf2, **sta, **cmp    ;
  FILE   *ps, *o_bitmap               ;

  char   date_stmp[64] ;
  time_t now           ; 
  
  nbcmp = int_alloc(*nsac)        ;
  buf   = char_alloc(9)           ;
  buf2  = char_alloc(32)          ;
  sta   = char_calloc2(*nsac, 9)  ;
  cmp   = char_calloc2(*nsac, 30) ;  
  
  /* PS FILE */
  ps = openfile_wt(opt->psfile) ;
  /* header */
  fprintf(ps,"%%!PS\n300 600 translate\n100 100 scale\n");
  fprintf(ps,"0 setlinewidth\n") ;

  /* Focal mechanism display - WPhase Solution */
  prad_pat(TMa, ps)              ;
  pnod_pat(s1a, d1a, ps)           ;
  pnod_pat(s2a, d2a, ps)           ;
  fprintf(ps,"-0.3 000 translate\n") ;

  /* title */
  fprintf(ps,"/Times-Roman findfont .2 scalefont setfont\n") ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., +1.3) ;
  fprintf(ps,"(%s) show\n", eq->evnm) ;
 
  /* Time stamp */
  time(&now) ;
  fprintf(ps,"/Times-Roman findfont .1 scalefont setfont\n") ;
  strftime(date_stmp,64,"Processed Date : %a %b %02d %02H:%02M:%02S %Y UT",gmtime(&now));
  date_stmp[44] = '\0';
  fprintf(ps,"%15.6f %15.6f moveto\n", -2.45, 1.7) ;
  fprintf(ps,"(%s) show\n", date_stmp) ;

  /* tensor elem., moment, planes, eigvalues */
  fprintf(ps,"/Times-Roman findfont .1 scalefont setfont\n") ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.45) ;
  fprintf(ps,"(Moment Tensor:") ;
  for(i=0 ; i<NM ; i++)
    fprintf(ps," %11.5f",eq->vm[0][i]) ;
  fprintf(ps,") show\n") ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.6) ;
  fprintf(ps,"(Scalar moment: %9.2e dyn cm    (Mo_12 = %9.2e dyn cm)) show\n"
	  , (*M0a), (*M0_12a)) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.75) ;

  fprintf(ps,"(Best nodal planes: %7.1f/%5.1f/%7.1f   %7.1f/%5.1f/%7.1f) show\n",
	  *s1a,*d1a,*r1a,*s2a,*d2a,*r2a) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.9) ;
  fprintf(ps,"(Eigenvalues:%13.5f %13.5f %13.5f (Mw = %4.2f)) show\n",
	  eval3a[0],eval3a[1],eval3a[2], *Mwa) ;

  /* fit quality */
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.1) ;
  fprintf(ps,"(Fit Quality:) show\n") ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -0.85, -2.25) ;
  fprintf(ps,"(WCMT - RMS: %9.5f mm (%6.3f),  Gap: %5.1f\\312,  C#%9.0f) show\n",
	  1000.*global_rms[0], global_rms[0]/global_rms[1], *gap, *Cond);
  if (flag == 2) 
    {
      fprintf(ps,"%15.6f %15.6f moveto\n", -0.85, -2.4) ;
      fprintf(ps,"(GCMT - RMS = %9.5f mm (%6.3f)) show\n", 1000.*global_rms[2], global_rms[2]/global_rms[3]);  
    }

  /* used stations */
  k = 0 ;
  for(i=0; i<*nsac; i++) /* Set list of channels per stations */
    {
      nb = nbchar(hd_synt[i].kstnm)   ;
      nb2 = nbchar(hd_synt[i].kcmpnm) ;      
      strncpy(buf,hd_synt[i].kstnm,nb) ;
      strncpy(buf2,hd_synt[i].kcmpnm,nb2) ;
      buf[nb]   = '\0' ;
      buf2[nb2] = '\0' ;
      for(j=0; j<k; j++)
	{
	  if (strcmp(buf,sta[j]) == 0) 
	    {
	      strcat(cmp[j],",")   ;
	      strcat(cmp[j], buf2) ;
	      nbcmp[j]++;
	      break ; 
	    }
	}
      if (j==k) 
	{
	  strcpy(sta[k], buf)  ;
	  strcpy(cmp[k], buf2) ;
	  nbcmp[k] = 1         ;
	  k++ ; 
	}
    }
  
  fprintf(ps,"/Courier-Bold findfont .1 scalefont setfont\n") ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.6)             ;
  fprintf(ps,"(Used stations (%d, %d channels) : ) show\n", k, *nsac ) ;
  fprintf(ps,"/Courier      findfont .07 scalefont setfont")  ;

  j  = 0;
  nb = 1;
  fprintf(ps,"%15.6f %15.6f moveto\n(", -1., -2.7) ;
  for(i=0; i<k; i++)
    {
      j += 1 + nbcmp[i] ;
      fprintf(ps, "  %s(%s)", sta[i], cmp[i]) ;
      if (j >= 16) {
	fprintf(ps,") show\n") ;
	fprintf(ps,"%15.6f %15.6f moveto\n(", -1., -2.7-((double)nb)/8.) ;
	j = 0 ;
	nb++; }
    }
  fprintf(ps,") show\n")   ;
  
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.7-((double)nb)/8.) ;
  fprintf(ps, "(WPWIN: %-8.2f %-8.2f %-8.2f %-8.2f ) show\n"
	  , eq->wp_win4[0], eq->wp_win4[1], eq->wp_win4[2], eq->wp_win4[3]) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.8-((double)nb)/8.) ;
  fprintf(ps, "(Dmin : %-8.2f Dmax :%-8.2f) show\n", opt->dmin, opt->dmax) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.9-((double)nb)/8.) ;
  fprintf(ps, "(wN   : %-8.2f wE   :%-8.2f wZ   :%-8.2f) show\n", opt->wN, opt->wE, opt->wZ);  

  /* filter parameters */
  fprintf(ps,"/Courier-Bold findfont .1 scalefont setfont") ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -3.8) ;
  fprintf(ps,"(Filter parameters: ) show\n") ;
  fprintf(ps,"/Courier     findfont .1 scalefont setfont")  ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -3.9)   ;
  fprintf(ps,"(filt_order: %-d) show\n", eq->filtorder) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.025) ;
  fprintf(ps,"(filt_cf1  : %-7.5f) show\n", eq->flow) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.150) ;
  fprintf(ps,"(filt_cf2  : %-7.5f) show\n", eq->fhigh) ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.275) ;
  fprintf(ps,"(filt_pass : %-d) show\n", eq->filtnpass) ;
  
  /* pde and reference solution */
  fprintf(ps,"/Courier-Bold findfont .07 scalefont setfont\n");
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.5) ;
  fprintf(ps,"(PDE and Centroid: ) show\n")       ;
  fprintf(ps,"/Courier      findfont .07 scalefont setfont\n");
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.6)    ;
  nb = nb_blank(eq->pdeline) ;
  eq->pdeline[60]='\0' ;
  fprintf(ps,"(%s) show\n", &eq->pdeline[nb])             ;
  fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.7)         ;
  fprintf(ps, "(Event id     : %s) show\n", eq->evid)     ;
  fprintf(ps, "%15.6f %15.6f moveto\n", -1., -4.8)        ;
  fprintf(ps, "(Time shift   : %-6.1f s)  show\n", eq->ts);
  fprintf(ps, "%15.6f %15.6f moveto\n", -1., -4.9)        ;
  fprintf(ps, "(Half duration: %-6.1f s)  show\n", eq->hd);
  fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.0)        ;
  
  fprintf(ps, "(Latitude     : %-8.3f) show\n", eq->evla) ;
  fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.1)        ;
  fprintf(ps, "(Longitude    : %-8.3f) show\n", eq->evlo) ;
  fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.2)        ;
  fprintf(ps, "(Depth        : %-8.3f) show\n", eq->evdp) ;
  fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.3)        ;  
  if (flag == 2) 
    {
      /* eq->pdeline[60]=' ' ; */
      /* fprintf(ps, "(%s) show\n", &eq->pdeline[61]) ;*/
      get_planes(eq->vm[1], &TMb, &eval3b, &s1b,&d1b,&r1b, &s2b,&d2b,&r2b) ;
      M0b = ((fabs(eval3b[0]) + fabs(eval3b[2])) * (double)POW) / 2. ; 
      Mwb = (log10(M0b) - 16.1) / 1.5 ;
      residual_moment(eq->vm, &ma, &mb, &mc) ;
      fprintf(ps,"1.8 -0.6 translate\n") ;
      fprintf(ps,"0.4  0.4 scale\n") ;
      fprintf(ps,"/Times-Roman findfont .2 scalefont setfont\n") ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -1.1, -1.4) ;
      fprintf(ps,"(GCMT, Mw= %5.2f ) show\n", Mwb)      ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -1.1, -1.6) ;
      /* fprintf(ps,"(ratio = %5.2f ;  epsilon = %6.3f) show\n",M0b/(*M0a),mc) ;*/
      fprintf(ps,"(ratio = %5.2f ;  epsilon = %6.3f) show\n",M0b/(*M0a),mc) ;
      fprintf(ps,".5 .5 .9 setrgbcolor\n") ;
      prad_pat(TMb, ps)        ;
      pnod_pat(&s1b, &d1b, ps) ;
      pnod_pat(&s2b, &d2b, ps) ;
    }

  
  fprintf(ps,"showpage\n") ;
  fclose(ps) ;


  /* STDOUT AND LOG */
  charplot(eq->vm[0], *s1a,*d1a, *s2a,*d2a, '-', '#', ' ', '\0','\0','\0', RX, RY, stdout)  ;
  printf("CMT: ") ;
  for(i =0 ; i<NM ; i++)
    printf("%15.4e",eq->vm[0][i]);
  printf("\nBest nodal planes: %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n",
	 *s1a, *d1a, *r1a, *s2a, *d2a, *r2a);
  printf("Eigenvalues: %-12.5f %-12.5f %-12.5f\n", eval3a[0], eval3a[1], eval3a[2]) ;
  printf("WCMT: RMS = %-9.5f mm (%-6.3f) Gap: %-6.1f deg, C#: %-10.0f\n",
	 1000.*global_rms[0], global_rms[0]/global_rms[1], *gap, *Cond);
  if (flag == 2) 
    { 
      printf("GCMT: RMS = %-9.5f mm (%-6.3f)\n", 1000.*global_rms[2], global_rms[2]/global_rms[3])  ;      
      diplow = d2b ;
      if (d1b < d2b)
	diplow = d1b ;
      M0_12b = M0b * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ; 
    }
  printf("Wmag: %-5.2f ; Wmom %-15.4e ; Wmom_12 %-15.4e\n",*Mwa,*M0a,*M0_12a) ;
  

  fprintf(o_log,"n_used_rec:         %-4d\n", *nsac)    ;
  fprintf(o_log,"Gap:                %-7.1f\n", *gap)   ;
  fprintf(o_log,"Cond_number:        %-10.0f\n", *Cond) ; 
  fprintf(o_log,"W_bestnodal planes: %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n",
	 *s1a, *d1a, *r1a, *s2a, *d2a, *r2a) ;
  fprintf(o_log,"W_eigenvalues:      %-12.5f %-12.5f %-12.5f\n", eval3a[0], eval3a[1], eval3a[2]) ;
  fprintf(o_log,"W_cmt_err:          %-12.8f %-12.8f\n", 1000.*global_rms[0], global_rms[0]/global_rms[1])    ;
  fprintf(o_log,"Wmag: %-5.2f ; Wmom %-15.4e ; Wmom_12 %-15.4e\n",*Mwa,*M0a,*M0_12a) ;

  if (flag == 2) 
    {
      printf("Rmag: %-5.2f ; Rmom %-15.4e ; Rmom_12 %-15.4e\n",Mwb,M0b,M0_12b) ;
      printf("ratio = %5.2f ;  epsilon = %6.3f\n",M0b/(*M0a),mc) ;
      fprintf(o_log,"R_bestnodal planes: %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n",s1b,d1b,r1b, s2b,d2b,r2b) ;
      fprintf(o_log,"R_eigenvalues:      %-12.5f %-12.5f %-12.5f\n", eval3b[0], eval3b[1], eval3b[2])           ;
      fprintf(o_log,"R_cmt_err:          %-12.8f %-12.8f\n", 1000.*global_rms[2], global_rms[2]/global_rms[3])  ;
      fprintf(o_log,"Rmag: %-5.2f ; Rmom %-15.4e ; Rmom_12 %-15.4e\n",Mwb,M0b,M0_12b)                           ;
      fprintf(o_log,"ratio = %12.8f ;  epsilon = %12.8f\n",M0b/(*M0a),mc) ;
    }

  /* BITMAP IMAGE */
  if (strlen(opt->wpbmfile) != 0)
	{
	  size = 401 ;
	  o_bitmap = openfile_wt(opt->wpbmfile)  ;
	  fprintf(o_bitmap, "P2\n#WPinversion focal mechanism\n%d %d\n9\n", size, size)                             ;
	  charplot(eq->vm[0], *s1a, *d1a, *s2a, *d2a,'9', '3', '9', '0',' ', '0', (size-1)/2,(size-1)/2, o_bitmap)  ;
	  fclose(o_bitmap) ;
	}

  if (strlen(opt->refbmfile) != 0 && flag == 2) 
	{
	  o_bitmap = openfile_wt(opt->refbmfile) ;
	  fprintf(o_bitmap, "P2\n#WPinversion focal mechanism\n%d %d\n9\n", size, size)                         ;
	  charplot(eq->vm[1], s1b, d1b, s2b, d2b,'9', '3', '9', '0',' ', '0', (size-1)/2,(size-1)/2, o_bitmap)  ;
	  fclose(o_bitmap)                       ; 
	}

  /* Memory Freeing */
  if (flag == 2)
    {
      for(i=0 ; i<3 ; i++) 
	free((void*)TMb[i]) ;
      free((void**)TMb)     ;
      free((void*)eval3b)   ;
    }
  free((void*)nbcmp)    ;
  free((void*)buf)      ;
  free((void*)buf2)     ;
  for(i=0 ; i<*nsac ; i++){
    free((void*)sta[i]) ;
    free((void*)cmp[i]) ; }
  free((void**)cmp)     ;
  free((void**)sta)     ;
}

/*******************************/
/*       prad_pat(TM, ps)      */
/*******************************/
/* Plot ps beach ball          */
/* Input: TM : moment tensor   */
/*        ps : ps FILE stream  */
void 
prad_pat(double **TM, FILE *ps)
{
  int    i, j, k, l, N = 50, M ;
  double x, y, ymax, rad ;
  double azi, ain, rpp, *r  ;

  r = double_alloc(3) ;

  fprintf(ps,"newpath\n") ;
  for(i=-N; i<=N; i++)
    {
      x    = ((double)i)/((double)N) ;
      ymax = sqrt(1.-x*x) ;
      M = (int) (((double)N)*ymax);
      if (M > 0) 
	{
	  for(j=-M; j<=M; j++)
	    {
	      y    = ((double)j)/((double)N) ;
	      rad  = sqrt(x*x+y*y)           ;
	      azi  = atan2(x,y)              ;
	      ain  = 2.*asin(rad/sqrt(2.))   ;
	      r[0] = sin(ain)*cos(azi)       ;
	      r[1] = sin(ain)*sin(azi)       ;
	      r[2] = cos(ain)                ;
	      rpp   = 0.                     ;

	      for(k=0; k<3; k++)
		for(l=0; l<3; l++)
		  rpp = rpp + r[k]*TM[k][l]*r[l] ;
	      if(rpp > 0.) {
		fprintf(ps,"%7.2f %7.2f .01 add exch moveto\n",y,x) ;
		fprintf(ps,"%7.2f %7.2f 0.01 0 360 arc\n",x,y)      ; }
	    }
	}
    }
  fprintf(ps,"1 0 moveto\n")      ;
  fprintf(ps,"0 0 1 0 360 arc\n") ;
  fprintf(ps,"stroke\n")          ;
  free((void*)r) ;
}

/***************************************/
/*           pnod_pat(s, d, ps)        */
/***************************************/
/* Plot nodal planes on the beach ball */
/* Input: s: strike (pointer)          */
/*        d: dip (pointer)             */
/*        ps: ps FILE stream           */
void 
pnod_pat(s, d, ps)
     double *s, *d ;
     FILE   *ps    ;
{
  fprintf(ps,"newpath\n" )                ;
  fprintf(ps, "/phi %6.1f def \n", (*s))     ;
  fprintf(ps,"/delta %6.1f def\n", (*d))     ;
  fprintf(ps,"phi sin phi cos moveto\n")  ;
  fprintf(ps,"0 1 180\n")                 ;
  fprintf(ps,"{/alpha exch def\n")        ;
  fprintf(ps,"/x alpha cos phi cos mul alpha sin phi sin delta cos mul mul sub def\n") ;
  fprintf(ps,"/y alpha cos phi sin mul alpha sin phi cos delta cos mul mul add def\n") ;
  fprintf(ps,"/z                       alpha sin         delta sin mul         def\n") ;
  fprintf(ps,"x abs 0.001 gt x abs 0.001 gt or {\n") ;
  fprintf(ps,"/PH y x atan def\n")                   ;
  fprintf(ps,"/I0   z dup dup mul neg 1 add sqrt exch atan def\n") ;
  fprintf(ps,"/R  I0 2 div sin 2 sqrt mul def\n")    ;
  fprintf(ps,"PH sin R mul PH cos R mul lineto\n")   ;
  fprintf(ps,"} if } for\n")                         ;
  fprintf(ps,"stroke\n")                             ;

}



/***************************************/
/*    get_gap(hd_synt, ns, gap)        */
/***************************************/
/* Compute gap                         */
/* Input: hd_synt: structs sac headers */
/*        ns     : nb of stations      */
/* Output: gap                         */
void 
get_gap(hd_synt, ns, gap)
     int    *ns      ; 
     double *gap     ;
     sachdr *hd_synt ;
{
  double *tmp, az, dist ;
  int i ;
  
  tmp = double_alloc(*ns) ;
  for( i=0 ; i<*ns ; i++)
    {
      az = hd_synt[i].az ;
      if(az < 0.) 
	tmp[i] = az + 360. ;
      else
	tmp[i] = az ;
    }
  
  sort(tmp, ns) ;
  
  *gap = tmp[0] + 360. - tmp[*ns-1] ;
  for( i=0 ; i<(*ns)-1 ; i++)
    {
      dist = tmp[i+1]-tmp[i] ;
      if(dist > *gap) *gap = dist ;
    }
  free((void*)tmp) ;
}


void 
w_o_saclst(o_saclst, ns, sacfiles, hd_synt, rms, opt, flag) 
     int    *ns, flag ;
     double **rms ;
     char *o_saclst, **sacfiles ;
     sachdr *hd_synt ; 
     structopt *opt  ;
{
  int  i, n0;
  FILE *o_sac ;

  o_sac = openfile_wt(o_saclst);
  n0 = 0;
  for(i=0; i<*ns; i++)
    {
      fprintf(o_sac,"%-65s %8.2f %8.2f %6d %6d %6d %6d %15.8f %15.8f %15.8f %15.8f %10.3f\n",
	      sacfiles[i], hd_synt[i].az, hd_synt[i].gcarc, n0, n0 + hd_synt[i].npts, 
	      (int)hd_synt[i].user[0], (int)hd_synt[i].user[1], rms[i][0], 
	      rms[i][0]/rms[i][1], opt->p2p[i], opt->avg[i], opt->wgt[i]);
      n0 += hd_synt[i].npts ;
    }
  fclose(o_sac);
  
  if (flag == 2)
    {
      o_sac = openfile_wt("_ref_o_wpinversion");
      n0 = 0;
      for(i=0; i<*ns; i++)
	{
	  fprintf(o_sac,"%-65s %8.2f %8.2f %6d %6d %6d %6d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
		  sacfiles[i], hd_synt[i].az, hd_synt[i].gcarc, n0, n0 + hd_synt[i].npts, 
		  (int)hd_synt[i].user[0], (int)hd_synt[i].user[1], rms[i][0], 
		  rms[i][0]/rms[i][1], opt->p2p[i], opt->avg[i],rms[i][2],rms[i][2]/rms[i][3]);
	  n0 += hd_synt[i].npts ;
	}
      fclose(o_sac);
    }

}


/***********************************************/
/*       calc_stat(npts, data, mi2ma, ave)     */
/***********************************************/
/* Compute average amplitude and max-min value */
/* Input: npts : nb of samples                 */
/*        data : data points (pointer)         */
/* Output: mi2ma : max-min value (pointer)     */
/*        ave    : average amplitude (pointer) */
void 
calc_stat(npts, data, mi2ma, ave) 
     int    npts ; 
     double *data, *mi2ma, *ave ;
{
  int    i    ;
  double mini ; 
  double maxi ;

  *ave  = 0. ;
  mini = data[0] ;
  maxi = data[0] ;
  for(i=0 ; i<npts ; i++) {
    *ave = *ave + data[i] ;
    if ( data[i] < mini ) mini = data[i] ;
    if ( data[i] > maxi ) maxi = data[i] ;}
  
  *mi2ma  = maxi - mini ;
  (*ave) /= npts        ;
}


/*************************************************************/
/* calc_rms(ns, hd_synt, data, dcalc, rms, global_rms, flag) */
/*************************************************************/
/* calculate rms amplitudes, rms deviations                  */
/* Input: ns : nb of channels (pointer)                      */
/*        hd_synt : array of sac headers                     */
/*        data    : array of data                            */
/*        flag    : if flag=1 no reference solution          */
/*        flag    : if flag=2 WPhase + reference solution    */
/* Output: dcalc  : array of predicted data                  */
/*         rms    : array of rms amplitudes and deviations   */
/*         global_rms : overall rms amp. and rms deviations  */
void 
calc_rms(ns, hd_synt, data, dcalc, rms, global_rms, opt, flag) 
     int    *ns, flag          ;
     double **data, ***dcalc   ;
     double *global_rms, **rms ;
     sachdr *hd_synt           ;
     structopt *opt            ;
{
  int    i, j, n0, f2;

  f2  = 2*flag;
  n0 = 0 ;
  for(i=0 ; i<*ns ; i++)
    {
      calc_rms_sub(hd_synt[i].npts, data[i], dcalc[i], rms[i], flag) ;
      n0 += hd_synt[i].npts ;
      for(j=0 ; j<f2 ; j++) 
	{
	  global_rms[j] += rms[i][j] ;
	  rms[i][j] = sqrt(rms[i][j]/(double)hd_synt[i].npts); 
	}
    }
  for(j=0; j<f2; j++)
    global_rms[j] = sqrt(global_rms[j]/n0) ;
}


void 
calc_rms_sub(npts, data, dcalc, rms, flag) 
     int    npts, flag   ;
     double *rms, *data, **dcalc ;
{
  int i;

  if (flag == 1)
    {
      rms[0] = 0. ;
      rms[1] = 0. ;
      for(i=0 ; i<npts ; i++){
	rms[0] += (data[i]-(*dcalc)[i]) * (data[i]-(*dcalc)[i]) ;
	rms[1] +=  (*dcalc)[i]*(*dcalc)[i];}
    }
  else if (flag == 2)
    {
      rms[0] = 0. ;
      rms[1] = 0. ;
      rms[2] = 0. ;
      rms[3] = 0. ;
      for(i=0 ; i<npts ; i++){
	rms[0] += (data[i]-dcalc[0][i])*(data[i]-dcalc[0][i]) ;
	rms[1] += dcalc[0][i]*dcalc[0][i] ;
	rms[2] += (data[i]-dcalc[1][i])*(data[i]-dcalc[1][i]) ;
	rms[3] += dcalc[1][i]*dcalc[1][i] ;}      
    }
  else{
    fprintf(stderr, "ERROR comp. rms (wrong flag value)\n");
    exit(1);}
}


/*************************************************************/
/*     calc_data(ns, hd_synt, G, vm, data, d, opt, flag)     */
/*************************************************************/
/* Compute predicted data                                    */
/* Input: ns : nb of channels (pointer)                      */
/*        hd_synt : array of sac headers                     */
/*        hd_synt : array of sac headers                     */
/*        G       : array of Green's functions               */
/*        vm      : pointer on moment tensor elements        */
/*        flag    : if flag=1 no reference solution          */
/*        flag    : if flag=2 WPhase + reference solution    */
/* Output: data   : array data samples                       */
/*         d      : array of predicted data samples          */
/*         opt    : structopt containing output filenames    */
void 
calc_data(nsac, hd_synt, G, vm, data, d, opt, flag)
     int    *nsac, flag      ;
     double **vm, ***G , **data, ***d ;
     sachdr *hd_synt         ;
     structopt *opt  ;
{
  int  i, j, k, s, N        ;
  char *fZ, *fL, *fT, *file ;
  FILE *o_dat, *o_Z, *o_N, *o_E, *ocmp ;

  /* Allocating memory */
  fZ = char_alloc(FSIZE) ;
  fL = char_alloc(FSIZE) ;
  fT = char_alloc(FSIZE) ;

  /* Set filenames     */
  strcpy(fZ,opt->p_data) ;
  strcpy(fL,opt->p_data) ;
  strcpy(fT,opt->p_data) ;
  strcat(fZ,"_LHZ") ;
  strcat(fL,"_LHN") ;
  strcat(fT,"_LHE") ;
  
  /* Opening files */
  o_dat = openfile_wt(opt->p_data) ;
  o_Z = openfile_wt(fZ) ;
  o_N = openfile_wt(fL) ;
  o_E = openfile_wt(fT) ;

  for( s=0 ; s< *nsac ; s++)
    {
      if (hd_synt[s].kcmpnm[2] == 'Z')
	ocmp = o_Z ;
      else if (hd_synt[s].kcmpnm[2] == 'N')
	ocmp = o_N ;
      else if (hd_synt[s].kcmpnm[2] == 'E')
	ocmp = o_E ;
      else 
	{
	  fprintf(stderr,"ERROR : invalid component\n") ;
	  fclose(o_dat) ;
	  fclose(o_Z)   ;
	  fclose(o_N)   ;
	  fclose(o_E)   ;
	  exit(1)       ; 
	}
      N = hd_synt[s].npts ;
      d[s] = double_calloc2(flag, N) ;
      for(i=0 ; i<N ; i++)
	{
	  fprintf(o_dat,"%15.6e ",data[s][i]) ;
	  fprintf( ocmp,"%15.6e ",data[s][i]); 
	  for(j=0; j<flag; j++)  
	    {
	      for(k=0; k<NM; k++)
		d[s][j][i] += G[s][k][i] * vm[j][k] ;
	      fprintf(o_dat,"%15.6e ",d[s][j][i]); 
	      fprintf( ocmp,"%15.6e ",d[s][j][i]); 
	    }
	  fprintf(o_dat,"\n") ;
	  fprintf(ocmp,"\n") ;
	}
      file = get_gf_filename(opt->osacdir, hd_synt[s].kstnm, hd_synt[s].knetwk, 
			     hd_synt[s].kcmpnm, "synth.sac") ;
      whdrsac(file, &hd_synt[s])          ;
      wdatsac(file, &hd_synt[s], d[s][0]) ;
      /* Memory Freeing */
      free((void*)file) ;
    }
  fclose(o_dat) ;
  fclose(o_Z)   ;
  fclose(o_N)   ;
  fclose(o_E)   ;
  free((void*)fZ)   ;
  free((void*)fL)   ;
  free((void*)fT)   ;
}

void 
write_cmtf(filename, eq, vm,flag)
     double *vm     ;     
     char *filename ;
     str_quake_params *eq ;
     int flag ;
{
  int   i ;
  char  Mcmp[6][4] = {"Mrr","Mtt","Mpp","Mrt","Mrp","Mtp"} ;
  FILE *cmtfile;
  
  cmtfile = openfile_wt(filename)     ;

  fprintf(cmtfile,"%s",eq->pdeline) ;
  fprintf(cmtfile,"event name:   %15s\n", eq->evid )  ;
  fprintf(cmtfile,"time shift:   %9.4f\n", eq->ts )   ;
  fprintf(cmtfile,"half duration:%9.4f\n", eq->hd )   ;
  fprintf(cmtfile,"latitude:     %9.4f\n", eq->evla ) ;
  fprintf(cmtfile,"longitude:    %9.4f\n", eq->evlo ) ;
  fprintf(cmtfile,"depth:        %9.4f\n", eq->evdp ) ;
  for (i=0 ; i<NM ; i++)
    fprintf(cmtfile,"%3s:     %14.6e\n", Mcmp[i],vm[i]*(double)POW ) ;
  fclose(cmtfile) ;
}



void 
comp_GtG( M, nsac, hd_synt, G, GtG, opt)
     int    *M, *nsac   ; 
     double ***G, **GtG ;
     sachdr *hd_synt    ; 
     structopt *opt     ;
{
  int i, j, k, s, N;

  if (GtG[0][0] != 0.)
    {
      fprintf(stderr,"WARNING GtG: memory must be set to zero (using calloc)\n");
      for( i=0 ; i<*M ; i++)
	for( j=0 ; j<*M ; j++)
	  GtG[i][j] = 0. ;
    }
  if (*M == NM-1)      /* Constrain null trace */
    {
      for( s=0 ; s<*nsac ; s++)
	{
	  N = hd_synt[s].npts ;
	  for( i=0 ; i<2 ; i++ )
	    {
	      for( j=0 ; j<2 ; j++ )
		for( k=0 ; k<N ; k++ )
		  GtG[i][j] += (G[s][i][k]-G[s][2][k])*opt->wgt[s]*(G[s][j][k]-G[s][2][k]) ;
	      for( j=2 ; j<*M ; j++ )
		for( k=0 ; k<N ; k++ )
		  GtG[i][j] += (G[s][i][k]-G[s][2][k])*opt->wgt[s]*G[s][j+1][k] ;
	    }
	  for( i=2 ; i<*M ; i++ )
	    {
	      for( j=0 ; j<2 ; j++ )
		for( k=0 ; k<N ; k++ )
		  GtG[i][j] += G[s][i+1][k]*opt->wgt[s]*(G[s][j][k]-G[s][2][k]) ;
	      for( j=2 ; j<*M ; j++ )
		for( k=0 ; k<N ; k++ )
		  GtG[i][j] += G[s][i+1][k]*opt->wgt[s]*G[s][j+1][k] ;
	    }
	}
    }
  else if (*M == NM) /* No constraints      */
    for( s=0 ; s<*nsac ; s++)
      for( i=0 ; i<*M ; i++ )
	for( j=0 ; j<*M ; j++ )
	  for( k=0 ; k<hd_synt[s].npts ; k++ )
	    GtG[i][j] += G[s][i][k]*opt->wgt[s]*G[s][j][k] ; 
}



void 
inversion(M, nsac, hd_synt, G, d, vma, Cond, opt, o_log)
     int *M, *nsac                 ; 
     double ***G, **d, *vma, *Cond ;
     sachdr *hd_synt ; 
     structopt *opt  ;
     FILE *o_log     ;
{
  int    i, j ,k, l, s, nk, nrot, N  ;
  double **GtG, *eigvals, **eigvects, **cov ;
  FILE   *o_cov ;

  /* Allocating memory */
  GtG      = double_calloc2 ( *M, *M )   ;
  eigvals  = double_alloc( *M )          ;
  eigvects = double_alloc2( *M, *M )     ;
  cov      = double_calloc2( *M, *M )    ;

  /* Azimuth ponderation */
  if (opt->azp > 0.)
    azpond(hd_synt,*nsac,opt->wgt) ;

  /* Constrains null trace (if M=5) and computes GtG */
  comp_GtG(M,nsac, hd_synt, G, GtG,opt) ;

  /* Computes eigenvalues and eigenvectors */
  jacobi(GtG, *M, *M, eigvals, eigvects, &nrot);
  eigsrt( eigvals, eigvects, *M) ;
  
  /* Conditioning */
  nk = *M;
  if (opt->cth_val <= 0.) /* Remove the eigenvalues smaller than max(eigval)/10000 if any */
    {
      *Cond = eigvals[0] / eigvals[nk-1] ;
      if (*Cond > 1.e4)
	for (i=1 ; i<nk ; i++) 
	  if(eigvals[i] < eigvals[0]/1.e4)  {
	    nk = i ;
	    break  ; }
    }
  else                   /* Damping factor applied to small eigenvalues */
    {
      if (eigvals[0] > opt->cth_val*eigvals[nk-1])
	{
	  fprintf(stderr,"##### Warning:   damping\n") ;
	  fprintf(stderr,"Conditioning number:    %e\n", eigvals[0]/eigvals[nk-1]) ;
	  fprintf(stderr,"Conditioning threshold: %e\n", opt->cth_val) ;
	  fprintf(stderr,"Damping factor:         %e\n", opt->df_val)  ;
	  for(i=1; i<=nk ; i++)
	    eigvals[i] = eigvals[i] + eigvals[0]*(opt->df_val) ;
	}
      *Cond = eigvals[0] / eigvals[nk-1] ;
    }
  if (opt->dts_step <= 0.)
    {
      printf("##############\n%d significant eigenvalues: \n", nk) ;
      for (i=0 ; i<nk; i++) 
	printf("  %13.5e",eigvals[i]);
      printf("\n");
    }
  fprintf(o_log,"Inversion: %d significant eigenvalues: \n", nk) ;  
  for (i=0 ; i<nk; i++) 
    fprintf(o_log,"%e\t",eigvals[i]) ;   
  fprintf(o_log,"\n") ;
  if(*M-nk>0)
    {
      printf(       "%d removed: \n", *M-nk) ;
      fprintf(o_log,"%d removed: \n", *M-nk) ;
      for (i=nk ; i<*M; i++) {
	printf("  %13.5e",eigvals[i])    ;
	fprintf(o_log,"%e\t",eigvals[i]) ;   }
      printf("\n")        ;
      fprintf(o_log,"\n") ;
      fflush(o_log)       ;
    }

  /* Posterior covariance matrix */
  o_cov    = openfile_wt(opt->o_covf) ;
  for(i=0 ; i<*M ; i++)
    {
      for(k=0 ; k<*M ; k++)
	for(j=0 ; j<nk ; j++)
	  cov[i][k] += eigvects[i][j]*eigvects[k][j]/eigvals[j] ;
      for(k=0 ; k<*M ; k++)
	fprintf( o_cov,"%16.2f ", cov[i][k]) ;
      fprintf( o_cov,"\n") ;
    }
  fclose(o_cov) ;

  /* Least-squares inversion */
  if (*M == NM-1)      /* Null trace     */
    {
      for( s=0 ; s<*nsac ; s++)
	{
	  N = hd_synt[s].npts ;
	  for(i=0 ; i<*M ; i++)
	    {
	      for(k=0 ; k<2 ; k++)
		for(l=0 ; l<N ; l++)
		  vma[i] += cov[i][k]*(G[s][k][l]-G[s][2][l])*opt->wgt[s]*d[s][l] ;
	      for(k=2 ; k<*M ; k++)
		for(l=0 ; l<N ; l++)
		  vma[i] += cov[i][k]*G[s][k+1][l]*opt->wgt[s]*d[s][l] ;
	    }
	}
      vma[5] = vma[4] ;
      vma[4] = vma[3] ;
      vma[3] = vma[2] ;
      vma[2] = -vma[0] -vma[1] ;
    }
  else if (*M == NM) /* No constraints */
    for( s=0 ; s<*nsac ; s++)
      for(i=0 ; i<*M ; i++)
	for(k=0 ; k<*M ; k++)
	  for(l=0 ; l<hd_synt[s].npts ; l++)
	    vma[i] += cov[i][k]*G[s][k][l]*opt->wgt[s]*d[s][l] ;
  else {
      fprintf(stderr,"ERROR : bad nb of parameters (M=%d)\n",*M);
      exit(1);  }

  /* Memory Freeing */
  free((void*)eigvals) ;
  for (i=0 ; i<*M ; i++) 
    { 
      free((void*)GtG[i])      ; 
      free((void*)eigvects[i]) ; 
      free((void*)cov[i])      ;
    }
  free((void**)GtG)      ; 
  free((void**)eigvects) ; 
  free((void**)cov)      ; 
}
  
/*********************************************************/
/*          residual_moment(vm, ma, mb, mc)              */
/*********************************************************/
/* Compute residual between moment tensor vm[0] and vm[1]*/
/****                                                *****/
/* Input: vm pointer on mom. tensor elements array (1x6) */
/* Output: ma : scalar moment of vm[0] (Dahlen & Tromp)  */
/*         mb : scalar moment of vm[1] (Dahlen & Tromp)  */
/*         mc : rms deviation of normalized mom. tensor  */
void 
residual_moment(vm, ma, mb, mc)
     double **vm, *ma, *mb, *mc ;
{
  int    i;
  double *vmc;
  
  vmc = double_alloc(NM) ;

  mt2sm(vm[0], ma) ;
  mt2sm(vm[1], mb) ;
  for(i=0 ; i<NM ; i++) 
    {
      vmc[i] = vm[1][i]/(*mb) - vm[0][i]/(*ma); 
    }
  mt2sm(vmc, mc) ;
  free((void*)vmc) ;
}

/******************************************/
/*               mt2sm(vm,sm)             */
/******************************************/
/* scalar moment (def. of Dahlen & Tromp) */
/* Input: vm : moment tensor              */
/* Output: sm : scalar moment             */
void 
mt2sm(double *vm,double *sm) 
{
  *sm = vm[0]*vm[0]+vm[1]*vm[1]+vm[2]*vm[2] ;
  *sm += 2.*(vm[3]*vm[3]+vm[4]*vm[4]+vm[5]*vm[5]) ;
  *sm = sqrt(*sm/2.) ;
}


double **
set_mt(double *vm)
{
  double **TM ;
  TM       = double_alloc2(3,3) ;
  TM[0][0] =   vm[1] ;  /* Rotation of pi around East-axis       */
  TM[1][1] =   vm[2] ;  /* The new system is (North, East, Down) */
  TM[2][2] =   vm[0] ;  /* (Aki's coordinates)                   */
  TM[0][1] =  -vm[5] ;
  TM[0][2] =   vm[3] ;
  TM[1][2] =  -vm[4] ;
  TM[1][0] = TM[0][1] ;
  TM[2][0] = TM[0][2] ;
  TM[2][1] = TM[1][2] ;
  return TM ;
}
  
void 
get_planes(vm, TM, eval3, s1,d1,r1, s2,d2,r2)
     double *vm, ***TM, **eval3, *s1, *d1, *r1, *s2, *d2, *r2;
{
  int    nrot, i ;
  double si, co  ;
  double **evec3, *vn, *vs ;

  
  /* Memory allocation */
  *eval3 = double_alloc(3)    ;
  evec3  = double_alloc2(3,3) ;
  vn     = double_alloc(3)    ;
  vs     = double_alloc(3)    ;
  
  /* Tensor representation */
  (*TM) = set_mt(vm) ;

  /* Get eigvalues and eigvectors*/
  jacobi((*TM),3,3,(*eval3),evec3,&nrot) ;
  eigsrt((*eval3),evec3,3) ;

  /* Check if eigenvectors are upwards */
  if(evec3[2][0] < 0.0) 
    for(i=0 ; i<3 ; i++)
      {
	evec3[i][0] = -evec3[i][0] ;
      }

  if(evec3[2][2] < 0.) 
    for(i=0 ; i<3 ; i++)
      evec3[i][2] = -evec3[i][2] ;

  /* Cross-product v2 = v1 x v3 */
  evec3[0][1] = evec3[1][0]*evec3[2][2]-evec3[1][2]*evec3[2][0] ;
  evec3[1][1] = evec3[2][0]*evec3[0][2]-evec3[2][2]*evec3[0][0] ;
  evec3[2][1] = evec3[0][0]*evec3[1][2]-evec3[0][2]*evec3[1][0] ;
     
  (*TM)[0][1] = (*TM)[1][0] ; /* useless */
  (*TM)[0][2] = (*TM)[2][0] ;
  (*TM)[1][2] = (*TM)[2][1] ;
  
  /* *** First nodal plane *** */
  for(i=0 ; i<3 ; i++)
    {
      vn[i] = (evec3[i][0]+evec3[i][2])/sqrt(2.) ;
      vs[i] = (evec3[i][0]-evec3[i][2])/sqrt(2.) ;
    }
  if (vn[2] > 0.)
    for(i=0; i<3; i++)
      {
	vn[i] = -vn[i] ;
	vs[i] = -vs[i] ;
      }
  
  *s1 = atan2(-vn[0], vn[1]) ;
  *d1 = acos(-vn[2]) ;
  si = sin(*s1) ;
  co = cos(*s1) ;
  *r1 = atan2((vs[0]*si - vs[1]*co),-(vs[0]*co + vs[1]*si)*vn[2]) ;

  /* *** Second nodal plane *** */
  for(i=0 ; i<3 ; i++)
    {
      vn[i] = (evec3[i][0]-evec3[i][2])/sqrt(2.) ;
      vs[i] = (evec3[i][0]+evec3[i][2])/sqrt(2.) ;
    }

  if (vn[2] > 0.) 
    for(i=0; i<3; i++)
      {
	vn[i] = -vn[i] ;
	vs[i] = -vs[i] ;
      }
  *s2 = atan2(-vn[0], vn[1]) ; /* strike */
  *d2 = acos(-vn[2]) ;         /* dip    */
  si = sin(*s2) ;             
  co = cos(*s2) ;
  *r2 = atan2((vs[0]*si - vs[1]*co),-(vs[0]*co + vs[1]*si)*vn[2]); /* rake */

  *s1 = (*s1)/((double)DEG2RAD) ;
  if ((*s1) < 0.) (*s1) += 360. ;
  *d1 = (*d1)/((double)DEG2RAD) ;
  *r1 = (*r1)/((double)DEG2RAD) ;
  
  *s2 = (*s2)/((double)DEG2RAD) ;
  if ((*s2) < 0.) (*s2) += 360. ;
  *d2 = (*d2)/((double)DEG2RAD) ;
  *r2 = (*r2)/((double)DEG2RAD) ;
  
  /* Memory Freeing */
  free((void*)vn) ;
  free((void*)vs) ;
  for(i=0 ; i<3 ; i++) 
    free((void*)evec3[i]) ;
  free((void**)evec3) ;
}

void
set_wgt(int ns, sachdr *hd_data,structopt *opt) 
{
  opt->wgt[ns] = 0. ;
  if (!strncmp(hd_data->kcmpnm, "LHZ",3))
    opt->wgt[ns] = opt->wZ;
  else if (!strncmp(hd_data->kcmpnm, "LHN",3))
    opt->wgt[ns] = opt->wN;
  else if (!strncmp(hd_data->kcmpnm, "LHE",3))
    opt->wgt[ns] = opt->wE;
}

void
azpond (sachdr *hd_synt, int ns, double *wgt) 
{
  int    i, j                          ;
  double moy, std, mok, mincov, maxwgt ;
  double onepisq, twopisq, twopidaz    ;
  double *aztab, *azcov                ;

  aztab = double_calloc(ns) ;
  azcov = double_calloc(ns) ;

  for (i=0; i<ns; i++)
    aztab[i] = hd_synt[i].az ;
  sort(aztab, &ns);
  for (i=0; i<ns-1; i++)
    {
      aztab[i] = aztab[i+1]-aztab[i] ;
      moy  += aztab[i] ;
    }
  moy /= ns ;
  for (i=0; i<ns-1; i++)
    std += (aztab[i] - moy)*(aztab[i] - moy) ;
  std = sqrt(std/ns) ;
  if (std == 0. )
    std = 0.0000001 ;
  
  mincov  = 1.e10  ;
  onepisq = 32400  ;
  twopisq = 129600 ;
  for(i=0; i<ns ; i++)
    {
      for(j=0; j<ns; j++)
	{
	  mok = (hd_synt[i].az-hd_synt[j].az) ;
	  if (mok > onepisq)
	    {
	      mok += twopisq ;
	      twopidaz  = 720.*(hd_synt[i].az-hd_synt[j].az) ;
	      if (twopidaz >= 0)
		mok -= twopidaz ;
	      else
		mok += twopidaz ;
	    }
	  azcov[i] += exp(-mok/std) ;
	}
      if (azcov[i] < mincov)
	mincov = azcov[i] ;
    }
  maxwgt = 0. ;
  for(i=0; i<ns; i++)
    {
      wgt[i] += mincov/azcov[i] ;
      if (wgt[i] > maxwgt)
	maxwgt = wgt[i]         ;
    }
  for(i=0; i<ns; i++)
    wgt[i] /= maxwgt ;
  printf("\n") ; 
  free((void*)aztab) ;
  free((void*)azcov) ;
}

void
free_G(double ***G)
{
  int i ;
  for(i=0; i<NM; i++)
    free((void *) (*G)[i]) ;
  free((void**) (*G)) ;
}

int
fill_G(char *gf_file, char *datafile, sachdr *hd_GF, sachdr *hd_data, int npts, 
       double Ptt, double twp_beg, double twp_end, double *buffer, double *G, 
       structopt *opt, FILE *o_log)
{
  int i, ierror = 1 ;
  int n1_GF, n2_GF  ;
  double t0         ;
  double *g = &G[0] ;

  /* Read Header */
  rhdrsac( gf_file, hd_GF, &ierror) ; 
  if (hd_GF->delta != hd_data->delta) 
    {
      fprintf( o_log,"**** Incorrect sampling period, rejected trace : %s\n", datafile) ;
      fprintf(stderr,"**** Incorrect sampling period, rejected trace : %s\n", datafile) ;
      return 1 ;
    }
  /* GF Time Window */
  t0 = Ptt + (double)hd_GF->o ;
  n1_GF = (int)((t0 + twp_beg - (double)hd_GF->b - opt->dts_val)  / ((double)hd_GF->delta)) ; /* first GF Sample (corrected) */
  n2_GF = n1_GF + (int)((twp_end - twp_beg) / ((double)hd_GF->delta)) ;                       /* Last GF Sample */
  if ( n2_GF >= hd_GF->npts ) /* GF Rejected */
    {
      fprintf(o_log,"stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", 
	      hd_data->kstnm, hd_data->knetwk, hd_data->kcmpnm, 
	      hd_data->gcarc, hd_data->az, hd_data->user[2], hd_data->user[3]) ;
      fprintf( o_log,"**** Incomplete GF, rejected : %s\n"  , gf_file)         ;  
      fprintf(stderr,"**** Incomplete GF, rejected : %s ", gf_file)            ; 
      fprintf(stderr,"((n2_GF=%d)>=(npts=%d))\n", n2_GF, hd_GF->npts)          ;
      return 1 ;
    }
  else if ( n1_GF < 0 )       /* Fill negative samples with zeros */
    {
      fprintf(o_log,"stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", 
	      hd_data->kstnm, hd_data->knetwk, hd_data->kcmpnm, 
	      hd_data->gcarc, hd_data->az, hd_data->user[2], hd_data->user[3]) ;
      fprintf( o_log,"**** Incomplete GF, filling with zeros : %s\n", gf_file) ;  
      fprintf(stderr,"**** Incomplete GF,")                                    ;
      fprintf(stderr,"filling with zeros : %s (n1_GF=%d)\n", gf_file, n1_GF )  ;
      for(i=n1_GF; i<0; i++)
	G[i-n1_GF] = 0. ;
      g = &G[-n1_GF] ;
      npts += n1_GF  ;
      hd_GF->b += n1_GF*hd_GF->delta ;
      n1_GF = 0 ;
    }
  /* Read GF samples */
  hd_GF->b += hd_GF->delta * (float)n1_GF  ;
  hd_GF->npts = n2_GF + 1                  ;
  rdatsac(gf_file, hd_GF, buffer, &ierror) ;
  memcpy (g,buffer+n1_GF,npts * sizeof(double));
  return 0;
}
/* Corrections, Notes, etc.:            */
/* hd_GF->b += (n1_GF-1)*hd_GF->delta ; */


void 
set_matrices (i_saclst, evdp, wp_win4, nsac, nsini, sacfiles, hd_synt,
	      data, G, opt, eq,o_log) 
     int    *nsac, *nsini;
     double *evdp, *wp_win4, ***data, ****G ;
     char   *i_saclst, ***sacfiles          ;  
     FILE    *o_log   ;
     sachdr **hd_synt ;
     structopt *opt   ;
     str_quake_params *eq;
{
  int    i, j, flag, flag2,ierror, npts  ;
  int    n1_data, n2_data  ;
  int    ns, nh = NDEPTHS, nd = NDISTAS     ;
  double gcarc, t0, Ptt, twp_beg, twp_end   ;
  double *tmparray, *dv, *tv                ;
  char   *datafile, *gf_file                ;
  char   *dum,*buf, *GF                     ; 
  char   gfdirs[6][7] = {"gf_rr/","gf_tt/", 
			 "gf_pp/","gf_rt/", 
			 "gf_rp/","gf_tp/"} ;
  FILE   *i_sac                             ;
  sachdr hd_data, hd_GF                     ;


  /* Opening data file list */
  i_sac = openfile_rt(i_saclst,nsini) ;
  *nsac = *nsini ;

  /* Allocating memory */
  dum      = char_alloc(32)    ;
  buf      = char_alloc(FSIZE) ;
  gf_file  = char_alloc(FSIZE) ;
  datafile = char_alloc(FSIZE) ;
  dv       = double_alloc(nd)  ;
  tv       = double_alloc(nd)  ;
  tmparray = double_alloc((int)__LEN_SIG__) ;  
  hdr_alloc(&hd_data) ;
  hdr_alloc(&hd_GF)   ;
  flag2 = 0 ;
  if ( *data == NULL    && *G == NULL       &&  opt->wgt == NULL &&   \
       opt->rms_in == NULL && opt->p2p == NULL && opt->avg == NULL && \
       *sacfiles == NULL ) 
    {
      *data     = double_alloc2p( *nsac )  ;
      *G        = double_alloc3p( *nsac )  ;
      opt->wgt  = double_alloc( *nsac )    ;
      *sacfiles = char_alloc2(*nsac,FSIZE) ;
      opt->rms_in = double_alloc(*nsac)    ;
      opt->p2p    = double_alloc(*nsac)    ;
      opt->avg    = double_alloc(*nsac)    ;
      hdr_tab(hd_synt, *nsac)              ; 
      flag2 = 1 ;
    }
  
  /* Set travel times */
  ierror = 1 ; /* error flag */
  trav_time_init(&nh, &nd, evdp, dv, tv, &ierror) ;  

  /* Loop on stations */
  ns = 0 ;
  opt->dmin = 2.e4 ;
  opt->dmax = 0.   ;  
  for(i=0; i<*nsac; i++)
    {
      /* Read data file list */
      if ( opt->op_pa <= 0 && opt->th_val <= 0)
	{ 
	  flag = fscanf (i_sac, "%s", datafile) ;
	  fgets(buf,LSIZE,i_sac); /* end of line */
	  check_scan(1, flag, i_saclst, i_sac)  ;
	}
      else 
	{
	  flag = fscanf (i_sac, "%s %f %f %s %s %s %s %s %lf %lf %lf %s",
			 datafile,&(*hd_synt)[ns].az,&(*hd_synt)[ns].gcarc,dum,dum,dum,dum,dum,
			 &opt->rms_in[ns],&opt->p2p[ns],&opt->avg[ns],dum) ;
	  strcpy(buf,i_saclst);
	  strcat(buf," (nb of columns may be incorrect)");
	  check_scan(12, flag, buf, i_sac) ;
	} 
      
      /* Read data header and weights */
      rhdrsac(datafile,  &hd_data, &ierror);
      set_wgt(ns, &hd_data, opt) ;
      if (opt->wgt[ns] <= 0.)
	{
	  fprintf(stderr,"**** null weight, rejected : %s\n", datafile) ;
	  continue ;
	}
      if ((float)(*evdp) != hd_data.evdp)	
	{
	  fprintf(stderr,"WARNING : depth %f in sac header is different of %f from pde in CMTFILE\n",(float)*evdp, hd_data.evdp);
	  fprintf(stderr," ...you should carefully re-check gcarc and evdp header variables in file %s\n",datafile); 
	}
       
      /* Data Time Window  */
      gcarc = (double) hd_data.gcarc ;
      if (opt->dmin > gcarc)
	opt->dmin = gcarc ;
      if (opt->dmax < gcarc)
	opt->dmax = gcarc ;
      fflush(stdout);
      trav_time(&gcarc, tv, dv, &nd, &Ptt, &ierror) ;
      wp_time_window(&gcarc, wp_win4, &twp_beg, &twp_end) ;
      t0 = Ptt + hd_data.o ; 
      n1_data = (int)((t0 + twp_beg - (double)hd_data.b) / ((double)hd_data.delta)) ; /* first data Sample (corrected)  */
      n2_data = n1_data + (int)((twp_end - twp_beg) / ((double)hd_data.delta))      ; /* Last data Sample               */
      npts    = n2_data - n1_data + 1 ;
      if ( (n1_data < 0) || (n2_data >= hd_data.npts) ) 
	{
	  fprintf(stderr,"**** Incomplete data, rejected : %s\n", datafile) ;
	  fprintf( o_log,"**** Incomplete data, rejected : %s\n", datafile) ;
	  continue ;                                      
	}
      
      /* Loop on GF components */
      flag = 0 ;
      if (flag2)
	(*G)[ns] = double_alloc2(NM,npts) ;
      for(j=0; j<NM; j++)  
	{
	  /* GF filename */
	  strcpy(gf_file,eq->gf_dir);
	  GF = get_gf_filename(gfdirs[j], hd_data.kstnm, hd_data.knetwk, hd_data.kcmpnm, ".SAC.sac.bp") ;
	  strcat(gf_file,GF);
	  free((void*) GF);
	  flag = fill_G(gf_file, datafile, &hd_GF, &hd_data, npts, Ptt, twp_beg, twp_end, 
			tmparray, (*G)[ns][j], opt, o_log) ;
	  if (flag)
	    {
	      free_G((*G)+ns) ;
	      break ;
	    }
	}
      if (flag) /* Error reading GF */
	continue ;
      
      
      /* Read data samples */
      hd_data.npts = n2_data + 1 ;
      rdatsac(datafile, &hd_data, tmparray, &ierror)   ;  
      if (flag2)
	(*data)[ns] = double_alloc(npts)               ;
      memcpy ((*data)[ns],tmparray+n1_data,npts * sizeof(double)) ;
      
      /* set synthetic header */
      (*hd_synt)[ns].delta  = hd_data.delta  ;
      (*hd_synt)[ns].npts   = npts           ;
      (*hd_synt)[ns].nzyear = hd_data.nzyear ;
      (*hd_synt)[ns].nzjday = hd_data.nzjday ;
      (*hd_synt)[ns].nzhour = hd_data.nzhour ;
      (*hd_synt)[ns].nzmin  = hd_data.nzmin  ;
      (*hd_synt)[ns].nzsec  = hd_data.nzsec  ;
      (*hd_synt)[ns].nzmsec = hd_data.nzmsec ;
      (*hd_synt)[ns].o      = hd_data.o      ;
      (*hd_synt)[ns].b      = hd_GF.b + hd_data.o - hd_GF.o ;
      
      /********************************  MAY BE MODIFIED *************************************/
      (*hd_synt)[ns].user[0] = (float)(n1_data + 1) ;
      (*hd_synt)[ns].user[1] = (float)(n2_data + 1) ;
      (*hd_synt)[ns].user[2] = twp_beg ;
      (*hd_synt)[ns].user[3] = twp_end ;     
      (*hd_synt)[ns].user[4] = (*hd_synt)[ns].b - hd_data.b - (n1_data - 1)*hd_data.delta ;
      /***************************************************************************************/

      strcpy((*hd_synt)[ns].kstnm, hd_data.kstnm)   ; 
      strcpy((*hd_synt)[ns].knetwk, hd_data.knetwk) ; 
      strcpy((*hd_synt)[ns].kcmpnm, hd_data.kcmpnm) ; 

      
      if (opt->op_pa <= 0.) 
	{
	  (*hd_synt)[ns].az    = hd_data.az    ;
	  (*hd_synt)[ns].gcarc = hd_data.gcarc ;
	  calc_stat( npts, (*data)[ns], &opt->p2p[ns], &opt->avg[ns]);  
	  opt->p2p[ns] *= 1000                 ;
	  opt->avg[ns] *= 1000                 ; 
	}

      if (opt->th_val <= 0. && opt->med_val <= 0.)
	fprintf( o_log,"stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", (*hd_synt)[ns].kstnm, 
		 (*hd_synt)[ns].knetwk, (*hd_synt)[ns].kcmpnm, (*hd_synt)[ns].gcarc, (*hd_synt)[ns].az, 
		 (*hd_synt)[ns].user[2], (*hd_synt)[ns].user[3]) ;
      strcpy( (*sacfiles)[ns], datafile) ;
      ns++ ;
    }
  *nsac  = ns   ;
  fclose(i_sac) ;
  /* Memory Freeing */
  free((void*)tmparray) ;
  free((void*)dv)       ;
  free((void*)tv)       ;
  free((void*)dum)      ;
  free((void*)buf)      ;
  free((void*)datafile) ;
  free((void*)gf_file)  ;
}
/* Corrections, Notes, etc.: */
// n1_data = (int)((t0 + twp_beg - (double)hd_data.b) / ((double)hd_data.delta)) - 1 ;           /* Error **** first data Sample  */
// n1_GF = (int)((t0 + twp_beg - (double)hd_GF.b - opt->dts_val)  / ((double)hd_GF.delta)) - 1 ; /* Error **** first GF Sample    */


void 
realloc_gridsearch(nsac, rms, global_rms, vm, dcalc, flag)
     int nsac, flag ;
     double ***rms, **global_rms, **vm, ****dcalc ;
{
  int i, j ; 

  free((void*)(*vm))  ;
  free((void*)(*rms)) ;
  free((void*)(*global_rms)) ;
  for(i=0 ; i<nsac ; i++)
    {
      for(j=0 ; j<flag ; j++)
	free((void*)(*dcalc)[i][j]) ;
      free((void**)(*dcalc)[i]);
    }
  *vm  = double_calloc(NM) ;
  *rms = double_calloc2(nsac, 2*flag) ;
  *global_rms = double_calloc(2*flag) ;
}

void
fast_ts_gridsearch(nsac, M, sacfiles, hd_synt, data, G, dcalc, rms, global_rms, opt, eq, 
		   rmsini, tsopt, rmsopt, o_log)
     int    nsac, M    ;
     char   **sacfiles ;
     sachdr *hd_synt   ;
     double **data, ***G, ****dcalc, ***rms         ;
     double **global_rms, *rmsini, *tsopt, *rmsopt ;
     structopt        *opt ;
     str_quake_params *eq  ; 
     FILE *o_log ;
{
  int    it, k, Nexp, flag, nsini;
  double Err, Cond, dt, ts ,dtmin, dtmax,tsopt2,rmsopt2; 
  FILE   *tmp, *o_gs   ;

  printf("FAST CENTROID TIME DELAY GRID SEARCH\n") ;
  fprintf(o_log,"Time-shift grid search is enabled (output : %s)\n",opt->gsfile) ;
  
  /* Initialize variables */
  dtmin = opt->dts_min - eq->ts ;
  dtmax = opt->dts_max - eq->ts ;
  dt    = opt->dts_step ;

  it   = 0 ;
  k    = 0 ;
  Nexp = 0 ;
  flag = 1 ;
  *tsopt  = 0. ;
  tsopt2  = dtmin  ;
  *rmsini = 1000.*(*global_rms)[0] ;
  *rmsopt = *rmsini ;
  rmsopt2 = 1.1e10 ;

  tmp     = openfile_wt("_tmp_ts_table") ;
  while( it < opt->Nit )
    {
      opt->dts_val = dtmin ;
      
      while ( opt->dts_val <= dtmax )
	{
	  /* Free memory */
	  realloc_gridsearch(nsac, rms, global_rms, &eq->vm[0], dcalc, flag) ;
	  /* Compute inversion for opt->dts_val */
	  set_matrices (opt->i_saclst, &eq->pde_evdp, eq->wp_win4, &nsac, &nsini,
			&sacfiles, &hd_synt, &data, &G, opt, eq, o_log) ;    
	  inversion(&M, &nsac, hd_synt, G, data, eq->vm[0], &Cond, opt, o_log) ;
	  calc_data(&nsac, hd_synt, G, eq->vm, data, (*dcalc), opt, flag)    ;   
	  calc_rms(&nsac, hd_synt, data, (*dcalc), (*rms), (*global_rms), opt, flag)    ;

	  /* Get RMS error */
	  Err = 1000.*(*global_rms)[0] ;
	  ts  = eq->ts+opt->dts_val    ;
	  fprintf( tmp,"%02d %8.2f %8.2f %8.2f %8.2f %12.8f %12.8f\n",k,ts,
		   eq->evla,eq->evlo,eq->evdp,Err,(*global_rms)[0]/(*global_rms)[1]);
	  printf("        ts = %5.1f rms = %12.7f mm\n",ts, Err) ;
	  if (Err < *rmsopt)
	    {
	      *rmsopt = Err ;
	      *tsopt   = opt->dts_val ;
	    }
	  else if (Err < rmsopt2)
	    {
	      rmsopt2 = Err ;
	      tsopt2  = opt->dts_val ;
	    }
	  opt->dts_val += dt ;
	  k++ ;
	}
      printf("Optimum values: time_shift = %5.1f rms = %12.7f mm \n",eq->ts+*tsopt,*rmsopt) ;
      if (dtmax <= *tsopt + dt && Nexp < 5)
	{
           printf("Optimum value on the maximum explored time-shift\n   ... extending the time-shift grid-search area\n");
	   dtmin = *tsopt + dt   ;
	   dtmax = *tsopt + 4*dt ;
	   Nexp++;
	   continue ;
	}
      if (it>0)
	dt = dt/2. ;
      if (tsopt2 <= *tsopt)
	{
	  dtmin = tsopt2 - dt/2. ;
	  dtmax = (*tsopt)  + dt/2. ;
	}
      else if (tsopt2 > *tsopt)
	{
	  dtmin = (*tsopt)  - dt/2. ;
	  dtmax = tsopt2 + dt/2. ;
	}
      if (dtmin < 1. - eq->ts)
	while (dtmin <  1. - eq->ts)
	  dtmin += dt ;
      it++;
    }
  fclose(tmp);
  /* Write output file */
  o_gs = openfile_wt(opt->gsfile) ;
  fprintf( o_gs,"%5.1f %12.8f\n",eq->ts+*tsopt,*rmsopt);
  fprintf( o_gs,"%5.1f %12.8f\n",eq->ts,*rmsini);
  tmp  = openfile_rt("_tmp_ts_table", &k) ;
  while ((k = getc(tmp)) != EOF)
    putc(k,o_gs) ;
  fclose( tmp) ;
  fclose(o_gs) ;
}

void 
screen_rms(nsac, data_name, data, G, hd_synt, opt, o_log)
     int    *nsac        ;
     double **data, ***G ;
     char   **data_name  ;
     sachdr *hd_synt     ; 
     structopt *opt      ;
     FILE *o_log         ;
{
  int j, newn;

  newn = 0 ;
  for (j=0; j<*nsac; j++)
    {
      if( opt->rms_in[j] < opt->th_val )
	{
	  fprintf( o_log,"%-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", hd_synt[j].kstnm, 
		   hd_synt[j].knetwk, hd_synt[j].kcmpnm, hd_synt[j].gcarc, hd_synt[j].az, 
		   hd_synt[j].user[2], hd_synt[j].user[3]) ; 
	  data_name[newn] = data_name[j] ;
	  data[newn]      = data[j]     ;
	  G[newn]         = G[j]        ;
	  hd_synt[newn]   = hd_synt[j]  ;
	  opt->p2p[newn]  = opt->p2p[j] ;
	  opt->avg[newn]  = opt->avg[j] ;
	  opt->wgt[newn]  = opt->wgt[j] ;
	  newn++ ;	  
	}
      else
	{
	  fprintf(stderr, "**** Rejected trace (rms exceed the threshold): %s\n", 
		  data_name[j]) ;
	  free((void*)data_name[j]) ;
	  free((void*)data[j])      ;
	  free_G(G+j)               ;
	}
    }
  *nsac = newn ;
}


void 
screen_med(nsac, data_name, data, G, hd_synt, opt, o_log)
     int    *nsac        ;
     double **data, ***G ;
     char   **data_name  ;  
     sachdr *hd_synt     ; 
     structopt *opt      ;
     FILE *o_log         ;
{
  int    j, newn ;
  double min, max, val ;
  
  min = 0.1 * (opt->p2p_med) ;
  max = 3.0 * (opt->p2p_med) ;
  
  fprintf(o_log,"screen_med:\n") ;
  fprintf(o_log,"   p2p_med: %15.8f\n",opt->p2p_med) ;
  fprintf(o_log,"   reject p2p < : %15.8f or > %15.8f\n",min,max) ;
  fprintf(o_log,"   reject avg > : %15.8f \n",opt->p2p_med/2) ;
  newn = 0 ;
  for (j=0;j<*nsac;j++)
    {
      val = opt->p2p[j];
      if ( (min < val) && (val < max) && (fabs(opt->avg[j]) < (opt->p2p_med)/2.) )
	{
	  if (opt->th_val <= 0.)
	    fprintf( o_log,"%-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", hd_synt[j].kstnm, 
		     hd_synt[j].knetwk, hd_synt[j].kcmpnm, hd_synt[j].gcarc, hd_synt[j].az, 
		     hd_synt[j].user[2], hd_synt[j].user[3]) ; 
	  data_name[newn] = data_name[j] ;
	  data[newn]     = data[j]     ;
	  G[newn]        = G[j]        ;
	  hd_synt[newn]  = hd_synt[j]  ;
	  opt->p2p[newn] = opt->p2p[j] ;
	  opt->avg[newn] = opt->avg[j] ;
	  opt->wgt[newn] = opt->wgt[j] ;
	  newn++ ;
	}
      else 
	{
	  fprintf(stderr,"**** Rejected trace (p2p or avg out of bounds): %s\n",
		  data_name[j])    ; 
	  free((void*)data_name[j]);
	  free((void*)data[j])     ;
	  free_G(G+j)              ;
	}
    }
 *nsac = newn ;
}


void 
median(nsac, opt)
     int       *nsac ; 
     structopt *opt  ;
{
  int    i    ;
  double *tmp ;
  
  tmp = double_alloc( *nsac ) ;
  for (i=0;i<*nsac;i++)
    tmp[i] = opt->p2p[i]    ;

  sort(tmp, nsac);

  opt->p2p_med  = tmp[(int)    ((*nsac)/2.)-1] ;
  opt->p2p_low  = tmp[(int)    ((*nsac)/4.)-1] ;
  opt->p2p_high = tmp[(int) ((*nsac)*3./4.)-1] ;

  free((void *)tmp);
}
  

void 
sort(tab, ns)
     int    *ns  ;
     double *tab ; 
{
  int    i, j ;
  double swp  ;

  for (i=0;i<*ns;i++)
      for(j=i+1;j<*ns;j++)
	{
	  if (tab[i] > tab[j]) {
	    swp    = tab[i] ;
	    tab[i] = tab[j] ;
	    tab[j] = swp ;     }
	}
}


void 
w_log_header(argv, opt, eq, wp_win4, o_log) 
     double *wp_win4 ;
     char   **argv   ; 
     FILE   *o_log   ;
     str_quake_params  *eq      ;
     structopt *opt  ;   
{
  int    i ;
  time_t t ;
  
  t = time(NULL);
  fprintf(o_log,"++++++++++\nDate-Time(UTC): %s\n",asctime(gmtime(&t)));  
  fprintf(o_log,"Command_line: %s\n",argv[0])     ;
  fprintf(o_log,"Ref_solution: %s\n",eq->cmtfile) ;
  fprintf(o_log,"rms_treshold: %f\ncond_thre: %f\ndamp_fac: %f\nmed_val: %f\n",
	 opt->th_val, opt->cth_val, opt->df_val, opt->med_val) ;  
  fprintf(o_log,"R_event_id: %s\nR_time_shift: %8.2f\nR_half_duration: %8.2f\n"
	  ,eq->evid, eq->ts, eq->hd)   ;
  fprintf(o_log,"R_latitude: %8.2f\nR_longitude: %8.2f\nR_depth: %8.2f\n",
	  eq->evla, eq->evlo,eq->evdp) ;
  fprintf(o_log,"R_Moment_tensor: ")   ;
  for (i=0; i<NM; i++)
    fprintf(o_log,"%12.4e ",eq->vm[1][i]);
  fprintf(o_log,"\n")       ;
  fprintf(o_log,"WP_WIN: ") ;
  for (i=0; i<4; i++)
    fprintf(o_log,"%7.1f ",wp_win4[i]) ;
  fprintf(o_log,"\n") ;  
  fflush(o_log)       ;
}

void get_param2(file, opt, eq,  flag)
     int    *flag   ;
     char   *file   ;
     structopt *opt ;
     str_quake_params *eq ;
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

  get_i_master(file, keys, nimas, eq) ;  
  *flag = get_cmtf(eq, *flag) ;
  if (*flag ==1)
    opt->ref_val = 0. ;    
  
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
  
  fprintf(stderr,"\nTime-shift grid-search :\n");
  fprintf(stderr,"  -dts real_value         apply a time-shift to Green's functions (no time-shift)\n");
  fprintf(stderr,"  -ts tsmin dts tsmax     fast time-shift grid-search (no grid-search)\n") ;
  fprintf(stderr,"  -Nit                    nb. of iteration for time_shift grid-search (%d)\n",opt->Nit) ;
  fprintf(stderr,"  -ogsf outputfile        grid-search output filename (%s)\n",opt->gsfile) ; 
  
  
  fprintf(stderr,"\n  -h, --help              display this help and exit\n\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>\n") ;
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
  strcpy(opt->i_saclst,"i_wpinversion") ;
  strcpy(opt->o_saclst,"o_wpinversion") ;
  strcpy(opt->log,"wpinversion.log")    ;
  strcpy(opt->o_covf,"o_covariance")    ;
  strcpy(opt->o_cmtf,"WCMTSOLUTION")    ;
  strcpy(opt->p_data,"fort.15")  ;
  strcpy(opt->psfile,"p_wpinversion")   ;
  strcpy(opt->wpbmfile,"")  ;
  strcpy(opt->refbmfile,"") ;
  strcpy(opt->osacdir,"./") ;
  strcpy(opt->gsfile,"grid_search_ts_out") ;
  eq->cmtfile[0] = '\0' ;
  opt->th_val    = 0. ;
  opt->cth_val   = 0. ;
  opt->df_val    = 0. ;
  opt->med_val   = 0. ;
  opt->op_pa     = 0. ;
  opt->ref_val   = 1. ;
  opt->ntr_val   = 1. ;
  opt->wZ        = 1. ;
  opt->wN        = 1. ;
  opt->wE        = 1. ;
  opt->dts_val   = 0. ;
  opt->dts_min   = 0. ;
  opt->dts_max   = 0. ;
  opt->dts_step  = 0. ;
  opt->azp       = 0. ;

  opt->Nit       = 1 ;

  k = 0 ;
  for( i = 0; i<numarg2; i++ )
    {
      j = i + numarg1 + 1 ;
      fflush(stdout);
      if (!strncmp(argv[j],"-imas",5)){
	get_char_arg(argv, j, i, numarg2, opt->i_master) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-ifil",5)){
	get_char_arg(argv, j, i, numarg2, opt->i_saclst) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-ofil",5)){
	get_char_arg(argv, j, i, numarg2, opt->o_saclst) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-log",4)){
	get_char_arg(argv, j, i, numarg2, opt->log) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-ocovf",6)){
	get_char_arg(argv, j, i, numarg2, opt->o_covf) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-ocmtf",6)){
	get_char_arg(argv, j, i, numarg2, opt->o_cmtf) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-pdata",6)){
	get_char_arg(argv, j, i, numarg2, opt->p_data) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-icmtf",6)){
	get_char_arg(argv, j, i, numarg2, eq->cmtfile) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-gfdir",6)){
	get_char_arg(argv, j, i, numarg2, eq->gf_dir) ;
	add_slash(eq->gf_dir);
	k+=2 ;}   
      if (!strncmp(argv[j],"-osyndir",8)){
	get_char_arg(argv, j, i, numarg2, opt->osacdir) ;
	k+=2 ;}      
      if (!strncmp(argv[j],"-ps",6)){
	get_char_arg(argv, j, i, numarg2, opt->psfile) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-wpbm",6)){
	get_char_arg(argv, j, i, numarg2, opt->wpbmfile) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-refbm",6)){
	get_char_arg(argv, j, i, numarg2, opt->refbmfile) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-th",3)){
	get_num_arg(argv, j, i, numarg2,"%lf",&opt->th_val) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-cth",4)){
	get_num_arg(argv, j, i, numarg2,"%lf",&opt->cth_val) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-df",3)){
	get_num_arg(argv, j, i, numarg2,"%lf",&opt->df_val) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-dts",4)){
	get_num_arg(argv, j, i, numarg2,"%lf",&opt->dts_val);
	k+=2 ;}
      if (!strncmp(argv[j],"-ts",3))
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
      if (!strncmp(argv[j],"-Nit",4)) 
	{
	  get_num_arg(argv, j, i, numarg2,"%d", &opt->Nit);
	  k+=2 ;
	}
      if (!strncmp(argv[j],"-ogsf",6)){
	get_char_arg(argv, j, i, numarg2, opt->gsfile) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-wz",3)){
	get_num_arg(argv, j, i, numarg2,"%lf", &opt->wZ) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-wn",3)){
	get_num_arg(argv, j, i, numarg2,"%lf", &opt->wN) ;
	k+=2 ;}
      if (!strncmp(argv[j],"-we",3)){
	get_num_arg(argv, j, i, numarg2,"%lf", &opt->wE) ;
	k+=2 ;}
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
	opt->ref_val = 0. ;
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
get_param1(argc, argv, M, opt, eq, flag) 
     int   argc, *M, *flag ;
     char  **argv   ;
     str_quake_params *eq      ;
     structopt *opt ; 
{
  int numarg1, numarg2 ;
  int max = 47 ;

  numarg1 = 0              ;
  numarg2 = argc-numarg1-1 ;

  if ( (argc < numarg1+1 ) || (numarg2 > max))
    error_syntax(argv," (nb of arguments) ");
  get_opt(numarg1, numarg2, argv, opt, eq);
  *M    = NM  ;
  *flag = 1   ;
  if (opt->ntr_val > 0.)
    *M = NM-1 ;
  if (opt->ref_val > 0.)
    *flag = 2 ;
}

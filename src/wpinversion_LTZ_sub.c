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

/* To avoid the warning with gcc (Debian 4.3.2-1.1) */
double round(double x);

void 
output_products(structopt *opt, str_quake_params *eq, double s1a, double d1a, 
		double r1a, double s2a, double d2a, double r2a, double **TMa, 
		double *eval3a, double M0a, double M0_12a, double Mwa, 
		double *global_rms, double gap, double Cond, int nsac, 
		sachdr *hd_synt, FILE *o_log) 
{
  int    i,j,k,nb,nb2,size,*nbcmp  ;
  double M0b=0.,Mwb=0.,M0_12b=0.,ma,mb,mc ;
  double s1b,s2b,d1b,d2b,r1b,r2b,tmp ;
  double diplow,**TMb,*eval3b    ;
  char   *buf,*buf2,**sta,**cmp  ;
  FILE   *ps,*o_bitmap           ;

  char   date_stmp[64] ;
  time_t now           ; 
  
  nbcmp  = int_alloc(nsac) ;
  buf    = char_alloc(9) ;
  buf2   = char_alloc(32) ;
  sta    = char_calloc2(nsac, 9) ;
  cmp    = char_calloc2(nsac,30) ;    
  eval3b = NULL;
  TMb    = NULL;

  /* Reference solution */
  if (opt->ref_flag)
    {
      eval3b = double_alloc(3)  ;
      TMb = double_alloc2(3,3) ;
      get_planes(eq->vm[1], TMb, eval3b, &s1b,&d1b,&r1b, &s2b,&d2b,&r2b) ;
      M0b = ((fabs(eval3b[0]) + fabs(eval3b[2])) * (double)POW) / 2. ; 
      Mwb = (log10(M0b) - 16.1) / 1.5 ;
      residual_moment(eq->vm, &ma, &mb, &mc) ;
    }

  /* PS FILE */
  if (opt->ps)
    {
      ps = openfile_wt(opt->psfile) ;
      /* header */
      fprintf(ps,"%%!PS\n200 610 translate\n100 100 scale\n");
      fprintf(ps,"0 setlinewidth\n") ;      
      /* Focal mechanism display - WPhase Solution */
      tmp=0.0;
      if(opt->ref_flag)
	tmp=0.5;
      fprintf(ps,"%5.2f %5.2f translate\n",1.1-tmp,0.12+tmp/10.0) ;
      prad_pat(TMa, ps) ;
      pnod_pat(&s1a,&d1a,ps) ;
      pnod_pat(&s2a,&d2a,ps) ;
      fprintf(ps,"/Times-Roman findfont .12 scalefont setfont\n") ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -0.50, -1.13) ;
      fprintf(ps,"(WCMT, Mw= %5.2f ) show\n", Mwa) ;
      if(opt->dc_flag)
	{
	  fprintf(ps, "%15.6f %15.6f moveto\n", -0.59, -1.25) ;
	  fprintf(ps,"(Double-Couple solution) show\n") ;
	}
      fprintf(ps,"%5.2f %5.2f translate\n",tmp-1.5,-tmp/10.0-0.12) ;
      /* title */
      fprintf(ps,"/Times-Roman findfont .2 scalefont setfont\n") ;
      fprintf(ps,"%15.6f %15.6f moveto\n", 1.5, +1.3) ;
      fprintf(ps,"(%s) dup stringwidth pop 2 div neg 0 rmoveto show\n", eq->evnm) ;      

      /* tensor elem., moment, planes, eigvalues */
      fprintf(ps,"/Times-Roman findfont .1 scalefont setfont\n") ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.3) ;
      fprintf(ps,"(Moment Tensor:") ;
      for(i=0 ; i<NM ; i++)
	fprintf(ps," %11.5f",eq->vm[0][i]) ;
      fprintf(ps,") show\n") ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.45) ;
      fprintf(ps,"(Scalar moment: %9.2e dyn cm    (Mo_12 = %9.2e dyn cm)) show\n",M0a,M0_12a) ;	      
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.60) ;
       if (opt->dc_flag)
	 fprintf(ps,"(N");
       else
	 fprintf(ps,"(Best N");
       fprintf(ps,"odal planes (strike/dip/rake):) show\n");
       fprintf(ps,"%15.6f %15.6f moveto\n", 0.5, -1.60) ;
       if (opt->ref_flag) 
	 {	   
	   fprintf(ps,"(WCMT: %7.1f/%5.1f/%7.1f   %7.1f/%5.1f/%7.1f) show\n",
		   s1a,d1a,r1a,s2a,d2a,r2a) ;	       
	   fprintf(ps,"%15.6f %15.6f moveto\n",0.5, -1.75) ;
	   fprintf(ps,"(RCMT : %7.1f/%5.1f/%7.1f   %7.1f/%5.1f/%7.1f) show\n",
		   s1b,d1b,r1b,s2b,d2b,r2b) ;
	 }
       else
	 {
	   fprintf(ps,"(NP1: %7.1f/%5.1f/%7.1f) show\n",s1a,d1a,r1a);
	   fprintf(ps,"%15.6f %15.6f moveto\n", 0.5, -1.75) ;
	   fprintf(ps,"(NP2: %7.1f/%5.1f/%7.1f) show\n",s2a,d2a,r2a);
	 }
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -1.9) ;
      fprintf(ps,"(Eigenvalues:%13.5f %13.5f %13.5f     (Mw = %4.2f)) show\n",
	      eval3a[0],eval3a[1],eval3a[2],Mwa) ;
      
      /* fit quality */
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.1) ;
      fprintf(ps,"(Fit Quality:) show\n") ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -0.4, -2.1) ;
      if (opt->ref_flag) 
	{
	  if (opt->dc_flag)
	    fprintf(ps,"(WCMT - RMS: %9.5f mm (%6.3f),  Gap: %5.1f\\312) show\n",
		    1000.*global_rms[0], global_rms[0]/global_rms[1],gap);
	  else
	    fprintf(ps,"(WCMT - RMS: %9.5f mm (%6.3f),  Gap: %5.1f\\312,  C#%9.0f) show\n",
		    1000.*global_rms[0], global_rms[0]/global_rms[1],gap,Cond);
	  
	  fprintf(ps,"%15.6f %15.6f moveto\n", -0.4, -2.25) ;
	  fprintf(ps,"(RCMT - RMS : %9.5f mm (%6.3f)) show\n", 1000.*global_rms[2], global_rms[2]/global_rms[3]);  
	}
      else
	{
	  fprintf(ps,"(RMS: %9.5f mm (%6.3f)) show\n",1000.*global_rms[0], global_rms[0]/global_rms[1]);
	  fprintf(ps,"%15.6f %15.6f moveto\n", -0.4, -2.25) ;
	  if (opt->dc_flag)
	    fprintf(ps,"(GAP: %5.1f\\312) show\n",gap);
	  else
	    fprintf(ps,"(GAP: %5.1f\\312,  C#%9.0f) show\n",gap,Cond);
	  
	}

      /* used stations */
      k = 0 ;
      for(i=0; i<nsac; i++) /* Set list of channels per stations */
	{
	  nb  = nbchar(hd_synt[i].kstnm)  ;
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
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.45)             ;
      fprintf(ps,"(Used stations (%d, %d channels) : ) show\n",k,nsac ) ;
      fprintf(ps,"/Courier      findfont .07 scalefont setfont")  ;
      
      j  = 0;
      nb = 1;
      fprintf(ps,"%15.6f %15.6f moveto\n(", -1., -2.55) ;
      for(i=0; i<k; i++)
	{
	  j += 1 + nbcmp[i] ;
	  fprintf(ps, "  %s(%s)", sta[i], cmp[i]) ;
	  if (j >= 21 && i<k-1) {
	    fprintf(ps,") show\n") ;
	    fprintf(ps,"%15.6f %15.6f moveto\n(", -1., -2.55-((double)nb)/8.) ;
	    j = 0 ;
	    nb++; }
	}
      fprintf(ps,") show\n")   ;
      
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.55-((double)nb)/8.) ;
      fprintf(ps, "(WPWIN: %-8.2f %-8.2f %-8.2f %-8.2f ) show\n"
	      , eq->wp_win4[0], eq->wp_win4[1], eq->wp_win4[2], eq->wp_win4[3]) ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.65-((double)nb)/8.) ;
      fprintf(ps, "(Dmin : %-8.2f Dmax :%-8.2f) show\n", opt->dmin, opt->dmax) ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -2.75-((double)nb)/8.) ;
      fprintf(ps, "(wL   : %-8.2f wT   :%-8.2f wZ   :%-8.2f) show\n", opt->wL, opt->wT, opt->wZ);  
      
      /* filter parameters */
      fprintf(ps,"/Courier-Bold findfont .1 scalefont setfont") ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.2) ;
      fprintf(ps,"(Filter parameters: ) show\n") ;
      fprintf(ps,"/Courier     findfont .1 scalefont setfont")  ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.3)   ;
      fprintf(ps,"(filt_order: %-d) show\n", eq->filtorder) ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.425) ;
      fprintf(ps,"(filt_cf1  : %-7.5f) show\n", eq->flow) ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.550) ;
      fprintf(ps,"(filt_cf2  : %-7.5f) show\n", eq->fhigh) ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.675) ;
      fprintf(ps,"(filt_pass : %-d) show\n", eq->filtnpass) ;
      
      /* pde  */
      fprintf(ps,"/Courier-Bold findfont .07 scalefont setfont\n");
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -4.9) ;
      fprintf(ps,"(PDE and Centroid: ) show\n")       ;
      fprintf(ps,"/Courier      findfont .07 scalefont setfont\n");
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -5.0)    ;
      nb = nb_blank(eq->pdeline) ;
      eq->pdeline[60]='\0' ;
      fprintf(ps,"(%s) show\n", &eq->pdeline[nb])             ;
      eq->pdeline[60]=' ' ;
      fprintf(ps,"%15.6f %15.6f moveto\n", -1., -5.1)         ;
      fprintf(ps, "(Event id     : %s) show\n", eq->evid)     ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.2)        ;
      fprintf(ps, "(Time shift   : %-6.1f s)  show\n", eq->ts);
      fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.3)        ;
      fprintf(ps, "(Half duration: %-6.1f s)  show\n", eq->hd);
      fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.4)        ;
      
      fprintf(ps, "(Latitude     : %-8.3f) show\n", eq->evla) ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.5)        ;
      fprintf(ps, "(Longitude    : %-8.3f) show\n", eq->evlo) ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.6)        ;
      fprintf(ps, "(Depth        : %-8.3f) show\n", eq->evdp) ;
      fprintf(ps, "%15.6f %15.6f moveto\n", -1., -5.7)        ;  

      /* Time stamp */
      time(&now) ;
      fprintf(ps,"/Times-Roman findfont .1 scalefont setfont\n") ;
      strftime(date_stmp,64,"Processed Date : %a %b %02d %02H:%02M:%02S %Y UT",gmtime(&now));
      date_stmp[44] = '\0';
      fprintf(ps,"%15.6f %15.6f moveto\n",  2.1, -5.8) ;
      fprintf(ps,"(%s) show\n", date_stmp) ;
      
      /* Reference solution */
      if (opt->ref_flag) 
	{
	  fprintf(ps,"2.4 -0.45 translate\n") ;
	  fprintf(ps,"0.4  0.4 scale\n") ;
	  fprintf(ps,"/Times-Roman findfont .2 scalefont setfont\n") ;
	  fprintf(ps, "%15.6f %15.6f moveto\n", -1.1, -1.4) ;
	  fprintf(ps,"(RCMT, Mw= %5.2f ) show\n", Mwb)      ;
	  fprintf(ps, "%15.6f %15.6f moveto\n", -1.1, -1.6) ;
	  fprintf(ps,"(ratio = %5.2f ;  epsilon = %6.3f) show\n",M0b/M0a,mc) ;
	  fprintf(ps,".5 .5 .9 setrgbcolor\n") ;
	  prad_pat(TMb, ps)        ;
	  pnod_pat(&s1b, &d1b, ps) ;
	  pnod_pat(&s2b, &d2b, ps) ;
	}
      
      
      fprintf(ps,"showpage\n") ;
      fclose(ps) ;
    }

  /* STDOUT AND LOG */
  charplot(eq->vm[0],s1a,d1a,s2a,d2a, '-', '#', ' ', '\0','\0','\0', RX, RY, stdout)  ;
  printf("CMT: ") ;
  for(i =0 ; i<NM ; i++)
    printf("%15.4e",eq->vm[0][i]);
  if (opt->dc_flag)
    printf("\nN");
  else
    printf("\nBest n");
  printf("odal planes: %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n",
	 s1a,d1a,r1a,s2a,d2a,r2a);
  printf("Eigenvalues: %-12.5f %-12.5f %-12.5f\n", eval3a[0], eval3a[1], eval3a[2]) ;
  if (opt->dc_flag)
    printf("WCMT: RMS = %-9.5f mm (%-6.3f) Gap: %-6.1f deg\n",
	   1000.*global_rms[0], global_rms[0]/global_rms[1],gap);
  else
    printf("WCMT: RMS = %-9.5f mm (%-6.3f) Gap: %-6.1f deg, C#: %-10.0f\n",
	   1000.*global_rms[0], global_rms[0]/global_rms[1],gap,Cond);
  if (opt->ref_flag) 
    { 
      printf("RCMT: RMS = %-9.5f mm (%-6.3f)\n", 1000.*global_rms[2], global_rms[2]/global_rms[3])  ;      
      diplow = d2b ;
      if (d1b < d2b)
	diplow = d1b ;
      M0_12b = M0b * sin(2.*diplow*(double)DEG2RAD) / sin(24.*(double)DEG2RAD) ; 
    }
  printf("Wmag: %-5.2f ; Wmom %-15.4e ; Wmom_12 %-15.4e\n",Mwa,M0a,M0_12a) ;
  
  fprintf(o_log,"n_used_rec:         %-4d\n",nsac) ;
  fprintf(o_log,"Gap:                %-7.1f\n",gap) ;
  fprintf(o_log,"WCMT: ") ;
  for(i =0 ; i<NM ; i++)
    fprintf(o_log,"%15.4e",eq->vm[0][i]);
  fprintf(o_log,"\nCond_number:        %-10.0f\n",Cond) ; 
  fprintf(o_log,"W_bestnodal planes: %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n",
	 s1a,d1a,r1a,s2a,d2a,r2a) ;
  fprintf(o_log,"W_eigenvalues:      %-12.5f %-12.5f %-12.5f\n", eval3a[0], eval3a[1], eval3a[2]) ;
  fprintf(o_log,"W_cmt_err:          %-12.8f %-12.8f\n", 1000.*global_rms[0], global_rms[0]/global_rms[1])    ;
  fprintf(o_log,"Wmag: %-5.2f ; Wmom %-15.4e ; Wmom_12 %-15.4e\n",Mwa,M0a,M0_12a) ;

  if (opt->ref_flag) 
    {
      printf("Rmag: %-5.2f ; Rmom %-15.4e ; Rmom_12 %-15.4e\n",Mwb,M0b,M0_12b) ;
      printf("ratio = %5.2f ;  epsilon = %6.3f\n",M0b/M0a,mc) ;
      fprintf(o_log,"R_bestnodal planes: %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n",s1b,d1b,r1b, s2b,d2b,r2b) ;
      fprintf(o_log,"R_eigenvalues:      %-12.5f %-12.5f %-12.5f\n", eval3b[0], eval3b[1], eval3b[2])           ;
      fprintf(o_log,"R_cmt_err:          %-12.8f %-12.8f\n", 1000.*global_rms[2], global_rms[2]/global_rms[3])  ;
      fprintf(o_log,"Rmag: %-5.2f ; Rmom %-15.4e ; Rmom_12 %-15.4e\n",Mwb,M0b,M0_12b)                           ;
      fprintf(o_log,"ratio = %12.8f ;  epsilon = %12.8f\n",M0b/M0a,mc) ;
    }

  /* BITMAP IMAGE */
  if (strlen(opt->wpbmfile) != 0)
	{
	  size = 401 ;
	  o_bitmap = openfile_wt(opt->wpbmfile)  ;
	  fprintf(o_bitmap, "P2\n#WPinversion focal mechanism\n%d %d\n9\n", size, size)                             ;
	  charplot(eq->vm[0],s1a,d1a,s2a,d2a,'9', '3', '9', '0',' ', '0', (size-1)/2,(size-1)/2, o_bitmap)  ;
	  fclose(o_bitmap) ;
	}

  if (strlen(opt->refbmfile) != 0 && opt->ref_flag) 
	{
	  o_bitmap = openfile_wt(opt->refbmfile) ;
	  fprintf(o_bitmap, "P2\n#WPinversion focal mechanism\n%d %d\n9\n", size, size)                         ;
	  charplot(eq->vm[1], s1b, d1b, s2b, d2b,'9', '3', '9', '0',' ', '0', (size-1)/2,(size-1)/2, o_bitmap)  ;
	  fclose(o_bitmap)                       ; 
	}

  /* Memory Freeing */
  if (opt->ref_flag)
    {
      for(i=0 ; i<3 ; i++) 
	free((void*)TMb[i]) ;
      free((void**)TMb)     ;
      free((void*)eval3b)   ;
    }
  free((void*)nbcmp)    ;
  free((void*)buf)      ;
  free((void*)buf2)     ;
  for(i=0 ; i<nsac ; i++)
    {
      free((void*)sta[i]) ;
      free((void*)cmp[i]) ; 
    }
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
  int    i, j, k, l, N = 20, M ;
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
pnod_pat(double *s, double *d, FILE *ps)
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
get_gap(sachdr *hd_synt, int ns, double *gap)
{
  double *tmp, az, dist ;
  int i ;
  
  tmp = double_alloc(ns) ;
  for( i=0 ; i<ns ; i++)
    {
      az = hd_synt[i].az ;
      if(az < 0.) 
	tmp[i] = az + 360. ;
      else
	tmp[i] = az ;
    }
  
  sort(tmp,ns) ;
  
  *gap = tmp[0] + 360. - tmp[ns-1] ;
  for( i=0 ; i<(ns)-1 ; i++)
    {
      dist = tmp[i+1]-tmp[i] ;
      if(dist > *gap) *gap = dist ;
    }
  free((void*)tmp) ;
}

void 
w_o_saclst(int ns, char **sacfiles, sachdr *hd_synt, double **rms, double *data_norm, 
	   structopt *opt) 
{
  int  i, n0;
  FILE *o_sac ;

  o_sac = openfile_wt(opt->o_saclst);
  n0 = 0;
  for(i=0; i<ns; i++)
    {
      fprintf(o_sac,"%-65s %8.2f %8.2f %6d %6d %6d %6d %14.8f %14.8f %14.8f %14.8f %14.8f %10.3f\n",
	      sacfiles[i], hd_synt[i].az, hd_synt[i].gcarc, n0, n0 + hd_synt[i].npts, 
	      (int)hd_synt[i].user[0], (int)hd_synt[i].user[1], rms[i][0], rms[i][0]/rms[i][1],
	      data_norm[i]/rms[i][1], opt->p2p[i], opt->avg[i], opt->wgt[i]);
      n0 += hd_synt[i].npts ;
    }
  fclose(o_sac);
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
calc_stat(int npts, double *data, double *mi2ma, double *ave) 
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
/* calc_rms(ns, hd_synt, data, dcalc, rms, global_rms)       */
/*************************************************************/
/* calculate rms amplitudes, rms deviations                  */
/* Input: ns : nb of channels                                */
/*        hd_synt : array of sac headers                     */
/*        data    : array of data                            */
/* Output: dcalc  : array of predicted data                  */
/*         rms    : array of rms amplitudes and deviations   */
/*         global_rms : overall rms amp. and rms deviations  */
void 
calc_rms(int ns, sachdr *hd_synt, double **data, double ***dcalc, double **rms, double *global_rms, structopt *opt) 
{
  int    i, j, n0, f2;

  f2  = 2*(opt->ref_flag+1);
  n0 = 0 ;
  for(i=0;i<ns;i++)
    {
      calc_rms_sub(hd_synt[i].npts, data[i], dcalc[i], rms[i],opt->ref_flag);
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
calc_rms_sub(int npts, double *data, double **dcalc, double *rms, int flag)
{
  int i;

  if (flag)
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
  else
    {
      rms[0] = 0. ;
      rms[1] = 0. ;
      for(i=0 ; i<npts ; i++){
	rms[0] += (data[i]-(*dcalc)[i]) * (data[i]-(*dcalc)[i]) ;
	rms[1] +=  (*dcalc)[i]*(*dcalc)[i];}
    }
}

void
calc_vec_norm(double *v, int npts, double *norm)
{
  int i;
  *norm = 0;
  for(i=0;i<npts;i++)
    *norm += v[i]*v[i] ;
  *norm = sqrt(*norm/(double) npts) ;  
}

void 
calc_data_norm(double **data, sachdr *hd_synth, int nsac, double *data_norm)
{
  int s;
  for(s=0;s<nsac;s++)
    calc_vec_norm(data[s],hd_synth[s].npts,data_norm+s);
}

/*************************************************************/
/*     calc_data(ns, hd_synt, G, vm, data, d, opt)           */
/*************************************************************/
/* Compute predicted data                                    */
/* Input: nsac : nb of channels                              */
/*        hd_synt : array of sac headers                     */
/*        hd_synt : array of sac headers                     */
/*        G       : array of Green's functions               */
/*        vm      : pointer on moment tensor elements        */
/* Output: data   : array data samples                       */
/*         d      : array of predicted data samples          */
/*         opt    : structopt containing output filenames    */
void 
calc_data(int nsac, sachdr *hd_synt, double ***G, double **vm, 
	  double **data, double ***d, structopt *opt, FILE *o_log)
{
  int  i, j, k, s, N = 0    ;
  char *fZ, *fL, *fT, *file ;
  FILE *o_dat, *o_Z, *o_L, *o_T, *ocmp=NULL;

  /* Computing data */
  for(s=0;s<nsac;s++)
    {
      if (hd_synt[s].kcmpnm[2] != 'Z' && hd_synt[s].kcmpnm[2] != 'L' && hd_synt[s].kcmpnm[2] != 'T')
	{
	  fprintf(stderr,"ERROR calc_data : invalid component\n") ;
	  exit(1) ;       
	}
      N = hd_synt[s].npts ;
      d[s] = double_calloc2(opt->ref_flag+1, N) ;
      for(i=0 ; i<N ; i++)
	for(j=0; j<opt->ref_flag+1; j++)  
	  for(k=0; k<NM; k++)
	    d[s][j][i] += G[s][k][i] * vm[j][k] ;
    }
  /* Write output files */
  if (o_log != NULL)
    {
      /* Allocating memory */
      fZ = char_alloc(FSIZE) ;
      fL = char_alloc(FSIZE) ;
      fT = char_alloc(FSIZE) ;
      /* Set filenames     */
      strcpy(fZ,opt->p_data) ;
      strcpy(fL,opt->p_data) ;
      strcpy(fT,opt->p_data) ;
      strcat(fZ,"_LHZ") ;
      strcat(fL,"_LHL") ;
      strcat(fT,"_LHT") ;
      
      /* Opening files */
      o_dat = openfile_wt(opt->p_data) ;
      o_Z = openfile_wt(fZ) ;
      o_L = openfile_wt(fL) ;
      o_T = openfile_wt(fT) ;
      for(s=0;s<nsac;s++)
	{
	  if (hd_synt[s].kcmpnm[2] == 'Z')
	    ocmp = o_Z ;
	  else if (hd_synt[s].kcmpnm[2] == 'L')
	    ocmp = o_L ;
	  else if (hd_synt[s].kcmpnm[2] == 'T')
	    ocmp = o_T ;
	  N = hd_synt[s].npts ;
	  for(i=0 ; i<N ; i++)
	    {
	      fprintf(o_dat,"%15.6e ",data[s][i]) ;
	      fprintf( ocmp,"%15.6e ",data[s][i]) ; 	     
	      for(j=0; j<opt->ref_flag+1; j++)  
		{
		  fprintf(o_dat,"%15.6e ",d[s][j][i]); 
		  fprintf( ocmp,"%15.6e ",d[s][j][i]); 		  
		}
	      fprintf(o_dat,"\n") ;
	      fprintf(ocmp,"\n") ;
	    }
	  file = get_gf_filename(opt->osacdir, hd_synt[s].kstnm, 
				 hd_synt[s].knetwk, hd_synt[s].kcmpnm[2], 
				 "synth.sac") ;
	  wsac(file, &hd_synt[s], d[s][0]) ;
	  /* Memory Freeing */
	  free((void*)file) ;
	}
      fclose(o_dat) ;
      fclose(o_Z)   ;
      fclose(o_L)   ;
      fclose(o_T)   ;
      free((void*)fZ)   ;
      free((void*)fL)   ;
      free((void*)fT)   ;
    }
}

void 
write_cmtf(char *filename, str_quake_params *eq, double *vm)
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
comp_GtG(int M,int nsac, sachdr *hd_synt, double ***G, double **GtG, structopt *opt)
{
  int i, j, k, s, N;

  if (GtG[0][0] != 0.)
    {
      fprintf(stderr,"WARNING GtG: memory must be set to zero (using calloc)\n");
      for( i=0;i<M;i++)
	for( j=0;j<M;j++)
	  GtG[i][j] = 0. ;
    }
  if (M==NM-1)    /* Constrain null trace */
    {
      for(s=0;s<nsac;s++)
	{
	  N = hd_synt[s].npts ;
	  for(i=0;i<2;i++)
	    {
	      for(j=i;j<2;j++ )
		for(k=0;k<N;k++)
		  GtG[i][j] += (G[s][i][k]-G[s][2][k])*(G[s][j][k]-G[s][2][k])*opt->wgt[s]/opt->wgtnorm ;
	      for(j=2;j<M;j++)
		for(k=0;k<N;k++)
		  GtG[i][j] += (G[s][i][k]-G[s][2][k])*G[s][j+1][k]*opt->wgt[s]/opt->wgtnorm ;
	    }
	  for(i=2;i<M;i++)
	    for(j=i;j<M;j++)
	      for(k=0;k<N;k++)
		GtG[i][j] += G[s][i+1][k]*G[s][j+1][k]*opt->wgt[s]/opt->wgtnorm ;
	}
    }
  else if (M==NM) /* No constraints      */
    for(s=0;s<nsac;s++)
      for(i=0;i<M;i++)
	for(j=i;j<M;j++)
	  for( k=0 ; k<hd_synt[s].npts ; k++ )
	    GtG[i][j] += G[s][i][k]*G[s][j][k]*opt->wgt[s]/opt->wgtnorm ;
  for(i=1;i<M;i++) /* Finish the work */
    for(j=0;j<i;j++)
      GtG[i][j] = GtG[j][i] ;
}


void 
inversion(int M, int nsac, sachdr *hd_synt, double ***G, double **d, 
	  double *vma, double *Cond, structopt *opt, FILE *o_log)
{
  int    i, j ,k, l, s, nk, nrot, N  ;
  double **GtG, *eigvals, **eigvects, **cov ;
  FILE   *o_cov ;

  /* Allocating memory */
  GtG      = double_calloc2 (M,M) ;
  eigvals  = double_alloc(M)      ;
  eigvects = double_alloc2(M,M)   ;
  cov      = double_calloc2(M,M)  ;
  /* Azimuth ponderation */
  if (opt->azp > 0. && opt->dts_step <= 0. && opt->xy_dx <= 0.)
    azpond(hd_synt,nsac,opt->wgt)    ;
  /* Weight normalization */
  opt->wgtnorm = 0.;
  for(i=0;i<nsac;i++)
    opt->wgtnorm += opt->wgt[i] ;
  opt->wgtnorm /= (double)nsac ;  
  /* Constrains null trace (if M=5) and computes GtG */
  comp_GtG(M,nsac,hd_synt,G,GtG,opt) ;
  /* Computes eigenvalues and eigenvectors */
  jacobi(GtG,M,M,eigvals,eigvects,&nrot);
  eigsrt(eigvals,eigvects,M) ;
  /* Conditioning */
  nk = M;
  if (opt->cth_val <= 0.) /* Remove the eigenvalues smaller than max(eigval)/10000 if any */
    {
      *Cond = eigvals[0] / eigvals[nk-1] ;
      if (*Cond > 1.e4)
	for(i=1;i<nk;i++) 
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
	  for(i=1;i<=nk;i++)
	    eigvals[i] = eigvals[i] + eigvals[0]*(opt->df_val) ;
	}
      *Cond = eigvals[0] / eigvals[nk-1] ;
    }  
  /* Display inversion details */
  if (o_log != NULL)
    {
      o_cov    = openfile_wt("Gd") ;
      for(i=0;i<nsac;i++)
	for(k=0;k<hd_synt[i].npts;k++)
	  {
	    for(j=0;j<NM;j++)
	      fprintf(o_cov,"%16.8e ",G[i][j][k]);
	    fprintf(o_cov,"%16.8e\n",d[i][k]);
	  }
      fclose(o_cov) ;

      if (opt->dts_step <= 0. && opt->xy_dx <= 0.)
	{
	  printf("##############\n%d significant eigenvalues: \n", nk) ;
	  for(i=0;i<nk;i++) 
	    printf("  %13.5e",eigvals[i]);
	  printf("\n");
	}
      fprintf(o_log,"Inversion: %d significant eigenvalues: \n", nk) ;  
      for(i=0;i<nk;i++) 
	fprintf(o_log,"%e\t",eigvals[i]) ;   
      fprintf(o_log,"\n") ;
      if(M-nk>0)
	{
	  printf(       "%d removed: \n",M-nk) ;
	  fprintf(o_log,"%d removed: \n",M-nk) ;
	  for (i=nk;i<M;i++) {
	    printf("  %13.5e",eigvals[i])    ;
	    fprintf(o_log,"%e\t",eigvals[i]) ;   }
	  printf("\n")        ;
	  fprintf(o_log,"\n") ;
	  fflush(o_log)       ;
	}
    }
  /* Posterior covariance matrix */
  for(i=0;i<M;i++)
    for(k=0;k<M;k++)
      for(j=0;j<nk;j++)
	cov[i][k] += eigvects[i][j]*eigvects[k][j]/eigvals[j] ;
  if (o_log != NULL)
    {
      o_cov    = openfile_wt(opt->o_covf) ;
      for(i=0;i<M;i++)
	{
	  for(k=0;k<M;k++)
	    fprintf( o_cov,"%16.2f ", cov[i][k]) ;
	  fprintf( o_cov,"\n") ;
	}
      fclose(o_cov) ;
    }
    
  /* Least-squares inversion */
  for(i=0;i<NM;i++)  /* Initialize MT components */
    vma[i] = 0. ;
  if (M==NM-1)    /* Null trace     */
    {
      for(s=0;s<nsac;s++)
	{
	  N = hd_synt[s].npts ;
	  for(i=0;i<M;i++)
	    {
	      for(k=0;k<2;k++)
		for(l=0;l<N;l++)
		  vma[i] += cov[i][k]*(G[s][k][l]-G[s][2][l])*d[s][l]*opt->wgt[s]/opt->wgtnorm ;
	      for(k=2;k<M;k++)
		for(l=0;l<N;l++)
		  vma[i] += cov[i][k]*G[s][k+1][l]*d[s][l]*opt->wgt[s]/opt->wgtnorm ;
	    }
	}
      vma[5] = vma[4] ;
      vma[4] = vma[3] ;
      vma[3] = vma[2] ;
      vma[2] = -vma[0] -vma[1] ;
    }
  else if (M==NM) /* No constraints */
    for(s=0;s<nsac;s++)
      for(i=0;i<M;i++)
	for(k=0;k<M;k++)
	  for(l=0 ; l<hd_synt[s].npts ; l++)
	    vma[i] += cov[i][k]*G[s][k][l]*d[s][l]*opt->wgt[s]/opt->wgtnorm ;
  else 
    {
      fprintf(stderr,"ERROR : bad nb of parameters (M=%d)\n",M);
      exit(1);  
    }
  /* Memory Freeing */
  free((void*)eigvals) ;
  for(i=0;i<M;i++) 
    { 
      free((void*)GtG[i])      ; 
      free((void*)eigvects[i]) ; 
      free((void*)cov[i])      ;
    }
  free((void**)GtG)      ; 
  free((void**)eigvects) ; 
  free((void**)cov)      ; 
}
  

void
sdr2mt(double *vm,double Mo,double strike,double dip,double rake)
{  
  double ssi,sco,ssi2,sco2,dsi,dco,dsi2,dco2,rsi,rco;
  ssi  = sin(strike*DEG2RAD);
  sco  = cos(strike*DEG2RAD);
  ssi2 = sin(2*strike*DEG2RAD);
  sco2 = cos(2*strike*DEG2RAD);
  dsi  = sin(dip*DEG2RAD);
  dco  = cos(dip*DEG2RAD);
  dsi2 = sin(2*dip*DEG2RAD);
  dco2 = cos(2*dip*DEG2RAD);
  rsi  = sin(rake*DEG2RAD);
  rco  = cos(rake*DEG2RAD);
  vm[0] =  Mo * dsi2*rsi                             ;
  vm[1] = -Mo * (dsi*rco*ssi2 + dsi2*rsi*ssi*ssi)    ;
  vm[2] =  Mo * (dsi*rco*ssi2 - dsi2*rsi*sco*sco)    ;
  vm[3] = -Mo * (dco*rco*sco  + dco2*rsi*ssi)        ;
  vm[4] =  Mo * (dco*rco*ssi  - dco2*rsi*sco)        ;
  vm[5] = -Mo * (dsi*rco*sco2 + 0.5 * dsi2*rsi*ssi2) ;  
}


/* Warning: This routine has not been fully tested */
void 
inversion_dc(int nsac, sachdr *hd_synt, double ***G, double **d,
	     double *sdrM0,double *rmsopt, structopt *opt, FILE *o_log)
{
  int    i, j, k, l, npts, ref_flag ;
  int ip,ib[4],mp,mv,idvt,icon,iquit,iprnt;
  float  *data,*xlsq,b[4];
  double ***dcalc, **rms, *global_rms ;
  double **vm;
  
  ip = opt->ip;
  for(i=0;i<ip;i++)
    ib[i]=opt->ib[i];
  /* Allocate memory */
  vm    = double_alloc2p(2) ;
  vm[0] = double_calloc(NM) ;
  dcalc = double_alloc3p(nsac)   ;      
  rms   = double_calloc2(nsac,2) ;
  global_rms = double_calloc(2)  ;
  npts = 0; /* # of samples */
  for(i=0;i<nsac;i++)
   npts += hd_synt[i].npts ;
  data = float_alloc(npts)    ;
  xlsq = float_alloc(NM*npts) ;
  /* Reference flag */
  ref_flag = opt->ref_flag ;
  opt->ref_flag = 0        ;      
  /* Set the lsq parameters */
  mp = 4;    /* no. of param            */
  mv = 1;    /* no. of indep. variables */
  idvt  = 1; /* Estimated derivatives will be used           */  
  icon  = 1; /* Omit nonlinear confidence limits calculation */
  iquit = 0; /* No force off         */
  iprnt = 0; /* abbreviated printout */
  /* A priori RMS value */
/*   if (o_log != NULL) */
/*     { */
      sdr2mt(*vm,sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2])   ;
      calc_data(nsac,hd_synt,G,vm,d,dcalc,opt,NULL)     ;
      calc_rms(nsac,hd_synt,d,dcalc,rms,global_rms,opt) ;       
      for(i=0;i<2;i++)   
	{
	  rmsopt[i]=global_rms[i];
	  global_rms[i] = 0.;
	  for(j=0;j<nsac;j++)
	      rms[j][i] = 0.;
	}
      for(i=0;i<nsac;i++)
	{
	  free((void*) dcalc[i][0]) ;
	  free((void**)dcalc[i])    ;
	}
/*     } */
  /* Weight normalization */
  opt->wgtnorm = 0.;
  for(i=0;i<nsac;i++)
    opt->wgtnorm += opt->wgt[i] ;
  opt->wgtnorm /= (double)nsac ;  
  /* Populate data and wgt array */
  k = 0; 
  for(i=0;i<nsac;i++)
    for(j=0;j<hd_synt[i].npts;j++)
      data[k++] = (float)(d[i][j]*sqrt(opt->wgt[i]/opt->wgtnorm));
  /* Populate xlsq array  */
  l = 0;
  for(i=0;i<NM;i++)
    for(j=0;j<nsac;j++)
      for(k=0;k<hd_synt[j].npts;k++)
	xlsq[l++] = (float)(G[j][i][k]*sqrt(opt->wgt[j]/opt->wgtnorm));
  /* Populate model array */  
  for(i=0;i<4;i++)
    b[i] = (float)sdrM0[i]; 
  /* LSQENP         */
  if (o_log != NULL)
    {
      printf("Double-Couple inversion:\n         Strike  Dip   Rake  M0(dyn.cm) RMS(mm)\n");
      printf("APRIORI: %6.2f %5.2f %6.2f %9.2e %8.5f\n",
	     sdrM0[0],sdrM0[1],sdrM0[2],sdrM0[3]*POW,rmsopt[0]*1000) ;
      fprintf(o_log,"Double-Couple inversion:\n         Strike  Dip   Rake  M0(dyn.cm) RMS(mm)\n");
      fprintf(o_log,"APRIORI: %6.2f %5.2f %6.2f %9.2e %8.5f\n",
	      sdrM0[0],sdrM0[1],sdrM0[2],sdrM0[3]*POW,rmsopt[0]*1000) ;      
    }
  lsqenp_ (&npts,&mp,&mv,data,xlsq,b,&ip,ib,&idvt,&icon,&iquit,&iprnt);
  /* Populate sdrM0 */
  for(i=0;i<4;i++)
    sdrM0[i] = (double)b[i];
  /* Calculate rms  */
  sdr2mt(*vm,sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2])   ;
  calc_data(nsac,hd_synt,G,vm,d,dcalc,opt,NULL)     ;
  calc_rms(nsac,hd_synt,d,dcalc,rms,global_rms,opt) ;   
  for(i=0;i<2;i++)    
    rmsopt[i]=global_rms[i];
  if (o_log != NULL)
    {
      printf("OPTIMUM: %6.2f %5.2f %6.2f %9.2e %8.5f\n",
	     sdrM0[0],sdrM0[1],sdrM0[2],sdrM0[3]*POW,rmsopt[0]*1000) ;
      fprintf(o_log,"OPTIMUM: %6.2f %5.2f %6.2f %9.2e %8.5f\n",
	      sdrM0[0],sdrM0[1],sdrM0[2],sdrM0[3]*POW,rmsopt[0]*1000) ;
    }
  /* Reference flag */
  opt->ref_flag = ref_flag ;
  /* Free memory */
  free((void*)data);
  free((void*)xlsq);
  for(i=0 ; i<nsac ; i++)
    {
      free((void*) rms[i])      ;
      free((void*) dcalc[i][0]) ;
      free((void**)dcalc[i])    ;
    }
  free((void***)dcalc)    ;
  free((void**) rms   )   ;
  free((void*)global_rms) ; 
  free((void*)vm[0]) ;
  free((void**)vm  ) ;
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
residual_moment(double **vm, double *ma, double *mb, double *mc)
{
  int    i;
  double *vmc;
  
  vmc = double_alloc(NM) ;

  mt2sm(vm[0], ma) ;
  mt2sm(vm[1], mb) ;
  for(i=0;i<NM;i++) 
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
mt2sm(double *vm, double *sm) 
{
  *sm = vm[0]*vm[0]+vm[1]*vm[1]+vm[2]*vm[2] ;
  *sm += 2.*(vm[3]*vm[3]+vm[4]*vm[4]+vm[5]*vm[5]) ;
  *sm = sqrt(*sm/2.) ;
}


void
set_mt(double *vm, double **TM)
{
  TM[0][0] =   vm[1] ;  /* Rotation of pi around East-axis       */
  TM[1][1] =   vm[2] ;  /* The new system is (North, East, Down) */
  TM[2][2] =   vm[0] ;  /* (Aki's coordinates)                   */
  TM[0][1] =  -vm[5] ;
  TM[0][2] =   vm[3] ;
  TM[1][2] =  -vm[4] ;
  TM[1][0] = TM[0][1] ;
  TM[2][0] = TM[0][2] ;
  TM[2][1] = TM[1][2] ;
}
  
void 
get_planes(double *vm, double **TM, double *eval3, double *s1,double *d1,
	   double *r1, double *s2,double *d2,double *r2)
{
  int    nrot, i ;
  double si, co  ;
  double **evec3, *vn, *vs ;

  
  /* Memory allocation */
  evec3  = double_alloc2(3,3) ;
  vn     = double_alloc(3)    ;
  vs     = double_alloc(3)    ;
  
  /* Tensor representation */  
  set_mt(vm,TM) ;

  /* Get eigvalues and eigvectors*/
  jacobi(TM,3,3,eval3,evec3,&nrot) ;
  eigsrt(eval3,evec3,3) ;

  /* Check if eigenvectors are upwards */
  if(evec3[2][0] < 0.0) 
    for(i=0;i<3;i++)
      {
	evec3[i][0] = -evec3[i][0] ;
      }

  if(evec3[2][2] < 0.) 
    for(i=0;i<3;i++)
      evec3[i][2] = -evec3[i][2] ;

  /* Cross-product v2 = v1 x v3 */
  evec3[0][1] = evec3[1][0]*evec3[2][2]-evec3[1][2]*evec3[2][0] ;
  evec3[1][1] = evec3[2][0]*evec3[0][2]-evec3[2][2]*evec3[0][0] ;
  evec3[2][1] = evec3[0][0]*evec3[1][2]-evec3[0][2]*evec3[1][0] ;
     
  TM[0][1] = TM[1][0] ; /* useless */
  TM[0][2] = TM[2][0] ;
  TM[1][2] = TM[2][1] ;
  
  /* *** First nodal plane *** */
  for(i=0;i<3;i++)
    {
      vn[i] = (evec3[i][0]+evec3[i][2])/sqrt(2.) ;
      vs[i] = (evec3[i][0]-evec3[i][2])/sqrt(2.) ;
    }
  if (vn[2] > 0.)
    for(i=0;i<3;i++)
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
  for(i=0;i<3;i++)
    {
      vn[i] = (evec3[i][0]-evec3[i][2])/sqrt(2.) ;
      vs[i] = (evec3[i][0]+evec3[i][2])/sqrt(2.) ;
    }

  if (vn[2] > 0.) 
    for(i=0;i<3;i++)
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
  if (hd_data->kcmpnm[2] == 'Z')
    opt->wgt[ns] *= opt->wZ;
  else if (hd_data->kcmpnm[2] == 'L')
    opt->wgt[ns] *= opt->wL;
  else if (hd_data->kcmpnm[2] == 'T')
    opt->wgt[ns] *= opt->wT;
  else
    opt->wgt[ns] = 0. ;
}

/* Warning: This routine has not been fully tested */
void
azpond (sachdr *hd_synt, int ns, double *wgt) 
{
  int    i, j                                ;
  double moy, var, mok, normcov ;
  double onepisq, twopisq, twopidaz, *aztab, *azcov ;
  aztab = double_calloc(ns) ;
  azcov = double_calloc(ns) ;
  for (i=0; i<ns; i++)
    {
      if (hd_synt[i].az < 0.)
	hd_synt[i].az += 360. ;
      aztab[i] = hd_synt[i].az ;
    }
  sort(aztab, ns);
  moy = 0;
  for (i=0; i<ns-1; i++)
    {
      aztab[i] = aztab[i+1]-aztab[i] ;
      moy  += aztab[i] ;
    }
  aztab[ns-1] = 360. + aztab[0] - aztab[ns-1] ;
  moy = (moy + aztab[ns-1])/(double)ns ;
  var = 0.;
  for (i=0; i<ns; i++)
    var += (aztab[i] - moy)*(aztab[i] - moy) ;
  if (var <= VAREPS )
    var = VAREPS ;
  var   = var/ns   ;
  normcov  = 1.e10  ;
  onepisq = 32400  ;
  twopisq = 129600 ;
  for(i=0; i<ns ; i++)
    {
      for(j=0; j<ns; j++)
	{ 
	  mok = (hd_synt[i].az-hd_synt[j].az)*(hd_synt[i].az-hd_synt[j].az) ;
	  if (mok > onepisq)  
	    {
	      mok += twopisq ;
	      twopidaz  = 720.*(hd_synt[i].az-hd_synt[j].az) ;
	      if (twopidaz >= 0)
		mok -= twopidaz ;
	      else
		mok += twopidaz ;
	    }
	  azcov[i] += exp(-mok/var) ;
	}
      if (azcov[i] < normcov)
	normcov = azcov[i] ;
    }  
  for(i=0; i<ns; i++)    
    wgt[i] *= normcov/azcov[i] ;
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
      if (o_log!=NULL)
	fprintf( o_log,"**** Incorrect sampling period, rejected trace : %s\n", datafile) ;
      fprintf(stderr,"**** Incorrect sampling period, rejected trace : %s\n", datafile) ;
      return 1 ;
    }
  /* GF Time Window */
  t0 = (double)hd_GF->o ;
  if (TWPTT)
    t0 += Ptt ;
  else
    t0 += opt->ts;
  n1_GF = (int)((t0 + twp_beg - (double)hd_GF->b - opt->dts_val)  / ((double)hd_GF->delta)) ; /* first GF Sample (corrected) */
  n2_GF = n1_GF + (int)((twp_end - twp_beg) / ((double)hd_GF->delta)) ;                       /* Last GF Sample */
  if ( n2_GF >= hd_GF->npts ) /* GF Rejected */
    {
      if (o_log!=NULL)
	{
	  fprintf(o_log,"Stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", 
		  hd_data->kstnm, hd_data->knetwk, hd_data->kcmpnm, 
		  hd_data->gcarc, hd_data->az, hd_data->user[2], hd_data->user[3]) ;
	  fprintf( o_log,"**** Incomplete GF, rejected : %s\n"  , gf_file)         ;  
	}
      fprintf(stderr,"**** Incomplete GF, rejected : %s ", gf_file)            ; 
      fprintf(stderr,"((n2_GF=%d)>=(npts=%d))\n", n2_GF, hd_GF->npts)          ;
      return 1 ;
    }
  else if ( n1_GF < 0 )       /* Fill negative samples with zeros */
    {
      if (o_log!=NULL)
	{
	  fprintf(o_log,"Stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", 
		  hd_data->kstnm, hd_data->knetwk, hd_data->kcmpnm, 
		  hd_data->gcarc, hd_data->az, hd_data->user[2], hd_data->user[3]) ;
/* 	  fprintf( o_log,"**** Incomplete GF, filling with zeros : %s\n", gf_file) ;   */
	}
/*       fprintf(stderr,"**** Incomplete GF,")                                    ; */
/*       fprintf(stderr,"filling with zeros : %s (n1_GF=%d)\n", gf_file, n1_GF )  ; */
      for(i=n1_GF; i<((n2_GF<0)?n2_GF+1:0); i++)
	G[i-n1_GF] = 0. ;
      g = &G[-n1_GF] ;
      npts += n1_GF  ;
      hd_GF->b += n1_GF*hd_GF->delta ;
      n1_GF = 0 ;
    }
  /* Read GF samples */
  hd_GF->b += hd_GF->delta*(float)n1_GF  ;
  if (n2_GF>=0)
    {
      hd_GF->npts = n2_GF + 1                  ;
      rdatsac(gf_file, hd_GF, buffer, &ierror) ;
      memcpy (g,buffer+n1_GF,npts * sizeof(double));
    }
  //printf("%d %d %e\n",n1_GF,npts,buffer[244]);
  return 0;
}
/* Corrections, Notes, etc.:            */
/* hd_GF->b += (n1_GF-1)*hd_GF->delta ; */


void 
set_matrices (int *nsac, int *nsini,char ***sacfiles,sachdr **hd_synt,double ***data, 
	      double ****G,structopt *opt,str_quake_params *eq,FILE *o_log) 
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
  FILE   *i_sac             ;
  sachdr hd_data, hd_GF     ;


  /* Opening data file list */
  i_sac = openfile_rt(opt->i_saclst,nsini) ;
  *nsac = *nsini ;

  /* Allocating memory */
  dum      = char_alloc(32)    ;
  buf      = char_alloc(LSIZE) ;
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
       *sacfiles == NULL && opt->rms_r == NULL) 
    {
      *data     = double_alloc2p( *nsac )  ;
      *G        = double_alloc3p( *nsac )  ;
      opt->wgt  = double_alloc( *nsac )    ;
      *sacfiles = char_alloc2(*nsac,FSIZE) ;
      opt->rms_in = double_alloc(*nsac)    ;
      opt->rms_r  = double_alloc(*nsac)    ;
      opt->p2p    = double_alloc(*nsac)    ;
      opt->avg    = double_alloc(*nsac)    ;
      hdr_tab(hd_synt, *nsac)              ; 
      flag2 = 1 ;
    }
  /* Set travel times */
  ierror = 1 ; /* error flag */
  trav_time_init(nh,nd,eq->pde_evdp,dv,tv,&ierror) ;  
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
	  check_scan(1, flag, opt->i_saclst, i_sac)  ;
	  opt->wgt[ns] = 1.0;
	}
      else 
	{
	  flag = fscanf (i_sac, "%s %f %f %s %s %s %s %s %lf %lf %lf %lf %lf",
			 datafile,&(*hd_synt)[ns].az,&(*hd_synt)[ns].gcarc,dum,dum,dum,dum,dum,
			 &opt->rms_in[ns],&opt->rms_r[ns],&opt->p2p[ns],&opt->avg[ns],&opt->wgt[ns]) ;
	  strcpy(buf,opt->i_saclst);
	  strcat(buf," (nb of columns may be incorrect)");
	  check_scan(13, flag, buf, i_sac) ;
	} 
      /* Read data header and weights */
      rhdrsac(datafile,  &hd_data, &ierror);
      set_wgt(ns, &hd_data, opt) ;
      if (opt->wgt[ns] <= 0.)
	{
	  fprintf(stderr,"**** null weight, rejected : %s\n", datafile) ;
	  fflush(stderr);
	  continue ;
	}
      /* Data Time Window */
      if ((float)eq->pde_evdp != hd_data.evdp)	
	{
	  fprintf(stderr,"WARNING : depth %f in sac header is different from the pde depth %f in CMTFILE\n",eq->pde_evdp, hd_data.evdp);
	  fprintf(stderr," ... you should carefully re-check gcarc and evdp header variables in file %s\n",datafile); 
	  fprintf(stderr," ... (gcarc from %s is used for the time windowing)\n",datafile); 
	  fflush(stderr);
	}
      gcarc = (double) hd_data.gcarc ;
      if (opt->dmin > gcarc)
	opt->dmin = gcarc ;
      if (opt->dmax < gcarc)
	opt->dmax = gcarc ;
      fflush(stdout);
      trav_time(gcarc,tv,dv,nd,&Ptt,&ierror) ;
      wp_time_window(gcarc, eq->wp_win4, &twp_beg, &twp_end) ;
      t0 = (double)hd_data.o ;
      if (TWPTT)
	t0 += Ptt ;
      else
	t0 += opt->ts;
      n1_data = (int)((t0 + twp_beg - (double)hd_data.b) / ((double)hd_data.delta)) ; /* first data Sample (corrected)  */
      n2_data = n1_data + (int)((twp_end - twp_beg) / ((double)hd_data.delta))      ; /* Last data Sample               */
      npts    = n2_data - n1_data + 1 ;
      if ( (n1_data < 0) || (n2_data >= hd_data.npts) ) 
	{
	  fprintf(stderr,"**** Incomplete data, rejected : %s\n", datafile) ;
	  fprintf( o_log,"**** Incomplete data, rejected : %s\n", datafile) ;
	  fflush(stderr);
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
	  GF = get_gf_filename(gfdirs[j], hd_data.kstnm, hd_data.knetwk, 
			       hd_data.kcmpnm[2], ".SAC.sac.bp") ;
	  strcat(gf_file,GF);
	  free((void*) GF);
	  flag = fill_G(gf_file, datafile, &hd_GF, &hd_data, npts, Ptt, twp_beg, twp_end, 
			tmparray, (*G)[ns][j], opt, o_log) ;
	  if (flag)
	    break ;
	}
      if (flag) /* Error reading GF */
	{
	  for(j=0; j<NM; j++)
	    free((void*)(*G)[ns][j]);
	  free((void**)(*G)[ns]);
	  continue ;
	}
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
      /* Calculate seismogram peak-to-peak and average amplitude */
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
  fclose(i_sac) ;
  /* Memory Freeing */
  for(i=ns;i<*nsac;i++)
    free((void*)(*sacfiles)[i]) ;
  *nsac  = ns ;
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
screen_rms(int *nsac, char **data_name, double **data, double ***G, sachdr *hd_synt, structopt *opt, FILE *o_log)
{
  int j, newn;

  newn = 0 ;
  for (j=0; j<*nsac; j++)
    {
      if( opt->rms_in[j] < opt->th_val )
	{
	  if (o_log != NULL)
	    fprintf( o_log,"stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", hd_synt[j].kstnm, 
		     hd_synt[j].knetwk, hd_synt[j].kcmpnm, hd_synt[j].gcarc, hd_synt[j].az, 
		     hd_synt[j].user[2], hd_synt[j].user[3]) ; 
	  data_name[newn] = data_name[j] ;
	  data[newn]      = data[j]      ;
	  G[newn]         = G[j]         ;
	  hd_synt[newn]   = hd_synt[j]   ;
	  opt->p2p[newn]  = opt->p2p[j]  ;
	  opt->avg[newn]  = opt->avg[j]  ;
	  opt->wgt[newn]  = opt->wgt[j]  ;
	  opt->rms_r[newn]=opt->rms_r[j] ;
	  newn++ ;	  
	}
      else
	{
	  if (o_log != NULL)
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
screen_ratio(int *nsac,char **data_name,double **data,double ***G,sachdr *hd_synt, 
	     structopt *opt, FILE *o_log)
{
  int j, newn ;

  newn = 0 ;
  for (j=0; j<*nsac; j++)
    {
      if (opt->rms_r[j]<opt->rms_r_th && opt->rms_r[j]>1./opt->rms_r_th)
	{
	  fprintf( o_log,"stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", hd_synt[j].kstnm, 
		   hd_synt[j].knetwk, hd_synt[j].kcmpnm, hd_synt[j].gcarc, hd_synt[j].az, 
		   hd_synt[j].user[2], hd_synt[j].user[3]) ; 
	  data_name[newn] = data_name[j] ;
	  data[newn]      = data[j]      ;
	  G[newn]         = G[j]         ;
	  hd_synt[newn]   = hd_synt[j]   ;
	  opt->p2p[newn]  = opt->p2p[j]  ;
	  opt->avg[newn]  = opt->avg[j]  ;
	  opt->wgt[newn]  = opt->wgt[j]  ;
	  opt->rms_r[newn]=opt->rms_r[j] ;
	  newn++ ;	  
	}
      else
	{
	  fprintf(stderr, "**** Rejected trace (rms ratio exceed the threshold): %s\n", 
		  data_name[j]) ;
	  free((void*)data_name[j]) ;
	  free((void*)data[j])      ;
	  free_G(G+j)               ;
	}
    }
  *nsac = newn ;
}


void 
screen_med(int *nsac, char **data_name, double **data, double ***G, 
	   sachdr *hd_synt, structopt *opt, FILE *o_log)
{
  int    j, newn ;
  double min, max, val ;
  
  min = 0.1 * (opt->p2p_med) ;
  max = 3.0 * (opt->p2p_med) ;
  
  if (o_log != NULL)
    {
      fprintf(o_log,"screen_med:\n") ;
      fprintf(o_log,"   p2p_med: %15.8f\n",opt->p2p_med) ;
      fprintf(o_log,"   reject p2p < : %15.8f or > %15.8f\n",min,max) ;
      fprintf(o_log,"   reject avg > : %15.8f \n",opt->p2p_med/2) ;
    }
  newn = 0 ;
  for (j=0;j<*nsac;j++)
    {
      val = opt->p2p[j];
      if ( (min < val) && (val < max) && (fabs(opt->avg[j]) < (opt->p2p_med)/2.) )
	{
	  if (opt->th_val <= 0. && o_log != NULL)
	    fprintf( o_log,"stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", hd_synt[j].kstnm, 
		     hd_synt[j].knetwk, hd_synt[j].kcmpnm, hd_synt[j].gcarc, hd_synt[j].az, 
		     hd_synt[j].user[2], hd_synt[j].user[3]) ; 
	  data_name[newn] = data_name[j] ;
	  data[newn]      = data[j]      ;
	  G[newn]         = G[j]         ;
	  hd_synt[newn]   = hd_synt[j]   ;
	  opt->p2p[newn]  = opt->p2p[j]  ;
	  opt->avg[newn]  = opt->avg[j]  ;
	  opt->wgt[newn]  = opt->wgt[j]  ;
	  opt->rms_r[newn]=opt->rms_r[j] ;
	  newn++ ;
	}
      else 
	{
	  if (o_log != NULL)
	    fprintf(stderr,"**** Rejected trace (p2p or avg out of bounds): %s\n",
		    data_name[j])  ; 
	  free((void*)data_name[j]);
	  free((void*)data[j])     ;
	  free_G(G+j)              ;
	}
    }
 *nsac = newn ;
}


void 
median(int nsac, structopt *opt)
{
  int    i    ;
  double *tmp ;
  
  tmp = double_alloc(nsac) ;
  for (i=0;i<nsac;i++)
    tmp[i] = opt->p2p[i]    ;

  sort(tmp,nsac);

  opt->p2p_med  = tmp[(int)    (nsac/2.)-1] ;
  opt->p2p_low  = tmp[(int)    (nsac/4.)-1] ;
  opt->p2p_high = tmp[(int) (nsac*3./4.)-1] ;

  free((void *)tmp);
}
  

void 
sort(double *tab,int ns)
{
  int    i, j ;
  double swp  ;

  for (i=0;i<ns;i++)
      for(j=i+1;j<ns;j++)
	{
	  if (tab[i] > tab[j]) {
	    swp    = tab[i] ;
	    tab[i] = tab[j] ;
	    tab[j] = swp ;     }
	}
}


void 
w_log_header(char **argv,structopt *opt,str_quake_params *eq,double *wp_win4,FILE *o_log) 
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


void 
load_kernel(str_quake_params *eq,structopt *opt,sachdr *hd_synth,int nsac,
	    int nd,double *dv,double *tv, double ***G,FILE *o_log)
{
  int    i, j, flag, ierror ;
  double gcarc, Ptt, twp_beg, twp_end   ;
  double *tmparray  ;
  char   *gf_file, *GF ;
  char   gfdirs[6][7] = {"gf_rr/","gf_tt/", 
			 "gf_pp/","gf_rt/", 
			 "gf_rp/","gf_tp/"} ;
  sachdr hd_GF ;
  /* Allocating memory */
  gf_file  = char_alloc(FSIZE) ;
  tmparray = double_alloc((int)__LEN_SIG__) ;  
  hdr_alloc(&hd_GF)   ;

  /* Loop on channels */
  for(i=0; i<nsac; i++)
    {
      /* Time Window  */
      gcarc = (double) hd_synth[i].gcarc ;
      trav_time(gcarc,tv,dv,nd,&Ptt,&ierror) ;
      wp_time_window(gcarc, eq->wp_win4, &twp_beg, &twp_end) ;
      /* read GF */
      for(j=0; j<NM; j++)  
	{
	  strcpy(gf_file,eq->gf_dir);
	  GF = get_gf_filename(gfdirs[j], hd_synth[i].kstnm, 
			       hd_synth[i].knetwk, hd_synth[i].kcmpnm[2], 
			       ".SAC.sac.bp") ;
	  strcat(gf_file,GF);
	  free((void*) GF);
	  flag = fill_G(gf_file,gf_file,&hd_GF,hd_synth+i,hd_synth[i].npts,Ptt,twp_beg,twp_end, 
			tmparray,G[i][j],opt,o_log) ;
	  if (flag)
	    {
	      fflush(stdout);
	      fprintf(stderr,"*** load_kernel (wpinversion_xy_gs): Error while reading %s\n",gf_file);
	      fflush(stderr);
	      exit(1);
	    }
	}
    }

  /* Memory Freeing */
  free((void*)tmparray) ;
  free((void*)gf_file)  ;
}

void 
realloc_gridsearch(int nsac, double **rms, double *global_rms, 
		   double ***dcalc,int flag)
{
  int i, j ; 
  for(i=0 ; i<nsac ; i++)  /* Free predicted data vector */
    {
      for(j=0 ; j<flag ; j++)
	free((void*)dcalc[i][j]) ;
      free((void**)dcalc[i]);
    }
  for(j=0 ; j<2*flag ; j++) /* Initialize rms vectors      */
    {
      global_rms[j] = 0. ;
      for(i=0 ; i<nsac ; i++)
	rms[i][j] = 0.   ;
    }
}




void
fast_ts_gridsearch(int nsac,int M,int nd,double *dv,double *tv,sachdr *hd_synt,double **data, 
		   double ***G,double ***dcalc,double **rms,double *global_rms,structopt *opt,
		   str_quake_params *eq,double *tsopt,double *rmsopt, FILE *o_log)
{
  int    i,it, k, Nexp, ref_flag;
  double Err, Cond,dt,ts,dtmin,dtmax,rmsini,tsopt2,rmsopt2,sdrM0[4]; 
  FILE   *tmp, *o_gs   ;
  printf("FAST CENTROID TIME DELAY GRID SEARCH\n") ;
  fprintf(o_log,"Fast time-shift grid search is enabled (output : %s)\n",opt->tsgsfile) ;
  /* Initialize variables */
  dtmin = opt->dts_min - eq->ts ;
  dtmax = opt->dts_max - eq->ts ;
  dt    = opt->dts_step ;
  it   = 0 ;
  k    = 0 ;
  Nexp = 0 ;
  ref_flag = opt->ref_flag ;
  opt->ref_flag = 0;
  *tsopt  = 0. ;
  tsopt2  = dtmin  ;
  rmsini  = 1000.*(global_rms[0]) ;
  *rmsopt = rmsini ;
  rmsopt2 = 1.1e10 ;
  tmp     = openfile_wt("_tmp_ts_table") ;
  while( it < opt->ts_Nit )
    {
      opt->dts_val = dtmin ;
      while ( opt->dts_val <= dtmax )
	{	 
	  /* Compute inversion for opt->dts_val */
	  load_kernel(eq,opt,hd_synt,nsac,nd,dv,tv,G,o_log)      ;
	  if (opt->dc_flag) /* Double Couple inversion                 */
	    {               /* Warning: This has not been fully tested */
	      for(i=0;i<4;i++)
		sdrM0[i] = opt->priorsdrM0[i] ;   	      
	      inversion_dc(nsac,hd_synt,G,data,sdrM0,global_rms,opt,NULL) ;
	      sdr2mt(eq->vm[0],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2]) ;
	    }
	  else
	    {
	      /* Free memory */
	      realloc_gridsearch(nsac, rms, global_rms, dcalc,1) ;	      
	      inversion(M,nsac,hd_synt,G,data,eq->vm[0],&Cond,opt,NULL) ;
	      calc_data(nsac,hd_synt,G,eq->vm, data,dcalc,opt,NULL) ;   
	      calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms, opt) ;
	      /* Get RMS error */
	    }
	  Err = 1000.*(global_rms[0]) ;
	  ts  = eq->ts+opt->dts_val    ;
	  fprintf(tmp,"%03d %03d %10.4f %10.4f %10.4f %10.4f %10.4f %12.8f %12.8f\n",
		  k,it,ts,eq->hd,eq->evla,eq->evlo,eq->evdp,Err,(global_rms[0])/(global_rms[1])) ;
	  printf("        ts = %10.4f rms = %12.8f mm\n",ts, Err) ;
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
      printf("Optimum values: time_shift = %10.4f rms = %12.8f mm \n",eq->ts+*tsopt,*rmsopt) ;
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
  opt->ref_flag = ref_flag ;
  /* Write output file */
  o_gs = openfile_wt(opt->tsgsfile) ;
  fprintf( o_gs,"%10.4f %12.8f\n",eq->ts+*tsopt,*rmsopt);
  fprintf( o_gs,"%10.4f %12.8f\n",eq->ts,rmsini);
  tmp  = openfile_rt("_tmp_ts_table", &k) ;
  while ((k = getc(tmp)) != EOF)
    putc(k,o_gs) ;
  fclose( tmp) ;
  fclose(o_gs) ;
}



void
run_ts_gs(double *ts,int Ngrid, int nsac, int M, int nd, double *dv,
	  double *tv, sachdr *hd_synt, double **data, str_quake_params *eq,
	  structopt *opt, double **vrms,int verbose)
{
  int i,j;
  double Cond,***dcalc,**rms,*global_rms,*sdrM0,***G=NULL ;
  str_quake_params eq_gs ;  
#ifdef _OPENMP      
  int rang,ntaches;
  rang = omp_get_thread_num();
  ntaches = omp_get_num_threads();
#endif
  sdrM0 = NULL;
  /* Memory allocation */
  G = double_alloc3p( nsac ) ; 
  for(i=0;i<nsac;i++)
    G[i] = double_alloc2(NM,hd_synt[i].npts) ; 
  dcalc = double_alloc3p(nsac)       ;  
  if (opt->dc_flag)
    sdrM0    = double_alloc(4)  ;
  rms = double_calloc2(nsac,2) ;
  global_rms    = double_calloc(2)     ;
  eq_gs.wp_win4 = double_alloc(4)      ;   
  eq_gs.vm      = double_calloc2(2,NM) ;
  copy_eq(eq,&eq_gs) ;  
  /* Main loop */
  for(i=0;i<Ngrid;i++)
    {
      eq_gs.ts = ts[i];
      eq_gs.hd = ts[i];
      /* Compute G */
      calc_kernel(&eq_gs,opt,hd_synt,nsac,"l",nd,dv,tv,G,NULL) ;
      /* Inversion */
       if (opt->dc_flag) /* Double Couple inversion                 */
	{                /* Warning: This has not been fully tested */
	  for(j=0;j<4;j++)
	    sdrM0[j] = opt->priorsdrM0[j] ;      	  
	  inversion_dc(nsac,hd_synt,G,data,sdrM0,vrms[i],opt,NULL) ;
	  sdr2mt(eq_gs.vm[0],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2]) ;
	}
      else
	{
	  inversion(M,nsac,hd_synt,G,data,eq_gs.vm[0],&Cond,opt,NULL) ;
	  /* Predicted data  */
	  calc_data(nsac,hd_synt,G,eq_gs.vm,data,dcalc,opt,NULL) ;
	  /* RMS error       */
	  calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms,opt) ;
	  vrms[i][0] = global_rms[0] ;
	  vrms[i][1] = global_rms[1] ;
	  realloc_gridsearch(nsac, rms, global_rms, dcalc,1) ;
	}
       fflush(stderr);
#ifdef _OPENMP      
      if (verbose)
	printf("thread %d/%d %10.4f %12.8f %12.8f\n",rang+1,ntaches,
	       ts[i],vrms[i][0]*1000,vrms[i][0]/vrms[i][1]);
#else
      if (verbose)
	printf("%10.4f %12.8f %12.8f\n",ts[i],vrms[i][0]*1000,vrms[i][0]/vrms[i][1]);	
#endif
       fflush(stdout); 
    }
  /* Free memory */
  free((void*)eq_gs.wp_win4) ;
  free((void*)eq_gs.vm[0])   ;
  free((void*)eq_gs.vm[1])   ;
  free((void**)eq_gs.vm)     ;
  if (opt->dc_flag)
    free((void*)sdrM0);
  for(i=0 ; i<nsac ; i++)
    {
      free((void*)rms[i] ) ;
      free_G(&G[i]);
    }
  free((void***)dcalc)    ;
  free((void**)rms )      ;
  free((void*)global_rms) ;  
}

void
ts_gridsearch(int nsac,int M,int nd,double *dv,double *tv,sachdr *hd_synt,double **data, 
	      double *rmsini,structopt *opt,str_quake_params *eq,double *tsopt,
	      double *rmsopt, FILE *o_log)
{
  int    i=0,j=0,k=0,it=0,Nexp=0,ref_flag,Ngrid,MaxNgrid;
  double Err,NErr,ts,dt,tsmin,tsmax,tsopt2,rmsopt2      ; 
  double *new_ts,**vrms ;
  FILE   *tmp, *o_gs     ;

  printf("CENTROID TIMING GRID SEARCH\n") ;
  fprintf(o_log,"Centroid timing grid-search is enabled (output : %s)\n",opt->tsgsfile) ;
  /* Initialize variables */
  tsmin = opt->dts_min  ;
  tsmax = opt->dts_max  ;
  dt    = opt->dts_step ;  
  ref_flag = opt->ref_flag ;
  opt->ref_flag = 0;
  *tsopt  = eq->ts ;
  tsopt2  = tsmin  ;
  *rmsopt = 1000.*rmsini[0] ;
  rmsopt2 = 1.1e10 ;
  /* Memory allocation */
  MaxNgrid = (int)((tsmax-tsmin+20*dt)*pow(2,(double)opt->ts_Nit)/dt) ;
  new_ts = double_alloc(MaxNgrid)    ; /* Timings to be explored */
  vrms   = double_alloc2(MaxNgrid,2) ;
  /* Main loop */
  k = 0;
  tmp     = openfile_wt("_tmp_ts_table") ;
  while( it < opt->ts_Nit )
    {
      if (!Nexp)
	printf("Iteration %d\n",it+1); 	      
      ts = tsmin; Ngrid=0;
      while ( ts <= tsmax )
	{
	  new_ts[Ngrid] = ts;
	  Ngrid        += 1 ;
	  ts += dt;
	}
      /* Compute vrms for each new_ts */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(j)
      {
#pragma omp for schedule(dynamic)
	for(j=0;j<Ngrid;j++)
	  run_ts_gs(new_ts+j,1,nsac,M,nd,dv,tv,hd_synt,data,eq,opt,vrms+j,1) ;
      }
#else
      run_ts_gs(new_ts,Ngrid,nsac,M,nd,dv,tv,hd_synt,data,eq,opt,vrms,1) ;
#endif      
      /* Find optimal timings, Fill history */
      for(i=0;i<Ngrid;i++)
	{
	  Err  = 1000*vrms[i][0]       ;
	  NErr = vrms[i][0]/vrms[i][1] ;
	  fprintf(tmp,"%03d %03d %10.4f %10.4f %10.4f %10.4f %10.4f %12.8f %12.8f\n",
		  k,it,new_ts[i],new_ts[i],eq->evla,eq->evlo,eq->evdp,Err,NErr) ;
	  printf("        ts = hd = %10.4f rms = %12.8f mm\n",new_ts[i], Err) ;
	  /* Best timings */
	  if (Err < *rmsopt)
	    {
	      *rmsopt = Err ;
	      *tsopt  = new_ts[i] ;
	    }
	  else if (Err < rmsopt2)
	    {
	      rmsopt2 = Err ;
	      tsopt2  = new_ts[i] ;
	    }
	  k++;
	}
      printf("Optimum values: time_shift = half_duration = %10.4f ; rms = %12.8f mm \n",*tsopt,*rmsopt) ;
      if (tsmax <= *tsopt + dt && Nexp < 5)
	{
           printf("Optimum value on the maximum explored time-shift (=half-duration)\n   ... extending the timing grid-search area\n");
	   tsmin = *tsopt + dt   ;
	   tsmax = *tsopt + 4*dt ;
	   Nexp++;
	   continue ;
	}
      if (it)
	dt = dt/2. ;
      if (tsopt2 <= *tsopt)
	{
	  tsmin = tsopt2   - dt/2. ;
	  tsmax = (*tsopt) + dt/2. ;
	}
      else if (tsopt2 > *tsopt)
	{
	  tsmin = (*tsopt) - dt/2. ;
	  tsmax = tsopt2   + dt/2. ;
	}
      if (tsmin < 1. )
	while (tsmin <  1. )
	  tsmin += dt ;	      
      Nexp=0;
      it++;
    }
  fclose(tmp);
  opt->ref_flag = ref_flag ;
  /* Write output file */
  o_gs = openfile_wt(opt->tsgsfile) ;
  fprintf( o_gs,"%10.4f %12.8f\n",*tsopt,*rmsopt);
  fprintf( o_gs,"%10.4f %12.8f\n",eq->ts,1000.*rmsini[0]);
  tmp  = openfile_rt("_tmp_ts_table", &k) ;
  while ((k = getc(tmp)) != EOF)
    putc(k,o_gs) ;
  fclose( tmp) ;
  fclose(o_gs) ;
  /* Free memory */
  
}



int 
find_coor(double **coords,int Ngrid, double lat,double lon,double dep)
{
  int i;
  for(i=0;i<Ngrid;i++)
    if( (int)round(coords[i][0]*100.)==(int)round(lat*100.)
	&& (int)round(coords[i][1]*100.)==(int)round(lon*100.)
	&& (int)round(coords[i][2]*100.)==(int)round(dep*100.) )
      return 1 ;     
  return 0 ;
}

int
find_dep(double *depths, float cmt_dep, int nd)
{
  int jd, index=0 ;
  double best_depth = depths[0];
  for (jd=1; jd<nd; jd++)
    if(fabs(depths[jd]-cmt_dep)<fabs(best_depth-cmt_dep)) 
      {
	best_depth = depths[jd];
	index = jd;
      }
  return index;
}


void
search_emptyedges(double **coords,int Ngrid,double lat,double lon,double dep,
		  double dlat,double dlon,double **extcoords,int *Next)
{
  int i,N=8;
  double lats[] = {lat+dlat,lat+dlat,lat+dlat,lat,lat-dlat,lat-dlat,lat-dlat,lat} ;
  double lons[] = {lon-dlon,lon,lon+dlon,lon+dlon,lon+dlon,lon,lon-dlon,lon-dlon} ;      
  for(i=0;i<N;i++)
    {
      if(!find_coor(coords,Ngrid,lats[i],lons[i],dep) && 
	 !find_coor(extcoords,*Next,lats[i],lons[i],dep))
	{	  
	  FILLCOOR(extcoords[(*Next)],lats[i],lons[i],dep)
	  (*Next)++;
	}
    }
}

void
set_coor(double lat, double lon,double dep,int *Ngrid,double **coords)
{
  if (!find_coor(coords,*Ngrid,lat,lon,dep))
    {
      FILLCOOR(coords[*Ngrid],lat,lon,dep)
      (*Ngrid)++;
    }
}


void 
fill_grid(structopt *opt, double *LL, double dep, int *Ngrid, double **new_loc)
{
  int i,j;
  for(i=0;i<opt->xy_Nx*2+1;i++)
    for(j=0;j<opt->xy_Nx*2+1;j++)
      set_coor(LL[0]+opt->dlat*(double)i,LL[1]+opt->dlon*(double)j,dep,Ngrid,new_loc);
}


void 
fill_arround(double *loc,double dlat,double dlon,double dep,int *Ngrid,double **new_loc)
{
  set_coor(loc[0]+dlat,loc[1]-dlon,dep,Ngrid,new_loc) ;
  set_coor(loc[0]+dlat,loc[1]     ,dep,Ngrid,new_loc) ;
  set_coor(loc[0]+dlat,loc[1]+dlon,dep,Ngrid,new_loc) ;
  set_coor(loc[0]     ,loc[1]+dlon,dep,Ngrid,new_loc) ;
  set_coor(loc[0]-dlat,loc[1]+dlon,dep,Ngrid,new_loc) ;
  set_coor(loc[0]-dlat,loc[1]     ,dep,Ngrid,new_loc) ;
  set_coor(loc[0]-dlat,loc[1]-dlon,dep,Ngrid,new_loc) ;
  set_coor(loc[0]     ,loc[1]-dlon,dep,Ngrid,new_loc) ;
}


void
copy_eq(str_quake_params *i_eq, str_quake_params *o_eq)
{
  strcpy(o_eq->evnm,i_eq->evnm)       ;
  strcpy(o_eq->cmtfile,i_eq->cmtfile) ;
  memcpy(o_eq->wp_win4,i_eq->wp_win4,4*sizeof(double));
  o_eq->dmin = i_eq->dmin ;
  o_eq->dmax = i_eq->dmax ;
  o_eq->filtorder = i_eq->filtorder ;
  o_eq->filtnpass = i_eq->filtnpass ;
  o_eq->flow      = i_eq->flow  ;
  o_eq->fhigh     = i_eq->fhigh ;
  o_eq->idtr      = i_eq->idtr  ;
  o_eq->preevent  = i_eq->preevent ;
  o_eq->fend    = i_eq->fend    ;
  o_eq->fend    = i_eq->fend    ;
  strcpy(o_eq->pdeline,i_eq->pdeline);
  o_eq->ot_ye = i_eq->ot_ye       ;
  o_eq->ot_mo = i_eq->ot_mo       ;
  o_eq->ot_dm = i_eq->ot_dm       ;
  o_eq->ot_ho = i_eq->ot_ho       ; 
  o_eq->ot_mi = i_eq->ot_mi       ;
  o_eq->ot_se = i_eq->ot_se       ;
  o_eq->ot_ms = i_eq->ot_se       ;
  o_eq->pde_evla = i_eq->pde_evla ;
  o_eq->pde_evlo = i_eq->pde_evlo ;
  o_eq->pde_evdp = i_eq->pde_evdp ;
  strcpy(o_eq->evid,i_eq->evid)   ;
  o_eq->ts = i_eq->ts ;
  o_eq->hd = i_eq->hd ;
  o_eq->evla = i_eq->evla ;
  o_eq->evlo = i_eq->evlo ;
  o_eq->evdp = i_eq->evdp ;
}


void
copy_opt(structopt *i_opt, structopt *o_opt,int nsac)
{
  o_opt->ts_Nit   = i_opt->ts_Nit   ;
  o_opt->ps       = i_opt->ps       ;
  o_opt->xy_Nx    = i_opt->xy_Nx    ;
  o_opt->xy_Nopt  = i_opt->xy_Nopt  ;
  o_opt->xy_Nit   = i_opt->xy_Nit   ;
  o_opt->hdsafe   = i_opt->hdsafe   ;
  o_opt->dc_flag  = i_opt->dc_flag  ;
  o_opt->ref_flag = i_opt->ref_flag ;
  o_opt->ip       = i_opt->ip       ;
  memcpy(o_opt->ib,i_opt->ib,4*sizeof(int));
  o_opt->th_val   = i_opt->th_val   ;
  o_opt->cth_val  = i_opt->cth_val  ;
  o_opt->df_val   = i_opt->df_val   ;
  o_opt->med_val  = i_opt->med_val  ;
  o_opt->op_pa    = i_opt->op_pa    ;
  o_opt->p2p_med  = i_opt->p2p_med  ;
  o_opt->p2p_low  = i_opt->p2p_low  ;
  o_opt->p2p_high = i_opt->p2p_high ;
  o_opt->dts_min  = i_opt->dts_min  ;
  o_opt->dts_step = i_opt->dts_step ;
  o_opt->dts_max  = i_opt->dts_max  ;
  o_opt->dts_val  = i_opt->dts_val  ;
  o_opt->rms_r_th = i_opt->rms_r_th ;
  memcpy(o_opt->rms_in,i_opt->rms_in,nsac*sizeof(double)) ;
  memcpy(o_opt->rms_r,i_opt->rms_r,nsac*sizeof(double))   ;
  memcpy(o_opt->p2p,i_opt->p2p,nsac*sizeof(double))       ;
  memcpy(o_opt->avg,i_opt->avg,nsac*sizeof(double))       ;
  memcpy(o_opt->wgt,i_opt->wgt,nsac*sizeof(double))       ;
  o_opt->xy_dx    = i_opt->xy_dx    ;
  o_opt->dlat     = i_opt->dlat     ;
  o_opt->dlon     = i_opt->dlon     ;
  o_opt->dz       = i_opt->dz       ;
  o_opt->mindep   = i_opt->mindep   ;
  o_opt->ntr_val  = i_opt->ntr_val  ;  
  o_opt->wZ       = i_opt->wZ       ;
  o_opt->wL       = i_opt->wL       ;
  o_opt->wT       = i_opt->wT       ;
  o_opt->azp      = i_opt->azp      ;   
  o_opt->dmin     = i_opt->dmin     ;
  o_opt->dmax     = i_opt->dmax     ;
  memcpy(o_opt->priorsdrM0,i_opt->priorsdrM0,4*sizeof(double));
  strcpy(o_opt->i_master ,i_opt->i_master)  ;
  strcpy(o_opt->i_saclst ,i_opt->i_saclst)  ;   
  strcpy(o_opt->o_saclst ,i_opt->o_saclst)  ;
  strcpy(o_opt->log      ,i_opt->log)       ;
  strcpy(o_opt->o_cmtf   ,i_opt->o_cmtf)    ;
  strcpy(o_opt->p_data   ,i_opt->p_data)    ;
  strcpy(o_opt->osacdir  ,i_opt->osacdir)   ;
  strcpy(o_opt->psfile   ,i_opt->psfile)    ;   
  strcpy(o_opt->wpbmfile ,i_opt->wpbmfile)  ;
  strcpy(o_opt->refbmfile,i_opt->refbmfile) ;
  strcpy(o_opt->tsgsfile ,i_opt->tsgsfile)  ;
  strcpy(o_opt->xygsfile ,i_opt->xygsfile)  ;
  strcpy(o_opt->o_covf   ,i_opt->o_covf)    ;
}


void 
extend_xy_gs(double **old_loc,int Nold_loc,double *gfdep,int ngfdep,double **opt_loc,
	     structopt *opt,str_quake_params *eq,int it,int Nexp,double *LL,double *UR,
	     int *IZ, double **new_loc,int *Nnew_loc)
{
  int i,j,k,l,Next=0,flag ;
  double lat,lon;
  double NLL[]={LL[0],LL[1]},NUR[]={UR[0],UR[1]};

  if (opt->xy_Nit>1) /* Search isolated best locations */
    {
      search_emptyedges(old_loc,Nold_loc,eq->evla,eq->evlo,eq->evdp,
			opt->dlat,opt->dlon,new_loc,&Next) ;
      for(i=0;i<opt->xy_Nopt;i++)
	search_emptyedges(old_loc,Nold_loc,opt_loc[i][0],opt_loc[i][1],
			  opt_loc[i][2],opt->dlat,opt->dlon,new_loc,&Next);
    }
  
  if (!it) 
    {
      if (opt->dz > 0. && !Nexp) /* Extending depth   */
	{
	  k = find_dep(gfdep,opt->mindep,ngfdep);
	  for(i=0;i<opt->xy_Nopt;i++)
	    {	      
	      l = find_dep(gfdep,opt_loc[i][2],ngfdep);
	      if (!i)      
		{
		  flag = 0;
		  if (l<ngfdep-1 && l>=IZ[1]-2) /* Extend downward */
		    {
		      IZ[0] = IZ[1]+1;
		      IZ[1] = find_dep(gfdep,gfdep[IZ[1]]+opt->dz,ngfdep);		      		      
		      flag++;
		    }
		  else if (l>k && l<=IZ[0]+2)       /* Extend upward   */
		    {	
		      IZ[1] = IZ[0]-1;
		      IZ[0] = find_dep(gfdep,gfdep[IZ[0]]-opt->dz,ngfdep);
		      if (IZ[0]<k)
			IZ[0] = k;		      
		      flag++;
		    }
		  if(flag)
		    for(j=IZ[0];j<=IZ[1];j++)
		      fill_grid(opt,LL,gfdep[j],&Next,new_loc) ;
		}	      
	      if(l>k)                          /* Extend bellow optimums */
		fill_grid(opt,LL,gfdep[l-1],&Next,new_loc) ;
	      if(l<ngfdep-1)                   /* Extend above optimums  */
		fill_grid(opt,LL,gfdep[l+1],&Next,new_loc) ;	      
	    }
	}
      else if (Next)    /* Extending lat-lon */
	{
	  for(i=0;i<Next;i++)
	    {
	      if (new_loc[i][0]<NLL[0])
		NLL[0]=new_loc[i][0];
	      if (new_loc[i][1]<NLL[1])
		NLL[1]=new_loc[i][1];
	      if (new_loc[i][0]>NUR[0])
		NUR[0]=new_loc[i][0];
	      if (new_loc[i][1]>NUR[1])
		NUR[1]=new_loc[i][1];
	    }
	  k = (int)round((NUR[0]-NLL[0])/opt->dlat);
	  l = (int)round((NUR[1]-NLL[1])/opt->dlon);
	  for(i=0;i<=k;i++)
	    for(j=0;j<=l;j++)
	      {
		lat = NLL[0]+opt->dlat*(double)i;
		lon = NLL[1]+opt->dlon*(double)j;
		if (!find_coor(old_loc,Nold_loc,lat,lon,eq->evdp) && 
		    !find_coor(new_loc,Next,lat,lon,eq->evdp))
		  {
		    FILLCOOR(new_loc[Next],lat,lon,eq->evdp)
		    Next++;
		  }
		
	      }
	}
    }
  *Nnew_loc = 0;
  if (Next)
    {
      printf(" ... extending the spatial grid-search area ... \n");
      *Nnew_loc = Next;
    }
  LL[0] = NLL[0] ;
  LL[1] = NLL[1] ;
  UR[0] = NUR[0] ;
  UR[1] = NUR[1] ;
}   


void
run_xy_gs(double **coords,int Ngrid,int nsac,int M,int nd,double *dv,
	  double *tv,sachdr *hd_synt,double **data,str_quake_params *eq,
	  structopt *opt,double **vrms,double **vms,int verbose)
{
  int i,j;
  double Cond,***dcalc,**rms,*global_rms,*sdrM0,***G=NULL ;
  structopt opt_gs       ;
  str_quake_params eq_gs ;  
#ifdef _OPENMP      
  int rang,ntaches;
  rang = omp_get_thread_num();
  ntaches = omp_get_num_threads();
#endif
  sdrM0 = NULL;
  /* Memory allocation */
  G = double_alloc3p( nsac ) ; 
  for(i=0;i<nsac;i++)
    G[i] = double_alloc2(NM,hd_synt[i].npts) ; 
  dcalc = double_alloc3p(nsac)       ;  
  if (opt->dc_flag)
    sdrM0    = double_alloc(4)  ;
  rms = double_calloc2(nsac,2) ;
  global_rms    = double_calloc(2)     ;
  eq_gs.wp_win4 = double_alloc(4)      ;   
  eq_gs.vm      = double_calloc2(2,NM) ;
  copy_eq(eq,&eq_gs) ;    
  opt_gs.rms_in = double_alloc(nsac) ;
  opt_gs.rms_r  = double_alloc(nsac) ;
  opt_gs.p2p    = double_alloc(nsac) ;
  opt_gs.avg    = double_alloc(nsac) ;
  opt_gs.wgt    = double_alloc(nsac) ;
  copy_opt(opt,&opt_gs,nsac) ;
  /* Main loop */
  for(i=0;i<Ngrid;i++)
    {
      eq_gs.evla = coords[i][0] ;
      eq_gs.evlo = coords[i][1] ;
      eq_gs.evdp = coords[i][2] ;
      /* Compute G       */
      calc_kernel(&eq_gs,&opt_gs,hd_synt,nsac,"l",nd,dv,tv,G,NULL) ;
      /* Inversion       */
      if (opt->dc_flag) /* Double Couple inversion                 */
	{               /* Warning: This has not been fully tested */
	  for(j=0;j<4;j++)
	    sdrM0[j] = opt_gs.priorsdrM0[j] ;      
	  inversion_dc(nsac,hd_synt,G,data,sdrM0,vrms[i],&opt_gs,NULL) ;
	  sdr2mt(vms[i],sdrM0[3],sdrM0[0],sdrM0[1],sdrM0[2]) ;
	}
      else
	{
	  inversion(M,nsac,hd_synt,G,data,vms[i],&Cond,&opt_gs,NULL) ;
	  /* Predicted data  */
	  calc_data(nsac,hd_synt,G,vms+i,data,dcalc,&opt_gs,NULL) ;
	  /* RMS error       */
	  calc_rms(nsac,hd_synt,data,dcalc,rms,global_rms,&opt_gs) ;
	  vrms[i][0] = global_rms[0] ;
	  vrms[i][1] = global_rms[1] ;
	  realloc_gridsearch(nsac, rms, global_rms, dcalc,1) ;
	}
      fflush(stderr);
#ifdef _OPENMP      
      if (verbose)
	printf("thread %d/%d %10.4f %10.4f %10.4f %12.8f %12.8f\n",rang+1,
	       ntaches,coords[i][0],coords[i][1],coords[i][2],vrms[i][0]*1000,vrms[i][0]/vrms[i][1]);
#else
      if (verbose)
	printf("%10.4f %10.4f %10.4f %12.8f %12.8f\n",coords[i][0],coords[i][1],coords[i][2],
	       vrms[i][0]*1000,vrms[i][0]/vrms[i][1]);      
#endif
      fflush(stdout);      
    }
  /* Free memory */
  free((void*)eq_gs.wp_win4) ;
  free((void*)eq_gs.vm[0])   ;
  free((void*)eq_gs.vm[1])   ;
  free((void**)eq_gs.vm)     ;
  free((void*)opt_gs.rms_in) ;
  free((void*)opt_gs.rms_r)  ;
  free((void*)opt_gs.p2p)    ;
  free((void*)opt_gs.avg)    ;
  free((void*)opt_gs.wgt)    ;
  if (opt->dc_flag)
    free((void*)sdrM0);
  for(i=0 ; i<nsac ; i++)
    {
      free((void*)rms[i] ) ;
      free_G(&G[i]);
    }
  free((void***)dcalc)    ;
  free((void**)rms )      ;
  free((void*)global_rms) ;
}

void
xy_gridsearch(int nsac,int M, int nd,double *dv,double *tv, sachdr *hd_synt,
	      double **data,double ***G,double ***dcalc,double **rms,double *global_rms,
	      structopt *opt,str_quake_params *eq,double *rmsopt,double *latopt,double *lonopt,
	      double *depopt,FILE *o_log)
{
  int i,j,k,it=0,ncel=0,Nexp=0,MaxNexp=5,Ngrid,Ngridini,MaxNgrid,flag=0,ref_flag;
  int    l,IZ[2],ndep, n2dep ;
  double Z[2], gfdep[NDEPTHS_GF] ;
  double rmsini,LL[2],UR[2],err,Nerr,wx         ;
  double **new_loc,**vrms,**vms,**opt_loc, **old_loc ;
  char   gf_path[FSIZE];
  FILE   *tmp,*o_gs ;
  
  printf("\nCENTROID POSITION GRID SEARCH\n") ;
  fprintf(o_log,"centroid position grid search is enabled (output : %s)\n",opt->xygsfile) ;  
  /* Init lat/lon window */
  wx = opt->xy_dx*(double)opt->xy_Nx ;
  LL[0] = eq->evla-wx; LL[1] = eq->evlo-wx/cos(eq->evla*DEG2RAD) ;
  UR[0] = eq->evla+wx; UR[1] = eq->evlo+wx/cos(eq->evla*DEG2RAD) ;
  opt->dlat = opt->xy_dx;
  opt->dlon = opt->xy_dx/cos(eq->evla*DEG2RAD) ;
  /* Init Depths         */
  Z[0] = eq->evdp-opt->dz; Z[1] = eq->evdp+opt->dz ;
  if(getenv("GF_PATH")) 
    strncpy(gf_path, getenv("GF_PATH"), FSIZE);
  else
    {
      fprintf(stderr, "Error: GF_PATH environment variable not defined\n");
      exit(1);
    }
  get_depths(gf_path,gfdep,&ndep);
  sort(gfdep,ndep);
  IZ[0] = find_dep(gfdep,Z[0],ndep);  
  IZ[1] = find_dep(gfdep,Z[1],ndep);
  k     = find_dep(gfdep,eq->evdp,ndep);
  l     = find_dep(gfdep,opt->mindep,ndep);
  eq->evdp = gfdep[k]; /* Re-define centroid depth */
  if (IZ[0]<l)
    IZ[0] = l;
  if ((k%2) != (IZ[0]%2)) 
    IZ[0]+= 1;
  if (IZ[1]<=IZ[0]) 
    {
      IZ[1]   = IZ[0] ;
      opt->dz = 0     ;
    }
  n2dep = (IZ[1]-IZ[0])/2+1 ;
  printf("Depth range %.4f - %.4f\n",gfdep[IZ[0]],gfdep[IZ[1]]);
  /* Init Grid              */
  Ngridini = 4*n2dep*(opt->xy_Nx*opt->xy_Nx+opt->xy_Nx)+n2dep ; /* Grid size for the first iteration */
  MaxNgrid = 5*Ngridini*(MaxNexp+1)*opt->xy_Nit+opt->xy_Nopt*8*opt->xy_Nit ; /* Max number of points */
  fflush(stderr);
  new_loc = double_alloc2(MaxNgrid,3) ; /* Locations to be explored */
  old_loc = double_alloc2(MaxNgrid,4) ;   /* Location History         */  
  vrms = double_alloc2(MaxNgrid,2)  ;
  vms  = double_alloc2(MaxNgrid,NM) ;
  Ngrid = 0;
  for(i=IZ[0];i<=IZ[1];i+=2) /* Fill grid */
    fill_grid(opt,LL,gfdep[i],&Ngrid,new_loc) ;
  /* Init optimum locations */
  opt_loc = double_alloc2(opt->xy_Nopt,4) ; /* Best locations */
  for(i=0;i<opt->xy_Nopt;i++)
    {
      FILLCOOR2(opt_loc[i],0.,0.,0.,1.e+10)
    }
  rmsini = 1000.*(global_rms[0]) ;
  FILLCOOR2(opt_loc[0],eq->evla,eq->evlo,eq->evdp,rmsini)
  ref_flag = opt->ref_flag ;
  opt->ref_flag = 0        ;
  tmp = openfile_wt("_tmp_xy_table") ;
  while( it < opt->xy_Nit )
    {      
      if (!Nexp)
	{
	  printf("Iteration %d\n",it+1); 
	  if (it) 
	    {
	      opt->dlat /= 2. ;
	      opt->dlon /= 2. ;
	      Ngrid = 0;
	      for(i=0;i<opt->xy_Nopt;i++) 
		{
		  /* Search around best locations   */
		  fill_arround(opt_loc[i],opt->dlat,opt->dlon,opt_loc[i][2],&Ngrid,new_loc) ;
		  if (opt->dz>0.)
		    {
		      k = find_dep(gfdep,opt_loc[i][2],ndep);
		      if(k>l)      /* search bellow */
			{
			  fill_arround(opt_loc[i],opt->dlat,opt->dlon,gfdep[k-1],&Ngrid,new_loc) ;
			  if (!find_coor(old_loc,ncel,opt_loc[i][0],opt_loc[i][1],gfdep[k-1]))
			    set_coor(opt_loc[i][0],opt_loc[i][1],gfdep[k-1],&Ngrid,new_loc) ;
			}
		      if(k<ndep-1) /* search above  */
			{
			  fill_arround(opt_loc[i],opt->dlat,opt->dlon,gfdep[k+1],&Ngrid,new_loc) ;
			  if (!find_coor(old_loc,ncel,opt_loc[i][0],opt_loc[i][1],gfdep[k+1]))
			    set_coor(opt_loc[i][0],opt_loc[i][1],gfdep[k-1],&Ngrid,new_loc) ;
			}
		    }
		}
	    }
	}
      /* Compute vrms and vms for each location in new_loc */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(j)
      {
#pragma omp for schedule(dynamic)
	for(j=0;j<Ngrid;j++)
	  run_xy_gs(new_loc+j,1,nsac,M,nd,dv,tv,hd_synt,data,eq,opt,vrms+j,vms+j,1) ;
      }
#else
      run_xy_gs(new_loc,Ngrid,nsac,M,nd,dv,tv,hd_synt,data,eq,opt,vrms,vms,1) ;
#endif

      /* Find optimal points, Fill location history */
      for(i=0;i<Ngrid;i++)
	{
	  err  = 1000*vrms[i][0]       ;
	  Nerr = vrms[i][0]/vrms[i][1] ;
	  /* Write output file */
	  printf(     "%03d %03d %10.4f %10.4f %10.4f %10.4f %10.4f %12.8f %12.8f\n",
		      ncel,it,eq->ts,eq->hd,new_loc[i][0],new_loc[i][1],new_loc[i][2],err,Nerr) ;
	  fprintf(tmp,"%03d %03d %10.4f %10.4f %10.4f %10.4f %10.4f %12.8f %12.8f\n",
		  ncel,it,eq->ts,eq->hd,new_loc[i][0],new_loc[i][1],new_loc[i][2],err,Nerr) ;
	  /* Best locations    */
	  for(j=0;j<opt->xy_Nopt;j++) 
	    {
	      if (find_coor(opt_loc,opt->xy_Nopt,new_loc[i][0],new_loc[i][1],new_loc[i][2]))
		break;
	      if (err < opt_loc[j][3]) 
		{
		  for(k=opt->xy_Nopt-1;k>j;k--)
		    {
		      FILLCOOR2(opt_loc[k],opt_loc[k-1][0],opt_loc[k-1][1],opt_loc[k-1][2],opt_loc[k-1][3])
		    }
		  FILLCOOR2(opt_loc[j],new_loc[i][0],new_loc[i][1],new_loc[i][2],err)
		  break;
		}
	    }
	  /* Fill location history */
	  FILLCOOR2(old_loc[ncel],new_loc[i][0],new_loc[i][1],new_loc[i][2],err) 
	  ncel++; 
	  if (ncel > MaxNgrid)
	    {
	      fprintf(stderr,"*** Warning: Maximum number of locations (%d) reached\n",MaxNgrid);
	      fprintf(stderr,"             Stopping the grid-search\n");
	      flag=1;
	      break;
	    }	  
	}
      if(flag)
	break;
      Ngrid = 0;
      /* Extend grid-search */
      if (Nexp < MaxNexp)
	{	  
	  extend_xy_gs(old_loc,ncel,gfdep,ndep,opt_loc,opt,eq,it,Nexp,LL,UR,IZ,new_loc,&Ngrid) ;
	  if (Ngrid)
	    {
	      Nexp++;
	      continue;
	    }
	}
      Nexp = 0;
      it++;
    }
  fclose(tmp);
  opt->ref_flag = ref_flag ;
  /* Write output file */
  o_gs = openfile_wt(opt->xygsfile) ;
  *latopt = opt_loc[0][0] ;
  *lonopt = opt_loc[0][1] ;
  *depopt = opt_loc[0][2] ;
  *rmsopt = opt_loc[0][3]  ;
  fprintf( o_gs,"%10.4f %10.4f %10.4f %12.8f\n",*latopt,*lonopt,*depopt,*rmsopt);
  fprintf( o_gs,"%10.4f %10.4f %10.4f %12.8f\n",eq->evla,eq->evlo,eq->evdp,rmsini);
  tmp  = openfile_rt("_tmp_xy_table", &i) ;
  while ((i = getc(tmp)) != EOF)
    putc(i,o_gs) ;
  fclose( tmp) ;
  fclose(o_gs) ;
  /* Free Memory */
  for(i=0;i<opt->xy_Nopt;i++)
    free((void*)opt_loc[i]);
  free((void**)opt_loc);    
  for(i=0;i<MaxNgrid;i++)
    {
      free((void*)old_loc[i]);
      free((void*)new_loc[i]);
      free((void*)vrms[i]);
      free((void*)vms[i]);
    }
  free((void**)old_loc);
  free((void**)new_loc);  
  free((void**)vrms);
  free((void**)vms);
}




void 
set_data_vector(int nd,double *dv,double *tv,int *nsac,double ***data,char ***sacfiles,sachdr **hd_synt,
		str_quake_params *eq,structopt *opt,FILE *o_log)
{
  int i,j,ns,n1_data,n2_data,npts,ierror=1;
  double *tmparray    ;
  double gcarc, t0, Ptt, twp_beg, twp_end   ;
  char   datafile[FSIZE],dum[64],buf[LSIZE] ;
  FILE   *i_sac;  
  sachdr hd_data;

  /* Opening data file list */
  i_sac = openfile_rt(opt->i_saclst,nsac) ;
  /* Allocating memory */
  tmparray = double_alloc((int)__LEN_SIG__) ;  
  *data    = double_alloc2p(*nsac) ;
  hdr_alloc(&hd_data) ;
  opt->wgt  = double_alloc( *nsac ) ;
  *sacfiles = char_alloc2(*nsac,FSIZE) ;
  opt->rms_in = double_alloc(*nsac) ;
  opt->rms_r  = double_alloc(*nsac) ;
  opt->p2p    = double_alloc(*nsac) ;
  opt->avg    = double_alloc(*nsac) ;
  hdr_tab(hd_synt, *nsac) ; 
  /* Main loop (on channels) */
  ns = 0 ; /* Accepted channel counter */
  opt->dmin = 2.e4 ; 
  opt->dmax = 0.   ; 
  for(i=0; i<*nsac; i++)
    {
      if ( opt->op_pa <= 0 && opt->th_val <= 0) /* Read data file list */
	{ 
	  j = fscanf (i_sac,"%s",datafile) ;
	  fgets(buf,LSIZE,i_sac); /* end of line */
	  check_scan(1,j,opt->i_saclst,i_sac)  ;
	  opt->wgt[ns] = 1.0;
	}
      else 
	{
	  j = fscanf (i_sac, "%s %f %f %s %s %s %s %s %lf %lf %lf %lf %lf",
		      datafile,&(*hd_synt)[ns].az,&(*hd_synt)[ns].gcarc,dum,dum,dum,dum,dum,
		      &opt->rms_in[ns],&opt->rms_r[ns],&opt->p2p[ns],&opt->avg[ns],&opt->wgt[ns]) ;
	  strcpy(buf,opt->i_saclst);
	  strcat(buf," (nb of columns may be incorrect)");
	  check_scan(13, j, buf, i_sac) ;
	} 
      rhdrsac(datafile,&hd_data,&ierror) ; /* Read data header and weights */
      set_wgt(ns,&hd_data,opt)           ;
      if (opt->wgt[ns] <= 0.)
	{
	  fprintf(stderr,"**** null weight, rejected : %s\n", datafile) ;
	  fflush(stderr);
	  continue ;
	}
      /* Data Time Window */
      if ((float)(eq->pde_evdp) != hd_data.evdp)	
	{
	  fprintf(stderr,"WARNING : depth %f in sac header is different from the pde depth %f in CMTFILE\n",eq->pde_evdp, hd_data.evdp);
	  fprintf(stderr," ... you should carefully re-check gcarc and evdp header variables in file %s\n",datafile); 
	  fprintf(stderr," ... (gcarc from %s is used for the time windowing)\n",datafile); 
	  fflush(stderr);
	}
      gcarc = (double) hd_data.gcarc ; 
      if (opt->dmin > gcarc)
	opt->dmin = gcarc ;
      if (opt->dmax < gcarc)
	opt->dmax = gcarc ;
      trav_time(gcarc,tv,dv,nd,&Ptt, &ierror) ;
      wp_time_window(gcarc,eq->wp_win4,&twp_beg,&twp_end) ; 
      t0 = (double)hd_data.o ;
      if (TWPTT)
	t0 += Ptt ;
      else
	t0 += opt->ts;
      n1_data = (int)((t0 + twp_beg - (double)hd_data.b) / ((double)hd_data.delta)) ; /* first data Sample (corrected)  */
      n2_data = n1_data + (int)((twp_end - twp_beg) / ((double)hd_data.delta))      ; /* Last data Sample               */
      npts    = n2_data - n1_data + 1 ;
      if ( (n1_data<0) || (n2_data>=hd_data.npts) ) 
	{
	  fprintf(stderr,"**** Incomplete data, rejected : %s\n", datafile) ;
	  fprintf( o_log,"**** Incomplete data, rejected : %s\n", datafile) ;
	  fflush(stderr);
	  continue ;
	}
      /* Read data samples */
      hd_data.npts = n2_data + 1                     ;
      rdatsac(datafile, &hd_data, tmparray, &ierror) ;  
      (*data)[ns] = double_alloc(npts)               ;
      memcpy((*data)[ns],tmparray+n1_data,npts*sizeof(double)) ;
      /* Set data header */
      (*hd_synt)[ns].delta  = hd_data.delta  ;
      (*hd_synt)[ns].npts   = npts           ;
      (*hd_synt)[ns].nzyear = hd_data.nzyear ;
      (*hd_synt)[ns].nzjday = hd_data.nzjday ;
      (*hd_synt)[ns].nzhour = hd_data.nzhour ;
      (*hd_synt)[ns].nzmin  = hd_data.nzmin  ;
      (*hd_synt)[ns].nzsec  = hd_data.nzsec  ;
      (*hd_synt)[ns].nzmsec = hd_data.nzmsec ;
      (*hd_synt)[ns].o      = hd_data.o      ;
      (*hd_synt)[ns].b      = hd_data.b + hd_data.delta*(float)n1_data;
      (*hd_synt)[ns].user[0] = (float)(n1_data + 1) ;
      (*hd_synt)[ns].user[1] = (float)(n2_data + 1) ;
      (*hd_synt)[ns].user[2] = twp_beg ;
      (*hd_synt)[ns].user[3] = twp_end ;     
      (*hd_synt)[ns].user[4] = (*hd_synt)[ns].b - hd_data.b - (n1_data - 1)*hd_data.delta ;
      strcpy((*hd_synt)[ns].kstnm, hd_data.kstnm)   ; 
      strcpy((*hd_synt)[ns].knetwk, hd_data.knetwk) ; 
      strcpy((*hd_synt)[ns].kcmpnm, hd_data.kcmpnm) ; 
      (*hd_synt)[ns].az    = hd_data.az    ;
      (*hd_synt)[ns].gcarc = hd_data.gcarc ;
      (*hd_synt)[ns].stla  = hd_data.stla  ;
      (*hd_synt)[ns].stlo  = hd_data.stlo  ;
      (*hd_synt)[ns].evlo  = hd_data.evlo  ;
      (*hd_synt)[ns].evla  = hd_data.evla  ;
      (*hd_synt)[ns].evdp  = hd_data.evdp  ;
      /* Calculate seismogram peak-to-peak and average amplitude */
      if (opt->op_pa <= 0.) 
	{
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
}


int 
make_stat_list(sachdr *hd_synth,int nsac,char ***stats,float **stlats,float **stlons)
{
  int   i,j,k=0,nstat,flag;
  char  stat[9], conv[40];

  *stats  = char_calloc2(nsac,9) ;
  *stlats = float_calloc(nsac) ;
  *stlons = float_calloc(nsac) ;
  nstat = 0 ;
  for(i=0; i<nsac; i++)
    {
      j=sscanf(hd_synth[i].kstnm,"%s",stat);
      if (j!=1)
	{
	  fprintf(stderr,"ERROR: incorrect station or/and network name in sac header.\n");
	  exit(1);
	}
      flag = 0;
      for(j=0;j<nstat;j++)
	if (!strcmp(stat,(*stats)[j]))
	  {
	    flag = 1;
	    break;
	  }
      if (flag)
	continue;
      strcpy((*stats)[nstat],stat) ;
      /* Round lat/lon to 4 digits */
      sprintf(conv,"%12.4f %12.4f",hd_synth[i].stla,hd_synth[i].stlo);
      sscanf(conv,"%f %f",&(*stlats)[nstat],&(*stlons)[nstat]);
      nstat++;
    }
  for(i=nstat;i<k;i++)
    free((void*)(*stats)[i]);
  return nstat;
}


int
fill_kernel_G(sachdr *hd_GF, sachdr *hd_data, double Ptt, double twp_beg, 
	      double twp_end, double *elem_disp, double *G, structopt *opt, 
	      FILE *o_log)
{
  int i, npts   ;
  int n1_GF, n2_GF  ;
  double t0         ;
  double *g = &G[0] ;
  
  npts = hd_data->npts;
  /* Read Header    */
  if (hd_GF->delta != hd_data->delta) 
    {
      if (o_log!=NULL)
	fprintf( o_log,"**** Incorrect sampling period, rejected trace : %s\n",hd_data->kstnm) ;
      fprintf(stderr,"**** Incorrect sampling period, rejected trace : %s\n",hd_data->kstnm) ;
      return 1 ;
    }
  /* GF Time Window */
  t0 = (double)hd_GF->o ;
  if (TWPTT)
    t0 += Ptt ;
  else
    t0 += opt->ts;
  n1_GF = (int)((t0 + twp_beg - (double)hd_GF->b - opt->dts_val)  / ((double)hd_GF->delta)) ; /* first GF Sample (corrected) */
  n2_GF = n1_GF + (int)((twp_end - twp_beg) / ((double)hd_GF->delta)) ;                       /* Last GF Sample */
  if ( n2_GF >= hd_GF->npts ) /* GF Rejected */
    {
      if (o_log!=NULL)
	{
	  fprintf(o_log,"Stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", 
		  hd_data->kstnm, hd_data->knetwk, hd_data->kcmpnm, 
		  hd_data->gcarc, hd_data->az, hd_data->user[2], hd_data->user[3]) ;
	  fprintf( o_log,"**** Incomplete GF, rejected %s\n",hd_data->kstnm)       ;
	}
      fprintf(stderr,"**** Incomplete GF, rejected %s",hd_data->kstnm)             ;
      fprintf(stderr,"((n2_GF=%d)>=(npts=%d))\n", n2_GF, hd_GF->npts)              ;
      return 1 ;
    }
  else if ( n1_GF < 0 )       /* Fill negative samples with zeros */
    {
      if (o_log!=NULL)
	{
	  fprintf(o_log,"Stat: %-9s %-9s %-9s %8.1f %8.1f %8.1f %8.1f\n", 
		  hd_data->kstnm, hd_data->knetwk, hd_data->kcmpnm, 
		  hd_data->gcarc, hd_data->az, hd_data->user[2], hd_data->user[3])       ;
/* 	  fprintf( o_log,"**** Incomplete GF, filling with zeros : %s\n",hd_data->kstnm) ; */
	}
/*       fprintf(stderr,"**** Incomplete GF,")                                    ; */
/*       fprintf(stderr,"filling with zeros : %s (n1_GF=%d)\n",hd_data->kstnm,n1_GF )  ; */
      for(i=n1_GF; i<((n2_GF<0)?n2_GF+1:0); i++)
	G[i-n1_GF] = 0. ;
      g = &G[-n1_GF] ;
      npts += n1_GF  ;
      hd_GF->b += n1_GF*hd_GF->delta ;
      n1_GF = 0 ;     
    }
  /* fill g */
  if(n2_GF>=0)
    memcpy (g,elem_disp+n1_GF,npts*sizeof(double));  
  return 0;
}


/***************************************************/
/*           Set the kernel matrix G               */
/***************************************************/
void 
calc_kernel(str_quake_params *eq,structopt *opt,sachdr *hd_synth,int nsac,char *itype,
	    int nd,double *dv,double *tv, double ***G,FILE *o_log)
{
  int   i,j,k,l,ns,jstat,jsac,nstat,nsects,flag,flag2,ngfcomp=6,npts,ierror=1 ;
  long  int nerr ;
  char  **stats  ; 
  char  stacmp[]={'Z','L','T'}  ;
  float *stlats,*stlons,*dists,*azs,*bazs,*xdegs,b ;
  double **GFs,*Z,*TH,*PH,*E,*N,*x_conv,**WAV    ;
  double gcarc,Ptt,twp_beg,twp_end, *ref_vm=NULL ;
  double *b1,*b2,*a1,*a2,gain,dt=1.; 
  sachdr hdr; 
  nstat=make_stat_list(hd_synth,nsac,&stats,&stlats,&stlons);
  /* Memory Allocations */
  if (opt->ref_flag)
    {
      ref_vm = double_alloc(NM);
      for(i=0;i<NM;i++) 
	ref_vm[i] = eq->vm[1][i] ;
    }
  dists = float_calloc(nstat) ; 
  azs   = float_calloc(nstat) ;
  bazs  = float_calloc(nstat) ;
  xdegs = float_calloc(nstat) ;
  GFs   = double_alloc2(10,__LEN_SIG__) ;/* GFs: Rrr, Rtt, Rpp, Rrt  */
  Z     = double_alloc(__LEN_SIG__)     ;/*    Vertical components   */
  TH    = double_alloc(__LEN_SIG__)     ;/*    Radial components     */
  PH    = double_alloc(__LEN_SIG__)     ;/*    Transverse components */
  E     = double_alloc(__LEN_SIG__)     ;/*    East components       */
  N     = double_alloc(__LEN_SIG__)     ;/*    North components      */
  x_conv   = double_alloc(__LEN_SIG__)  ;
  WAV      = double_alloc2p(3) ;
  *WAV     = Z ;
  *(WAV+1) = N ;
  *(WAV+2) = E ;
  hdr_alloc(&hdr) ; /* SAC header allocation       */
  nsects = (eq->flow > 0.)? eq->filtorder : eq->filtorder/2 ;
  b1 = double_alloc(nsects) ; 
  b2 = double_alloc(nsects) ;
  a1 = double_alloc(nsects) ; 
  a2 = double_alloc(nsects) ;
  /* Distance, azimuth, back-azimuth, etc          */
  distaz(eq->evla, eq->evlo, stlats, stlons, nstat, dists, azs, bazs, xdegs, &nerr) ;
  flag = 0 ;  
  for(i=0;i<ngfcomp;i++) /* Main loop */
    {
      ns   = 0 ; /* Channel counter */
      for(j=0;j<ngfcomp;j++) /* Inititializing the MT components */
	eq->vm[1][j] = 0. ;
      eq->vm[1][i]   = 1. ;
      for (jstat=0;jstat<nstat;jstat++) /* Computing exitation kernels for MT component #i at each station */
	{ 
	  gcarc = (double)hd_synth[ns].gcarc     ;
	  trav_time(gcarc,tv,dv,nd,&Ptt,&ierror) ;
	  wp_time_window(gcarc,eq->wp_win4,&twp_beg,&twp_end) ; /* Data Time Window  */
	  flag2 = 0;
	  /* Computing Z, TH, PH  */ 
	  fast_synth_sub(azs[jstat],bazs[jstat],xdegs[jstat],tv,dv,nd,eq,&hdr,GFs,Z,TH,PH);
	  rotate_traces(TH, PH, bazs[jstat], hdr.npts, N, E)  ; /* Rotating TH, PH to N, E */
	  b    = hdr.b    ;
	  npts = hdr.npts ;
	  k    = strlen(stats[jstat]);
	  for(jsac=0;jsac<nsac;jsac++)
	    {
	      l = nbchar(hd_synth[jsac].kstnm);
	      if ( k==l && !strncmp(hd_synth[jsac].kstnm,stats[jstat],k) )
	         {
	          for(j=0;j<3;j++) /* Which component (either LHZ, LHL or LHT)? */
		    if (hd_synth[jsac].kcmpnm[2] == stacmp[j])
		      break;
	          if (j==3)
		    {
		      fprintf(stderr,"*** ERROR: Unknownk component %s for sta %s\n",hd_synth[jsac].kcmpnm,hd_synth[jsac].kstnm);
		      fprintf(stderr,"    -> Exiting\n") ;
		      fflush(stderr);
		      exit(1);
    		    }
	          hdr.b    = b    ; /* Re-initialize hdr */
	          hdr.npts = npts ;
	          conv_by_stf(eq->ts,eq->hd,itype,&hdr,WAV[j],x_conv) ;/* Perform convolution */
	          if (flag == 0) /* Set the butterworth sos (dt must be the same for all stations)   */
	    	    {
		      flag = 1 ; 
		      dt = (double)hdr.delta;
		      if (eq->flow>0.)
			bpbu2sos(eq->flow,eq->fhigh,dt,eq->filtorder,&gain,b1,b2,a1,a2);
		      else
			lpbu2sos(eq->fhigh,dt,eq->filtorder,&gain,b1,b2,a1,a2);		  		      
	    	    }
	          else if (dt != (double)hdr.delta) /* Check the sampling frequency (must be uniform) */
		    {
		      fprintf(stderr, "*** ERROR: non uniform samp. period between sac files (%s)\n",hd_synth[jsac].kstnm);
		      fprintf(stderr,"    -> Exiting\n") ;
		      fflush(stderr);
		      exit(1);
		    }	  
		  filter_with_sos(gain,b1,b2,a1,a2,nsects,x_conv,hdr.npts) ; /* Apply sos */
	          flag2 = fill_kernel_G(&hdr,&(hd_synth[jsac]),Ptt,twp_beg,twp_end,x_conv,G[jsac][i],opt,o_log);
	          if (flag2)
		    {
		      fprintf(stderr,"*** ERROR: Incomplete green function for %s, %s\n",stats[jstat],hd_synth[jsac].kstnm) ;
		      fprintf(stderr,"    -> Exiting\n") ;
		      fflush(stderr);
		      exit(1);
		    }
	          hd_synth[jsac].b  = hdr.b + hd_synth[jsac].o - hdr.o ;
	          ns++;
	        }
          }
	}
      if (nsac!=ns)
	{	
	  fprintf(stderr,"\n*** ERROR: Kernel G is incomplete (%d vs %d)\n",nsac,ns);  
	  fprintf(stderr,"    -> Exiting\n") ;
	  fflush(stderr);
	  exit(1);	  
	}
    }

  /* Free Memory */
  if (opt->ref_flag)
    {
      for(i=0;i<NM;i++) 
	eq->vm[1][i] = ref_vm[i] ;
      free((void*)ref_vm);
    }
  for(i=0;i<nstat;i++)
    free((void*)stats[i]);
  free((void**)stats) ;
  free((void*)stlats) ;
  free((void*)stlons) ;
  free((void*)dists)  ;
  free((void*)xdegs)  ;
  free((void*)bazs)   ;
  free((void*)azs)    ;
  for(i=0;i<10;i++)
    free((void*)GFs[i]) ;
  free((void**)GFs)     ;
  free((void*)x_conv)   ;
  free((void*)Z)    ;
  free((void*)N)    ;
  free((void*)E)    ;
  free((void*)TH)   ;
  free((void*)PH)   ;
  free((void**)WAV) ;
  free((void*)a1)   ;
  free((void*)a2)   ;
  free((void*)b1)   ;
  free((void*)b2)   ;
}

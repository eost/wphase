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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "proto_alloc.h"  /* proto_alloc.c  */
#include "read_i_files.h" /* read_i_files.c */

#ifndef RX
#define	RX 17
#endif

#ifndef RY
#define	RY 9
#endif


#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.)
#endif


double **set_mt(double *vm) ;
void   get_planes(double *vm, double **eval3, double ***evec3, double *s1, double *d1, 
		  double *r1, double *s2,double *d2,double *r2) ;
void   jacobi(double **a,int n, int np, double *d, double **v, int *nrot) ;
void   eigsrt(double *d, double **v, int n) ;
int    charplot(double *M, double s1, double d1, double s2, double d2, 
		char D, char P, char W, char B, char sep, char pnod, 
		int rx, int ry, FILE *stream) ;
void   format_latlon(double lat, double lon, char *slat, char *slon) ;
void   vn2sdr(double *vn, double *vs, double *s, double *d, double *r);
double round(double x);

int 
main(int argc, char *argv[])
{
  int    i, p;
  double s1, d1, r1, s2, d2, r2, plg[3], azm[3] ;
  double *eval3, **evec3, M0, Mw, scale         ;
  char pdela[8], pdelo[9], cenla[8],cenlo[9]    ;
  str_quake_params eq ;

  if (argc < 2)
    {
      fprintf(stderr,"Error syntax : %s CMTFILE\n",argv[0]);
      exit(1) ;
    }
  
  /* Allocate memory             */
  eq.vm    = double_alloc2p(2) ;
  eq.vm[1] = double_calloc(6)  ;
  
  /* Read cmtfile                */
  strcpy(eq.cmtfile,argv[1])  ;
  get_cmtf(&eq, 2) ;
  
  /* Best double couple solution */
  get_planes(eq.vm[1], &eval3, &evec3, &s1,&d1,&r1, &s2,&d2,&r2) ;
  M0     = ((fabs(eval3[0]) + fabs(eval3[2])) * (double)POW) / 2. ; 
  Mw     = (log10(M0) - 16.1) / 1.5 ;  
  
  /* Principal axes              */
  for (i=0; i<3 ; i++)
    {
      azm[i] = atan2(evec3[2][i],-evec3[1][i])/DEG2RAD+360.;
      scale  = evec3[1][i]*evec3[1][i] + evec3[2][i]*evec3[2][i] ;
      plg[i] = atan2(-evec3[0][i],sqrt(scale))/DEG2RAD;
      if(plg[i]<0.0)
		{
		  plg[i] *= -1. ;
		  azm[i] += 180.;
		}
	  azm[i] = fmod(azm[i], 360.);
	}
  
  /* Scale Seismic moment        */
  scale = 0. ;
  for (i=0;i<6;i++)
    if (fabs(eq.vm[1][i]) > scale)
      scale = fabs(eq.vm[1][i]) ;
  p = (int)log10(scale) ;
  scale = pow(10,(double)p)  ;
  for (i=0;i<6;i++)
    eq.vm[1][i] /= scale;
  for (i=0;i<3;i++)
    eval3[i] /= scale;
  p += (int)round(log10(POW));
  /* Display                      */
  format_latlon(eq.pde_evla, eq.pde_evlo, pdela, pdelo);
  format_latlon(eq.evla, eq.evlo, cenla, cenlo);
  printf("Moment mag.  : %5.2f\n",Mw)  ;  
  printf("PDE location : Lat=%7s; Lon=%8s; Dep=%5.1f km\n",pdela,pdelo,eq.pde_evdp) ; 
  printf("Centroid loc.: Lat=%7s; Lon=%8s; Dep=%5.1f km\n",cenla,cenlo,eq.evdp)        ;
  printf("Origin time  : %04d/%02d/%02d %02d:%02d:%05.2f\n", eq.ot_ye,eq.ot_mo,
	 eq.ot_dm,eq.ot_ho,eq.ot_mi, (double)eq.ot_se+(double)eq.ot_ms/1000.)  ;
  printf("Time delay   : %-5.1f sec\n",eq.ts) ;
  printf("Half duration: %-5.1f sec\n\n",eq.hd) ;
  printf("Moment tensor: scale= 1.0E+%02d dyn.cm\n",p) ;
  printf("rr=%6.3f ; tt=%6.3f ; pp=%6.3f \n",eq.vm[1][0],eq.vm[1][1],eq.vm[1][2]) ;
  printf("rt=%6.3f ; rp=%6.3f ; tp=%6.3f \n\n",eq.vm[1][3],eq.vm[1][4],eq.vm[1][5]) ;
  printf("Principal Axes: \n") ;
  printf("1.(T) Val=%6.3f ; Plg= %3.0f ; Azm=%3.0f\n",eval3[0],plg[0],azm[0]);
  printf("2.(N)     %6.3f ;      %3.0f ;     %3.0f\n",eval3[1],plg[1],azm[1]);
  printf("3.(P)     %6.3f ;      %3.0f ;     %3.0f\n\n",eval3[2],plg[2],azm[2]);
  printf("Best Double Couple: M0=%8.2E dyn.cm\n",M0) ;
  printf("NP1: Strike=%5.1f ; Dip=%4.1f ; Slip=%6.1f \n",s1,d1,r1) ;
  printf("NP2: Strike=%5.1f ; Dip=%4.1f ; Slip=%6.1f \n\n",s2,d2,r2) ;
  /* printf("NP1: Strike=%3.0f ; Dip=%.4f ; Slip=%3.0f \n",s1,d1,r1)   ; */
  /* printf("NP2: Strike=%3.0f ; Dip=%.4f ; Slip=%3.0f \n\n",s2,d2,r2) ; */
  charplot(eq.vm[1], s1,d1, s2,d2, '-', '#', ' ', '\0','\0','\0', RX, RY, stdout)  ;
  
  free((void*)eq.vm[1]) ; 
  free((void**)eq.vm) ; 
  for(i=0 ; i<3 ; i++) 
    free((void*)evec3[i]) ;
  free((void**)evec3)     ;
  free((void*) eval3)     ;

  return 0;
}


void format_latlon(double lat, double lon, char *slat, char *slon)
{
  if (lat < 0.)
    sprintf(slat,"%6.2fS",-lat) ;
  else
    sprintf(slat,"%6.2fN", lat) ;

  if (lon < 0.)
    sprintf(slon,"%7.2fW",-lon) ;
  else
    sprintf(slon,"%7.2fE", lon) ;
}

double **
set_mt(double *vm)
{
  double **TM ;
  TM       = double_alloc2(3,3) ;
  TM[0][0] =   vm[0] ;
  TM[1][1] =   vm[1] ;
  TM[2][2] =   vm[2] ; 
  TM[0][1] =   vm[3] ;
  TM[0][2] =   vm[4] ;
  TM[1][2] =   vm[5] ;
  TM[1][0] = TM[0][1] ;
  TM[2][0] = TM[0][2] ;
  TM[2][1] = TM[1][2] ;

  return TM ;
}

void 
get_planes(vm, eval3, evec3, s1,d1,r1, s2,d2,r2)
     double *vm, **eval3, ***evec3, *s1, *d1, *r1, *s2, *d2, *r2;
{
  int    nrot, i ;
  double **TM;
  double *vn1, *vn2 ;
  
  /* Memory allocation */
  *eval3 = double_alloc(3)    ;
  *evec3 = double_alloc2(3,3) ;
  vn1    = double_alloc(3)    ;
  vn2    = double_alloc(3)    ;
  
  /* Tensor representation */
  TM = set_mt(vm) ;

  /* Get eigvalues and eigvectors*/
  jacobi(TM,3,3,(*eval3),*evec3,&nrot) ;
  eigsrt((*eval3),*evec3,3) ;

  for(i=0 ; i<3 ; i++)
  {
     vn1[i] = ((*evec3)[i][0]+(*evec3)[i][2])/sqrt(2.) ;
     vn2[i] = ((*evec3)[i][0]-(*evec3)[i][2])/sqrt(2.) ;
  }
  vn2sdr(vn1, vn2, s1, d1, r1); 
  vn2sdr(vn2, vn1, s2, d2, r2); 
  
  /* Memory Freeing */
  free((void*)vn1) ;
  free((void*)vn2) ;
  for(i=0 ; i<3 ; i++) 
    free((void*)TM[i]) ;
  free((void**)TM)    ;
}

void vn2sdr(double *vn, double *vs, double *s, double *d, double *r)
{
   const float EPSI = 0.001;
   int   i;
   /* printf("%f %f %f\n", vn[0], vn[1], vn[2]); */
   if (vn[0] < 0.)              // Upwards normal
      for(i=0; i<3; i++)
      {
         vn[i] *= -1.;
         vs[i] *= -1.;
      }

   if ( vn[0] > 1. - EPSI )     // Horizontal plane
   {
      *s = 0.;
      *d = 0.;
      *r = atan2(-vs[2], -vs[1]);
   }

   else if ( vn[0] < EPSI )    // Vertical plane
   {
	  *s = atan2(vn[1], vn[2]);
      *d = M_PI/2.;
      *r = atan2(vs[0], -vs[1]*vn[2] + vs[2]*vn[1]);
   }

   else                        // Oblique plane
   { 
      *s = atan2(vn[1], vn[2]);
      *d = acos(vn[0]);
      *r = atan2((-vs[1]*vn[1] - vs[2]*vn[2]), (-vs[1]*vn[2] + vs[2]*vn[1])*vn[0]);
   }

   *s /= (double)DEG2RAD;
   if ((*s) < 0.) (*s) += 360.;
   *d /= (double)DEG2RAD;
   *r /= (double)DEG2RAD;
   return;
}

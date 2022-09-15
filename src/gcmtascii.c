/***************************************************************************
 * 
 *                      W phase source inversion package              
 *                                -------------
 * 
 *         Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
 *                       gcmtascii was provided by Anthony Jamelot
 *
 *  (c) California Institute of Technology and Universite de Strasbourg / CNRS 
 *                                   April 2013
 * 
 *     Neither the name of the California Institute of Technology (Caltech) 
 *     nor the names of its contributors may be used to endorse or promote 
 *     products derived from this software without specific prior written 
 *     permission
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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

/*
 * 

CENTROID-MOMENT-TENSOR  SOLUTION
WCMT EVENT:     M201103110546A  
DATA: II DK IU CU G  GE MN LD     
MANTLE WAVES:  100S, 271C, T=300  
TIMESTAMP:      Q-20110311234057  
CENTROID LOCATION:
ORIGIN TIME:      05:47:32.8 0.2 
LAT:37.52N 0.01;LON:143.05E 0.02  
DEP: 20.0  FIX;TRIANG HDUR: 70.0
MOMENT TENSOR: SCALE 10**29 D-CM
RR= 1.730 0.006; TT=-0.281 0.005
PP=-1.450 0.005; RT= 2.120 0.068
RP= 4.550 0.065; TP=-0.657 0.004
PRINCIPAL AXES:
1.(T) VAL=  5.305;PLG=55;AZM=295
2.(N)       0.014;     0;    205
3.(P)      -5.319;    35;    115
BEST DBLE.COUPLE:M0= 5.31*10**29
NP1: STRIKE=203;DIP=10;SLIP=  88
NP2: STRIKE= 25;DIP=80;SLIP=  90

            ---#######-           
        ---#############---       
      --################-----     
    --#################--------   
   --##################---------  
  --##################----------- 
  -#######   #########----------- 
 --####### T ########-------------
 -########   #######--------------
 -#################---------------
 -################--------   -----
  -##############--------- P ---- 
  -#############----------   ---- 
   -###########-----------------  
    -#########-----------------   
      -######----------------     
        -###---------------       
            -----------           
Example;
cc -O2 -Wall -DCENTER=CPPT -D__LEN_SIG__=5000 -D__FSIZE__=128 -D__IDSIZE__=16 -D__LSIZE__=256 -fopenmp -D__GFS_01D__=1 -c gcmtascii.c
cc -O2 -fsecond-underscore gcmtascii.o proto_alloc.o rwtextfiles.o jacobi.o charplot.o read_i_files.o -o ../bin/gcmtascii -lm

*/
double **set_mt(double *vm) ;
void   get_planes(double *vm, double **eval3, double ***evec3, double *s1, double *d1, 
		  double *r1, double *s2,double *d2,double *r2) ;
void   jacobi(double **a,int n, int np, double *d, double **v, int *nrot) ;
void   eigsrt(double *d, double **v, int n) ;
int    charplot(double *M, double s1, double d1, double s2, double d2, 
		char D, char P, char W, char B, char sep, char pnod, 
		int rx, int ry, FILE *stream) ;
void   format_latlon(double lat, double lon, char *slat, char *slon, int res) ;
void   vn2sdr(double *vn, double *vs, double *s, double *d, double *r);

int main(int argc, char *argv[])
{
  int    buildtime, i, p;
  int stat_count=0;
  int phase_count=0;
  double s1, d1, r1, s2, d2, r2, plg[3], azm[3] ;
  double *eval3, **evec3, M0, Mw, scale         ;
  char pdela[8], pdelo[9], cenla[8],cenlo[9]    ;
  char tRegion[60] ;
  str_quake_params eq ; //params from cmt_file
  str_quake_params im ; //params from i_master
  // time variables
  time_t now = time(NULL);
  struct tm *t ;  
  struct tm evt_time;  
  char buffer_time [60],buffer_data [60];
  // Get Stations list variables
  char Network[][3] = { "" , "" , "" , "" , "", "" , "", "", ""};
  char tmp[3], str[80];
  char tmp2[6],stat[6];
  FILE *f;
  char AgencyID[5];
  int count,k,n,nlength;
  char **keys ;
  
  if (argc < 4)
    {
      fprintf(stderr,"Error syntax, should be : %s PATH/xy_WCMTSOLUTION PATH/xy_o_wpinversion AgencyID Minutes_after_OT\n",argv[0]);
      exit(1) ;
    }
  if (argc == 5)
    {
      sscanf(argv[4],"%d",&buildtime);
    }
  else
    {
      buildtime=30;
    }
  strncpy(AgencyID,argv[3],4);
  AgencyID[4]='\0';
  
  /* Allocate memory             */
  eq.vm    = double_alloc2p(2) ;
  eq.vm[1] = double_calloc(6)  ;
  
  /* Read cmtfile                */
  strcpy(eq.cmtfile,argv[1])  ;
  get_cmtf(&eq, 2) ;
  
  /* Read i_master file to get filt_cf2 */
  keys = char_alloc2(1, 8) ;
  strcpy(keys[0],"filt_cf2")   ;
  get_i_master("i_master", keys, 1, &im) ;
  free((void*)keys[0]) ;
  free((void**) keys ) ;
  /* Create time oject */
  //  -- created time -- 
  t=gmtime(&now);
  //  -- seism time --
  evt_time.tm_year = eq.ot_ye-1900;
  evt_time.tm_mon = eq.ot_mo-1; // should test over other OS
  evt_time.tm_mday = eq.ot_dm;
  evt_time.tm_hour = eq.ot_ho;
  evt_time.tm_min = eq.ot_mi;
  evt_time.tm_sec = eq.ot_se;
  evt_time.tm_isdst = 0;
  // GCMT LIKE
  // strftime (buffer_time,60,"%B %d, %Y",&evt_time);
  // CPPT LIKE
  strftime (buffer_time,60,"%F %H:%M:%S UTC",&evt_time);
	
  /* Best double couple solution */
  // Moment tensor values are divide by POW (defined at 1E+28) when read from cmtfile
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


  /* Get Network of stations used for moment tensor solution    */
  if (( f=fopen(argv[2],"r"))==NULL)
    {
      fprintf (stderr, "ERROR (read) : opening file: %s \n", argv[2]) ;
      exit (1) ;
    }
  i=0;
  n=0;
  
  nlength=(int)(sizeof(Network) / 3);

  while ( fscanf(f,"%s %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f\n",str) != EOF ) {
    sscanf(str+28,"%[a-zA-Z0-9].%[a-zA-Z0-9].",tmp,tmp2);
    if ( i == 0 ) { 
      strncat(Network[0],tmp,3);
      strcpy(stat,tmp2);
      stat_count++;
      n++;
    }
    else if (strcmp(stat,tmp2) != 0) {
      stat_count++ ;
      strcpy(stat,tmp2);
    }
    if ( i != 0 && n < nlength) {
      count=0;
      for ( k=0 ; k<nlength ; k++ ) 
	if (strncmp(Network[k],tmp,3) == 0)  count++;
      
      if ( count == 0 ) { 
	strncat(Network[n],tmp,3);
	n++;
      }
    }
    i++;
  }	
  phase_count = i;
   
  fclose(f); 
  sprintf(buffer_data,"DATA:");
  
  for ( k=0 ; k<nlength ; k++ ) {
    strcat(buffer_data," ");
    strcat(buffer_data,Network[k]);
  }
  
  /* Display as GCMT bulletin format     */
  format_latlon(eq.pde_evla, eq.pde_evlo, pdela, pdelo,3);
  format_latlon(eq.evla, eq.evlo, cenla, cenlo,3);
  sscanf(&eq.pdeline[5],"%*d %*d %*d %*d %*d %*f %*f %*f %*f %*f    %59[a-zA-Z ,.]s",tRegion);
  // GCMT LIKE
  // printf("%s %s, Mw=%4.2f\n\n",buffer_time,tRegion,Mw) ;                       
  // CPPT LIKE
  printf("Mw=%4.2f, %s\n",Mw,tRegion) ;
  printf("Event Time: %s\n",buffer_time) ;
  printf("Epicenter : %s %s %4.1fkm\n",pdela,pdelo,eq.pde_evdp);
  printf("Centroid  : %s %s %4.1fkm\n\n",cenla,cenlo,eq.evdp);
  /* Display as GCMT bulletin format     */
  format_latlon(eq.pde_evla, eq.pde_evlo, pdela, pdelo,2);
  format_latlon(eq.evla, eq.evlo, cenla, cenlo,2);
  printf("CENTROID-MOMENT-TENSOR  SOLUTION\n") ;
  printf("AUTO WPHASE INVERSION - 0T+%dMIN \n",buildtime) ;
  printf("%-4s EVENT:     C%04d%02d%02d%02d%02dA\n",AgencyID,eq.ot_ye,eq.ot_mo,
	 eq.ot_dm,eq.ot_ho,eq.ot_mi);
  // END CPPT LIKE	 
  printf("%s\n",buffer_data);
  printf("TIMESTAMP:      Q-%04d%02d%02d%02d%02d%02d\n",t->tm_year+1900, t->tm_mon+1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec) ;
  printf("WPHASE WAVES: %4dS,%4dC, T=%3.0f\n",stat_count,phase_count,floor(1./im.fhigh)) ; 
  printf("CENTROID LOCATION:\n") ;
  printf("ORIGIN TIME:      %02d:%02d:%04.1f 0.0\n",eq.ot_ho,eq.ot_mi, (double)eq.ot_se+(double)eq.ot_ms/1000.) ; 
  printf("LAT:%s 0.00;LON:%s 0.00\n",cenla,cenlo) ;  
  printf("DEP:%5.1f  0.0;TRIANG HDUR:%5.1f\n",eq.evdp,eq.hd) ;  
  printf("MOMENT TENSOR: SCALE 10**%02d D-CM\n",(int)log10(POW)+p) ;
  printf("RR=%6.3f 0.000; TT=%6.3f 0.000\n",eq.vm[1][0],eq.vm[1][1]) ;
  printf("PP=%6.3f 0.000; RT=%6.3f 0.000\n",eq.vm[1][2],eq.vm[1][3]) ;
  printf("RP=%6.3f 0.000; TP=%6.3f 0.000\n",eq.vm[1][4],eq.vm[1][5]) ;
  printf("PRINCIPAL AXES:\n") ;
  printf("1.(T) VAL=%7.3f;PLG=%2.0f;AZM=%3.0f\n",eval3[0],plg[0],azm[0]) ;
  printf("2.(N)     %7.3f;    %2.0f;    %3.0f\n",eval3[1],plg[1],azm[1]) ;
  printf("3.(P)     %7.3f;    %2.0f;    %3.0f\n",eval3[2],plg[2],azm[2]) ;
  printf("BEST DBLE.COUPLE:M0= %4.2f*10**%d \n",M0 / (double)pow(10,(int)log10( M0 )),(int)log10(M0));
  printf("NP1: STRIKE=%3.0f;DIP=%2.0f;SLIP=%4.0f\n",s1,d1,r1) ;
  printf("NP2: STRIKE=%3.0f;DIP=%2.0f;SLIP=%4.0f\n\n",s2,d2,r2) ;
  charplot(eq.vm[1], s1,d1, s2,d2, '-', '#', ' ', '\0','\0','\0', 16, 9, stdout)  ;
  
  //free((void*)eq.vm[0]) ; ERROR NOT NECESSARY, THIS ARRAY IS NOT INITIALIZED
  free((void*)eq.vm[1]) ; 
  free((void**)eq.vm) ; 
  for(i=0 ; i<3 ; i++) 
    free((void*)evec3[i]) ;
  free((void**)evec3)     ;
  free((void*) eval3)     ;
  
  
  
   
  return 0;
}


void format_latlon(double lat, double lon, char *slat, char *slon, int res)
{
  if (res == 3) {
    if (lat < 0.)
      sprintf(slat,"%6.3fS",-lat) ;
    else
      sprintf(slat,"%6.3fN", lat) ;

    if (lon < 0.)
      sprintf(slon,"%7.3fW",-lon) ;
    else
      sprintf(slon,"%7.3fE", lon) ;
  }
  else if ( res == 2 ) {
    if (lat < 0.)
      sprintf(slat,"%5.2fS",-lat) ;
    else
      sprintf(slat,"%5.2fN", lat) ;

    if (lon < 0.)
      sprintf(slon,"%6.2fW",-lon) ;
    else
      sprintf(slon,"%6.2fE", lon) ;
  }
}

double **set_mt(double *vm)
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


void get_planes(vm, eval3, evec3, s1,d1,r1, s2,d2,r2)
    double *vm, **eval3, ***evec3, *s1, *d1, *r1, *s2, *d2, *r2;
{
    int    nrot, i ;
    double **TM;
    double *vn1, *vn2 ;
    double tmp;
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
  
    if (*d1 > *d2)
    {
         tmp = *s1; *s1 = *s2; *s2 = tmp;
         tmp = *d1; *d1 = *d2; *d2 = tmp;
         tmp = *r1; *r1 = *r2; *r2 = tmp;
    }
  
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
    if ((*s) < 0.) 
        (*s) += 360.;
    *d /= (double)DEG2RAD;
    *r /= (double)DEG2RAD;
    return;
}

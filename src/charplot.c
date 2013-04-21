/****************************************************************
*	W phase package - ASCII beach balls
*                                           
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proto_alloc.h"

int 
pnodal(char **figure, double s, double d, int rx, int ry, char pnod)
{
  int    th, jx, jy;
  double x, y, z, thr, azi, ain, r;
  for (th=0; th<=360; th++)
    {
      thr = ((double)th) * M_PI/360. ;
      y  = cos(s) * cos (thr) - sin(s) * sin(thr) * cos(d) ;
      x  = sin(s) * cos (thr) + cos(s) * sin(thr) * cos(d) ;
      z  = sin(thr) * sin(d)      ;
      azi = atan2(y,x)            ;
      ain = atan2(sqrt(1-z*z),z)  ;
      r   = sqrt(2) * sin(ain/2.) ;
      jx = (int) (r * cos(azi) * (double)rx) ;
      jy = (int) (r * sin(azi) * (double)ry) ;
      figure[ry-jy][rx+jx] = pnod ;
    }
  return 1;
}

int 
charplot(double *M, double s1, double d1, double s2, double d2, 
	     char D, char C, char W, char B, char sep, char pnod, int rx, 
	     int ry, FILE *stream)
{
  int 	jx, jy, rint       ;
  int 	jxP, jxT, jyP, jyT ;
  double  radius, aP, aT   ;
  double  ain, azi         ;
  double  x, y, z, amp     ;
  char    **figure         ;
  
  figure = char_calloc2(2*ry+1,2*rx+1) ;
  jxP = 0; jyP=0; jxT=0; jyT=0;
  aP = M[0]; aT = M[0]; 	
  for(jy = ry; jy >= -ry; jy--)
    {
      for(jx = -rx; jx <= rx; jx++)
	{
	  radius = hypot((double)jx, (double)jy*(double)rx/(double)ry);
	  rint   = (int) (radius+0.5) ;
	  figure[ry+jy][rx+jx] = W;
	  if (rint <= rx)
	    {
	      figure[ry+jy][rx+jx] = D;
	      ain = 2.*asin(radius/(double)rx/M_SQRT2);
	      azi= atan2(jx, jy*(double)rx/(double)ry);
	      x =  sin(ain)*sin(azi);
	      y =  sin(ain)*cos(azi);
	      z = -cos(ain);
	      amp  =      M[0]*z*z + M[1]*y*y + M[2]*x*x;
	      amp += 2.*( M[3]*z*y + M[4]*z*x + M[5]*x*y);
	      if( amp > 0. ) 
		figure[ry+jy][rx+jx] = C;
	      if( amp < aP ) {
		aP = amp; 
		jxP = jx; 
		jyP = jy; }
	      if( amp > aT ) {
		aT = amp; 
		jxT = jx; 
		jyT = jy; }
	    }
	  else if (B != '\0' && radius >= rx && radius < rx+1) 
	    figure[ry+jy][rx+jx] = B;

	}
    }
  
  if (sep == '\0')  /* T and P axis */
    {
      for(jx = jxP-1; jx <= jxP+1; jx++)
	for(jy = jyP-1; jy <= jyP+1; jy++)
	  if(jy >= -ry  &&  jy <= ry && jx >= -rx  &&  jx <= rx) 
	    figure[ry+jy][rx+jx] = W;
      for(jx = jxT-1; jx <= jxT+1; jx++)
	for(jy = jyT-1; jy <= jyT+1; jy++)
	  if(jy >= -ry  &&  jy <= ry && jx >= -rx  &&  jx <= rx) 
	    figure[ry+jy][rx+jx] = W;
      figure[ry+jyT][rx+jxT] = 'T';
      figure[ry+jyP][rx+jxP] = 'P';
      /* 	Printing    */
      for(jy = -ry; jy <= ry; jy++)
	{
	  for(jx = -rx; jx <= rx; jx++)
	    fprintf(stream,"%c", figure[ry+jy][rx+jx]);
	  fprintf(stream,"\n");
	}
    }
  else            /* Nodal planes */
    {  
      pnodal(figure, s1*M_PI/180. , d1*M_PI/180. , rx, ry, pnod) ;
      pnodal(figure, s2*M_PI/180. , d2*M_PI/180. , rx, ry, pnod) ; 
      /* 	Printing    */
      for(jy = -ry; jy <= ry; jy++)
	{
	  for(jx = -rx; jx <= rx; jx++)
	    fprintf(stream,"%c%c", figure[ry+jy][rx+jx],sep);
	  fprintf(stream,"\n");
	}
    }
  for(jx=0;jx<2*ry+1;jx++)
    free((void*)figure[jx]);
  free((void**)figure);
  return(1);
}


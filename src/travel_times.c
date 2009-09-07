#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "proto_alloc.h"    /* proto_alloc.c    */
#include "rwtextfiles.h" /* rwtextfiles.c */
#include "travel_times.h" /* rwtextfiles.c */


void 
make_table(char* name, char* table)
{
  int  nc ;
  char local_name[256];
  strcpy(local_name, name);
  if(getenv(local_name))
    strncpy(table, getenv(local_name), 256);
  else
    {
      fprintf(stderr, "Error: %s environment variable not defined\n", name);
      exit(1);
    }
  nc = strlen(table);
  if (table[nc-1] == '/')
    table[nc-1] = '\0' ;
  strcat(table, "/P_tt_JB") ;
}


void 
trav_time_init(nh, nd, h, dv, tv, ierror)
int    *nh, *nd ;
double *h, *dv, *tv ;
int    *ierror ;
{
  int    i, j, nl, flag, co;
  double *depths, *tmp, h1, h2, t1, t2;
  char   *file;
  FILE   *P_ttfile ;
 
  tmp    = double_alloc((*nh+1) * (*nd+1));
  depths = double_alloc(*nh);
  file  = char_alloc(256) ;

  make_table("WPHASE_HOME", file) ;
  P_ttfile = openfile_rt(file,&nl);
  
  /* Read file */
  i=0;
  while((flag = fscanf(P_ttfile, "%lf", &tmp[i++]))!=EOF)
	check_scan(1, flag, file, P_ttfile);
  fclose(P_ttfile); 

  /* Set clother depths */
  for (i=1; i<*nh; i++) 
    {
      if ((*h >= tmp[i]) && (*h <= tmp[i+1])) 
	{
	  flag = 0;
	  break;
	}
    }
  if (flag) 
    {
      if (*ierror == 1)
	{
	  fprintf(stderr, "Error: depth: %f out of range\n", *h);
	  exit(1);
	}
      fprintf(stderr, "WARNING : depth: %f out of range; extrapolating ... \n", *h); 
      *ierror = 1;
      i--;
    }

  j  = i;
  h1 = tmp[j];
  h2 = tmp[j+1];
  co = 0;
  /* Set distance and travel time arrays */
  for (i = 0; i < *nd; i++)
    {
      co   += (*nh)+1;
      /* Set epicentral distance */
      dv[i] = tmp[co];
      /* Set surrounding travel-times for this epicentral distance */
      t1  = tmp[co+j] ;
      t2  = tmp[co+j+1] ;
      /* Linear Interpolation */
      tv[i] =  t1 + (t2 - t1) * (*h-h1)/(h2-h1) ;
    }
  free((void*)depths);
  free((void*)tmp);
}
 



void 
trav_time(d, tv, dv, nd, P_tt, ierror)
double *d, *tv, *dv;
int    *nd;
double *P_tt; 
int    *ierror;
{
  int    i, flag=1;

  /* Set clother distances */
  for (i=0; i<*nd; i++) 
    {
      if ((*d >= dv[i]) && (*d <= dv[i+1])) 
	{
	  flag = 0;
	  break;
	}
    }
  if (flag) 
    {
      if (*ierror == 1)
	{
	  fprintf(stderr, "Error: distance: %f out of range\n", *d);
	  exit(1);
	}
      fprintf(stderr, "WARNING : distance: %f out of range; extrapolating ... \n", *d); 
      i--;
    }
  *P_tt = tv[i] + (tv[i+1] - tv[i]) * (*d - dv[i])/(dv[i+1] - dv[i]) ;
} 


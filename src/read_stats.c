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
#include "proto_alloc.h"
#include "rwtextfiles.h"

/* Read the filelist of station and networks */
int 
read_stats(char* file, char*** stats, char*** nets, float** stlats, float** stlons)
{
  FILE* fp;  
  int   ns, jstat, nstat;
  char* line;

  line = char_alloc(64);
  fp = openfile_rt(file,&nstat);
  *stats  = char_calloc2(nstat,5) ;
  *nets   = char_calloc2(nstat,3) ;
  *stlats = float_calloc(nstat) ;
  *stlons = float_calloc(nstat) ;
  for(jstat=0; jstat<nstat; jstat++)
    {
      fgets(line,64,fp);
      ns = sscanf(line, "%s%s%f%f", (*stats)[jstat], (*nets)[jstat],
		  &(*stlats)[jstat], &(*stlons)[jstat]);
      check_scan(4, ns, file, fp) ;
    }
  fclose(fp);
  return(nstat);
}

int 
r_scr_dat_fil_list(char* file, char*** stats, char*** nets, float** stlats, float** stlons)
{
  FILE* fp;  
  int   ns, jstat, nstat,i,k,flag;
  char  dum1[FSIZE],dum2[FSIZE],stat[FSIZE],net[FSIZE];
  float dum3,dum4,dum5,lat,lon;
  fp = openfile_rt(file,&k);
  *stats  = char_calloc2(k,5) ;
  *nets   = char_calloc2(k,3) ;
  *stlats = float_calloc(k) ;
  *stlons = float_calloc(k) ;
  nstat = 0 ;
  for(jstat=0; jstat<k; jstat++)
    {
      ns = fscanf(fp,"%s %s %s %s %f %f %f %f %f\n",
		  dum1,stat,net,dum2,&lat,&lon,&dum3,
		  &dum4,&dum5);
      check_scan(9, ns, file, fp) ;
      flag = 0;
      for(i=0;i<nstat;i++)
	if ((!strcmp(stat,(*stats)[i])) && (!strcmp(net,(*nets)[i])))
	  {
	    flag = 1;
	    break;
	  }
      if (flag)
	continue;
      strcpy((*stats)[nstat],stat) ;
      strcpy((*nets)[nstat],net)   ;
      (*stlats)[nstat] = lat ;
      (*stlons)[nstat] = lon ;
      nstat++;
    }
  fclose(fp);
  for(i=nstat;i<k;i++)
    {
      free((void*)(*stats)[i]);
      free((void*)(*nets)[i]);
    }
  return nstat ;
}

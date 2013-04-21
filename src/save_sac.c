/***************************************************************************
*
*	              W phase source inversion package 	            
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Université de Strasbourg / CNRS 
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
#include "rwsacs.h"


/* save displacement and header in sac file */
void 
save_sac(char *stnm, char *netwk, char *chan, float *lat, float *lon, sachdr *hdr,  double*   depval)
{
  int  i,nbc;
  char *sac_filename;
  
  sac_filename = char_alloc(FSIZE);

  strcpy(sac_filename, stnm);
  strcat(sac_filename, ".");
  strcat(sac_filename, netwk);
  strcat(sac_filename, ".");
  strcat(sac_filename, chan);
  strcat(sac_filename, ".SAC");

  nbc = strlen(stnm);
  strncpy(hdr->kstnm, stnm, nbc);
  for (i=nbc; i<8; i++)
    hdr->kstnm[i] = ' ';
  hdr->kstnm[8] = '\0';

  nbc = strlen(netwk);
  strncpy(hdr->knetwk,netwk,nbc);
  for (i=nbc; i<8; i++)
    hdr->knetwk[i] = ' ';
  hdr->knetwk[8] = '\0';  
  nbc = strlen(chan);
  strncpy(hdr->kcmpnm,chan, nbc);
  for (i=nbc; i<8; i++)
    hdr->kcmpnm[i] = ' ';
  hdr->kcmpnm[8] = '\0';
  
  hdr->stla = *lat;
  hdr->stlo = *lon;

  wsac(sac_filename, hdr, depval);
  
  free((void *)sac_filename);
}

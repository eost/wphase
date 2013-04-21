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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "proto_alloc.h"


/* External subroutine */
off_t lseek(int fildes, off_t offset, int whence);


/* Return the nb of samples of sacfiles at segm1/??/segm2 */
int 
get_length(char *segm1, char *segm2)
{
  char *sac_GF;
  int  npts;
  int  fd;

  sac_GF = char_alloc(FSIZE);
  strncpy(sac_GF, segm1, FSIZE);
  strcat(sac_GF, "/RR/");
  strcat(sac_GF, segm2);
  strcat(sac_GF, "Z.SAC");

  fd = open(sac_GF, O_RDONLY, "r");
  if(fd == -1)
	{
	 fflush(stdout);
         fprintf(stderr, "\n*** ERROR (get_length): sac_GF file %s not accessible\n", sac_GF);
         fprintf(stderr, "*** ... Exiting the program ... ***\n");
         fflush(stderr);
         exit(1);
	}

  lseek(fd, 316, SEEK_SET);
  read(fd, (void *)(&npts), sizeof(int));
  close(fd);
  free((void *)sac_GF);

  return(npts);
}

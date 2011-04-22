#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "../proto_alloc.h"


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
	  fprintf(stderr, "\n*** ERROR (get_length): sac_GF file %s not accesible\n", sac_GF);
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

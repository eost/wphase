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

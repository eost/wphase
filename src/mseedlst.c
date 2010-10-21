/************************************************************
*	List the sta net locid channel header fields in the 
*       first block of a miniseed record
*	Usage
*		mseedlst mseed_file
*       History
*             2010  Original Coding             Zacharie Duputel
*       License
*             Distributed under the same License as the libmseed source
*             and binaries, used here under permission by Zacharie Duputel,
*             Luis Rivera and Hiroo Kanamori
*
*************************************************************/

#include <stdio.h>

int get_channel_name(char *msfil,char *netwk,char *stanm,char *loc, char *cmpnm);

int 
main(int argc, char **argv)
{
  char netwk[9],stanm[9],locid[9],cmpnm[9];
  if(argc < 2) 
    {
      fprintf(stderr,"Usage: mseedlst mseed_file\n");
      return -1;
    }
  get_channel_name(argv[1],netwk,stanm,locid,cmpnm);
  printf("%9s %9s %9s %9s\n",stanm,netwk,locid,cmpnm);
  return 0;
}

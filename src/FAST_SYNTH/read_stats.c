#include <stdio.h>
#include <stdlib.h>
#include "../proto_alloc.h"
#include "../rwtextfiles.h"

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
  free((void*)line);
  return(nstat);
}

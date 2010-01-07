#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../proto_alloc.h"
#include "../rwtextfiles.h"
#include "../rwsacs.h"
#include "../read_i_files.h"


int   yyyymmdd2jjj(int year, int month, int day);

/* Read the cmt file */
void 
read_cmt(str_quake_params *eq, sachdr *hdr)
{
  get_cmtf(eq, 2) ;
  
  hdr->nzyear   = eq->ot_ye    ;
  hdr->nzhour = eq->ot_ho    ;
  hdr->nzmin  = eq->ot_mi    ;

/*   hdr->evla   = eq->pde_evla ; */
/*   hdr->evlo   = eq->pde_evlo ; */
/*   hdr->evdp   = eq->pde_evdp ; */

  hdr->evla   = eq->evla ;
  hdr->evlo   = eq->evlo ;
  hdr->evdp   = eq->evdp ;

  hdr->nzsec  = floor(eq->ot_se) ;
  hdr->nzmsec = (int) (1000.*(eq->ot_se - hdr->nzsec) + .5) ;

  hdr->nzjday = yyyymmdd2jjj(hdr->nzyear, eq->ot_mo, eq->ot_dm);
}

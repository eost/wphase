#include <stdio.h>
#include <string.h>

#include "proto_alloc.h"    /* proto_alloc.c    */
#include "rwtextfiles.h"

typedef struct
{
  int    n_repeat ;
  double dist_min, dist_max, an_sig, amw ;
  double *deltar, *ampr ;
  char   fcamp[FSIZE], fiamp[FSIZE], foamp[FSIZE] ;
  char   fdeln[FSIZE], foamp[FSIZE], frsta[FSIZE] ;
  char   fomwa[FSIZE], fsaclst[FSIZE], *saclst    ;
} structopt ;


void
get_param(structopt *opt)
{
  int  i,nc,tmp;
  FILE *ifile ;
  strcpy(opt->fcamp,"c_amp_azm_reg")     ;
  strcpy(opt->fiamp,"i_amp_azm_reg")     ;
  strcpy(opt->foamp,"o_amp_azm_reg")     ;
  strcpy(opt->fdeln,"f_delta_norm")      ;
  strcpy(opt->foamp,"o_amp_azm_reg_2")   ;
  strcpy(opt->frsta,"rejected_stations") ;
  strcpy(opt->fomwa,"o_mw-amp")          ;
  
  ifile = openfile_rt(opt->fcamp,&nc)  ;
  tmp   = fscanf(ifile,"%ld %ld %d %ld",opt->dist_min,\
		 opt->dist_max,opt->n_repeat,opt->an_sig) ;
  check_scan(4, tmp, opt->fcamp, ifile);
  fclose(ifile);

  ifile = openfile_rt(opt->fdeln,&nc)  ;
  opt->deltar = double_alloc(nc) ;
  opt->ampr   = double_alloc(nc) ;
  for(i=0; i<nc; i++)
    {
      tmp = fscanf(ifile,"%ld %ld", \
		   &(opt->deltar[i]),&(opt->ampr(i]));
      check_scan(2, tmp, opt->fdeln, ifile) ;
    }
  fclose(ifile) ;

  ifile = openfile_rt(opt->fiamp,&nc)     ;
  tmp   = fscanf(ifile,"%s",opt->fsaclst) ;
  check_scan(1, tmp, opt->fiamp, ifile)   ;
  tmp   = fscanf(ifile,"%ld",opt->amw)    ;
  check_scan(1, tmp, opt->fiamp, ifile)   ;
  fclose(ifile) ;

  ifile = openfile_rt(opt->fsaclst,&nc)     ;
  for(i=0; i<nc; i++)
    {
      *** CONTINUE HERE ***
    }
}

int 
main(int argc, char *argv[])
{
  structopt opt;
  
  /* Set input params */

  get_param(&opt);
  
}

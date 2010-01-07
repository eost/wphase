/***************************************************************************
 * Mon premier devoir utilisant libmseed
 *   ** lire un fichier seed et le convertir en sac
 ***************************************************************************/
//#include <stdlib.h>
//#include <time.h>

#include <stdio.h>
#include <libmseed.h>
#include <string.h>

#include "proto_alloc.h"
#include "rwsacs.h"



void fill_hdr_char(char *o_c, char* i_c,int N);

int check_hdr(sachdr *hdr, MSRecord *msr);

void fill_sac(sachdr *hdr, double *o_x, MSRecord *msr, int first, int npts);



int
main (int argc, char **argv)
{
  MSRecord *msr = 0;
  

  int   dataflag   = 1 ;
  int   totalrecs  = 0 ;
  int   totalsamps = 0 ; 
  int   retcode        ;
  
  short int verbose   =  0 ;
  short int basicsum  =  0 ;
  int   reclen        = -1 ;
  char *i_file     =  0 ;
  char *o_file     =  0 ;
  
  double *o_x;
  sachdr o_hd;

  /* Input parameters */
  i_file = argv[1];  
  o_file = argv[2];  
  
  
  /* Allocate memory  */
  o_x = double_alloc((int)__LEN_SIG__);
  hdr_alloc(&o_hd);
    
  /* Loop over the input file */
  while ( (retcode = ms_readmsr (&msr, i_file, reclen, NULL, NULL, 1,
				 dataflag, verbose)) == MS_NOERROR )
    {
      fill_sac(&o_hd, o_x, msr, 0, msr->numsamples);
      totalrecs++;
      totalsamps += msr->samplecnt;
    }
  
  if ( retcode != MS_ENDOFFILE )
    ms_log (2, "Cannot read %s: %s\n", i_file, ms_errorstr(retcode));
  
  /* Make sure everything is cleaned up */
  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);
  
  if ( basicsum )
    ms_log (1, "Records: %d, Samples: %d\n", totalrecs, totalsamps);
  whdrsac(o_file, &o_hd);
  whdrsac(o_file, &o_hd);
  wdatsac(o_file, &o_hd, o_x);
  
  free((void*)o_x) ;
    
  return 0;
}  /* End of main() */



void 
fill_hdr_char(char *o_c, char* i_c,int N)
{
  int  i,n;
  n = strlen(i_c);
  strncpy(o_c,i_c,n);
  for(i=n; i<(N-1); i++)
    o_c[i] = ' ';
  o_c[N] = '\0' ;
}


int 
check_hdr(sachdr *hdr, MSRecord *msr)
{
  int N;
  BTime    bt_beg;
  hptime_t hpt_beg;
  

  /* Check time-related fields */
  if (hdr->delta != (float)(1./(msr->samprate)))
    return 1;
  
  hpt_beg = msr->starttime -  MS_EPOCH2HPTIME((hdr->npts)/msr->samprate);
  ms_hptime2btime (hpt_beg,&bt_beg);

  if (bt_beg.year  != (uint16_t) hdr->nzyear )
    return 1;
  if (bt_beg.day   != (uint16_t) hdr->nzjday )
    return 1;
  if (bt_beg.hour  !=  (uint8_t) hdr->nzhour )
    return 1;
  if (bt_beg.min   !=  (uint8_t) hdr->nzmin  )
    return 1;
  if (bt_beg.sec   !=  (uint8_t) hdr->nzsec  )
    return 1;
  if ((int) ((double)bt_beg.fract/10.) != (hdr->nzmsec) )
    return 1;

  /* Check stream identifiers  */
  N = strlen(msr->station);
  if (strncmp(hdr->kstnm,msr->station,N) != 0)
    return 1;
  N = strlen(msr->network);
  if (strncmp(hdr->knetwk,msr->network,N) != 0)
    return 1;
  N = strlen(msr->location);
  if (strncmp(hdr->khole,msr->location,N) != 0)
    return 1;
  N = strlen(msr->channel);
  if (strncmp(hdr->kcmpnm,msr->channel,N) != 0)
    return 1;

  return 0;
}


void 
fill_sac(sachdr *hdr, double *o_x, MSRecord *msr, int first, int npts)
{
  int  i ;
  BTime t_beg;

  ms_hptime2btime(msr->starttime,&t_beg); 

  if (hdr->npts == -12345)
    {
      hdr->delta  = (float)(1./(msr->samprate));
      hdr->b      = (float)0.;
      hdr->e      = (float)((double)(msr->numsamples-1)/msr->samprate);
      hdr->nzyear = (int) t_beg.year ;
      hdr->nzjday = (int) t_beg.day  ;
      hdr->nzhour = (int) t_beg.hour ;
      hdr->nzmin  = (int) t_beg.min  ;
      hdr->nzsec  = (int) t_beg.sec  ;
      hdr->nzmsec = (int) ((double)t_beg.fract/10.)  ;
      hdr->npts   = (int) msr->numsamples ;
      hdr->iztype = 9 ;
      fill_hdr_char(hdr->kstnm,msr->station,9) ;
      fill_hdr_char(hdr->knetwk,msr->network,9);
      fill_hdr_char(hdr->khole,msr->location,9);
      fill_hdr_char(hdr->kcmpnm,msr->channel,9);
      for(i=first; i<msr->numsamples; i++)
	o_x[i] = (double)((int*)msr->datasamples)[i] ;
    }
  else if (first != 0)
    {
      fprintf(stderr,"ERROR...first must be 0");
      exit(1);
    }
  else
    {
      if (check_hdr(hdr, msr))
	{
	  fprintf(stderr,"ERROR...");
	  exit(1);
	}
      for(i=0; i<msr->numsamples; i++)
	o_x[i+hdr->npts] = (float)((int*)msr->datasamples)[i] ;
      hdr->npts += msr->numsamples;
      hdr->e     = (float)((double)(hdr->npts-1)/msr->samprate);
    }
}

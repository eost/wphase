#include <stdio.h>
#include <libmseed.h>
#include <string.h>

#include <time.h>
#include <locale.h>

#include "proto_alloc.h"
#include "rwsacs.h"
#include "read_i_files.h"
#include "rwtextfiles.h"

int read_ms_file(str_quake_params *eq, int year, int yday, hptime_t t0_in, \
		 hptime_t t1_in, sachdr *o_hd, double *o_x);
void get_msname(sachdr *hdr, str_quake_params *eq, int year, int day, char *msfil);
void fill_hdr_char(char *o_c, char* i_c,int N);
int  check_hdr(sachdr *hdr, MSRecord *msr);
int  fill_sac(sachdr *hdr, double *o_x, MSRecord *msr);

int
rmseed(str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1, \
       double *o_x, sachdr *o_hd)
{
  int      flag   ;
  hptime_t t0_in, t1_in ;
  

  /* Initialize variables */
  t0_in = MS_EPOCH2HPTIME(i_t0);
  t1_in = MS_EPOCH2HPTIME(i_t1);
  o_hd->npts = -12345 ;

  /* Loop over the input file */
  flag = read_ms_file(eq, tm0->tm_year+1900, tm0->tm_yday+1, t0_in, t1_in, o_hd, o_x) ;

  /* printf("%s : Records: %d, Samples: %d\n", msfil, totalrecs, totalsamps); */
  /* printf("RMSEED FLAG==%d\n",flag); */
  return flag;
} 


int
read_ms_file(str_quake_params *eq, int year, int yday, hptime_t t0_in, \
	     hptime_t t1_in, sachdr *o_hd, double *o_x)
{
  int  flag = 1       ;
  int  retcode        ;
  int  totalrecs  = 0 ;
  int  totalsamps = 0 ;
  char msfil[FSIZE]   ;
  hptime_t t0 ;
  MSRecord *msr = 0  ;
  FILE     *fid ;
  
  get_msname(o_hd, eq, year, yday, msfil);
  if ((fid = fopen (msfil,"r"))==NULL)
    {
      fprintf(stderr,"Warning: (rmseed) Cannot read %s (file rejected)\n",msfil);
      flag = 1 ;
    }
  else
    {
      fclose(fid);
      while ( (retcode = ms_readmsr (&msr, msfil, -1, NULL, NULL, 1, 1, 0)) == MS_NOERROR )
	{
	  t0 = msr->starttime+(double)msr->numsamples*(double)HPTMODULUS/msr->samprate;
	  if (t0_in <= t0 && t1_in >= msr->starttime)
	    {
	      if (fill_sac(o_hd, o_x, msr))
		{
		  fprintf(stderr," Warning: (rmseed) file %s rejected\n",msfil);
		  break ;
		}
	    }
	  else if ( t1_in < msr->starttime )
	    {
	      if (o_hd->npts ==  -12345) /* Data missing */
		flag = 1;
	      else
		{
		  retcode = MS_ENDOFFILE ;
		  flag = 0 ;
		}
	      break    ;
	    }
	  totalrecs++ ;
	  totalsamps += msr->samplecnt ;
	}
      if ( retcode != MS_ENDOFFILE)
	{
	  fprintf(stderr,"Warning: (rmseed) Cannot read %s (file rejected)\n",msfil);	  
	  flag = 1 ;
	}
    }
  /* Make sure everything is cleaned up */
  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);

  if (retcode == MS_ENDOFFILE  && flag)
    {
      fprintf(stderr,"Warning: (rmseed) file %s is not complete, reading next file...\n",msfil);           
      if (yday == 365)
	flag = read_ms_file(eq, year+1, 0, t0_in, t1_in, o_hd, o_x);
      else
	flag = read_ms_file(eq, year, yday+1, t0_in, t1_in, o_hd, o_x);
      if (flag == 0)
	{
	  fprintf(stderr,"             ... No problem encountered while reading next file\n");           
	}
    }

  return flag ;
}
void 
get_msname(sachdr *hdr, str_quake_params *eq, int year, int day, char *msfil)
{
  int i;
  char knetwk[9], kstnm[9], kcmpnm[9], khole[9];
  
  i = nbchar(hdr->knetwk) ;
  strncpy(knetwk,hdr->knetwk,i);
  knetwk[i] = '\0';
  i = nbchar(hdr->kstnm) ;
  strncpy(kstnm,hdr->kstnm,i);
  kstnm[i] = '\0';
  i = nbchar(hdr->kcmpnm) ;
  strncpy(kcmpnm,hdr->kcmpnm,i);
  kcmpnm[i] = '\0';
  i = nbchar(hdr->khole) ;
  strncpy(khole,hdr->khole,i);
  khole[i] = '\0';

  if (!strcmp(hdr->khole,"--"))
    sprintf(msfil,"%s%s/%s/%s.%s.%s.%s.%d.%d",eq->seed,knetwk,kstnm,
	    knetwk,kstnm,"",kcmpnm,year,day);
  else
    sprintf(msfil,"%s%s/%s/%s.%s.%s.%s.%d.%d",eq->seed,knetwk,kstnm,
	    knetwk,kstnm,khole,kcmpnm,year,day);
}

void 
fill_hdr_char(char *o_c, char* i_c,int N)
{
  int  i,n;
  n = strlen(i_c);
  strncpy(o_c,i_c,n);
  for(i=n; i<(N-1); i++)
    o_c[i] = ' ';
  o_c[N-1] = '\0' ;
}


void
error_mseed(char *field)
{
  fprintf(stderr,"Warning: (rmseed) discontinuous data (%s does not match)",field);
}

int 
check_hdr(sachdr *hdr, MSRecord *msr)
{
  int N;
  BTime    bt_beg;
  hptime_t hpt_beg;
  

  /* Check time-related fields */
  if ((int)(hdr->delta * 1000.+0.5) != (int)((float)(1./(msr->samprate)) * 1000.+0.5))
    {
      error_mseed("sample period");
      return 1;
    }
  
  
  hpt_beg = msr->starttime -  MS_EPOCH2HPTIME((hdr->npts)/msr->samprate);
  ms_hptime2btime (hpt_beg,&bt_beg);

  if (bt_beg.year  != (uint16_t) hdr->nzyr )
    {
      error_mseed("year");
      return 1;
    }
  if (bt_beg.day   != (uint16_t) hdr->nzjday )
    {
      error_mseed("day");
      return 1;
    }
  if (bt_beg.hour  !=  (uint8_t) hdr->nzhour )
    {
      error_mseed("hour");
      return 1;
    }
  if (bt_beg.min   !=  (uint8_t) hdr->nzmin  )
    {
      error_mseed("min.");
      /*       printf("min\n"); */
      /*       printf("bt_beg = %d\n",bt_beg.min); */
      /*       printf("hdr    = %d\n",hdr->nzmin); */
      /*       ms_hptime2btime (msr->starttime,&bt_beg); */
      /*       printf("npts = %d\n",hdr->npts); */
      /*       printf("bt_ini = %d %d %d:%d:%d.%d\n",bt_beg.year,bt_beg.day,bt_beg.hour,bt_beg.min,bt_beg.sec,bt_beg.fract/10.); */
      /*       printf("hdr    = %d %d %d:%d:%d.%d\n",hdr->nzyr,hdr->nzjday,hdr->nzhour,hdr->nzmin,hdr->nzsec,hdr->nzmsec); */
      return 1;
    }
  if (bt_beg.sec   !=  (uint8_t) hdr->nzsec  )
    {
      error_mseed("sec");
      return 1;
    }
  if ((int) ((double)bt_beg.fract/10.) != (hdr->nzmsec) )
    {
      error_mseed("msec");
      return 1;
    }

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


int 
fill_sac(sachdr *hdr, double *o_x, MSRecord *msr)
{
  int  i ;
  BTime t_beg;


  if (hdr->npts == -12345)
    {
      ms_hptime2btime(msr->starttime,&t_beg); 
      hdr->delta  = (float)(1./(msr->samprate));
      hdr->b      = (float)0.;
      hdr->e      = (float)((double)(msr->numsamples-1)/msr->samprate);
      hdr->nzyr   = (int) t_beg.year ;
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
      for(i=0; i<msr->numsamples; i++)
	o_x[i] = (double)((int*)msr->datasamples)[i] ;
    }
  else
    {
      if (check_hdr(hdr, msr))
	return 1;
      for(i=0; i<msr->numsamples; i++)
	o_x[i+hdr->npts] = (float)((int*)msr->datasamples)[i] ;
      hdr->npts += msr->numsamples;
      hdr->e     = (float)((double)(hdr->npts-1)/msr->samprate);
    }
  return 0;
}

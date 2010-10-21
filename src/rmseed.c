/****************************************************************
*	W phase package - subroutines to read demultiplexed 
*                         miniseed records 
*       History
*             2010  Original Coding          Zacharie Duputel
*       License
*             Distributed under the same License as the libmseed 
*             source and binaries, used here under permission by 
*             Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#include <stdio.h>
#include <libmseed.h>
#include <string.h>

#include <time.h>
#include <locale.h>

#include "proto_alloc.h"
#include "rwsacs.h"
#include "read_i_files.h"
#include "rwtextfiles.h"

int get_channel_name(char *msfil,char *netwk,char *stanm,char *loc, char *cmpnm);
int check_hdr(sachdr *hdr,MSRecord *msr);
int fill_sac(sachdr *hdr,double *o_x,MSRecord *msr);
int read_mseed_file(char *msfil,str_quake_params *eq,int year,int yday, 
		     hptime_t t0_in,hptime_t t1_in,sachdr *o_hdr,double *o_x);
int rmseed_file(char *mseedfil,str_quake_params *eq,struct tm *tm0,int i_t0, 
		int i_t1,double *o_x,sachdr *o_hdr);
int read_mseed(str_quake_params *eq,int year,int yday,hptime_t t0_in,
	       hptime_t t1_in,sachdr *o_hd,double *o_x);
int rmseed(str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1,
	   double *o_x, sachdr *o_hdr);


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
  if ((int)(hdr->delta*10000.+0.5)!=(int)((float)(1./(msr->samprate))*10000.+0.5))
    {
      error_mseed("sample period");
      return 1;
    }
  hpt_beg = msr->starttime -  MS_EPOCH2HPTIME((hdr->npts)/msr->samprate);
  ms_hptime2btime (hpt_beg,&bt_beg);
  if (bt_beg.year != (uint16_t) hdr->nzyear )
    {
      error_mseed("year");
      return 1;
    }
  if (bt_beg.day  != (uint16_t) hdr->nzjday )
    {
      error_mseed("day");
      return 1;
    }
  if (bt_beg.hour !=  (uint8_t) hdr->nzhour )
    {
      error_mseed("hour");
      return 1;
    }
  if (bt_beg.min  !=  (uint8_t) hdr->nzmin  )
    {
      error_mseed("min.");
      return 1;
    }
  if (bt_beg.sec  !=  (uint8_t) hdr->nzsec  )
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


/* Fill hdr and o_x from msr */
int 
fill_sac(sachdr *hdr, double *o_x, MSRecord *msr)
{
  int  i ;
  BTime t_beg;


  if (hdr->npts == -12345) /* First assignation */
    {
      ms_hptime2btime(msr->starttime,&t_beg); 
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



/* Miniseed filename */
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
    sprintf(msfil,"%s%s/%s/%s.%s.%s.%s.%04d.%03d",eq->seed,knetwk,kstnm,
	    knetwk,kstnm,"",kcmpnm,year,day);
  else
    sprintf(msfil,"%s%s/%s/%s.%s.%s.%s.%04d.%03d",eq->seed,knetwk,kstnm,
	    knetwk,kstnm,khole,kcmpnm,year,day);
}

/*******************************************/
/* Extract channel name in the first block */
int 
get_channel_name(char *msfil,char *netwk,char *stanm,char *loc, char *cmpnm)
{
  int tmp,retcode,flag=0;
  MSRecord *msr = 0 ;
  FILE *fid;
  if((fid=fopen(msfil,"r"))==NULL)
    {
      fprintf(stderr,"Warning: (rmseed) Cannot read %s (file rejected)\n",msfil);
      flag = 1;
    }
  else
    {
      fclose(fid);
      if((retcode=ms_readmsr (&msr, msfil, -1, NULL, NULL, 1, 0, 0))!=MS_NOERROR)
	{ 
	  fprintf(stderr,"Warning: (rmseed) Cannot read %s (file rejected)\n",msfil);
	  flag = 1;
	}
    }
  if(!flag)
    {
      strcpy(netwk, msr->network) ;
      strcpy(stanm, msr->station) ;
      tmp = nbchar(msr->location) ;
      if ( tmp == 0 )
	strcpy(loc, "--") ;
      else
	strcpy(loc,msr->location) ;
      strcpy(cmpnm, msr->channel) ;
    }
  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);
  return flag;
}

/***********************************************************/
/* Main subroutine                                         */
/* miniseed filename is specified in the first argument    */

/* Read mseed file, fill o_hdr struct and o_x to have data 
 *  between t0_in and t1_in                                */
int
read_mseed_file(char *msfil, str_quake_params *eq, int year, int yday, hptime_t t0_in, \
		hptime_t t1_in, sachdr *o_hdr, double *o_x)
{
  int  flag = 1       ;
  int  retcode        ;
  int  totalrecs  = 0 ;
  int  totalsamps = 0 ;
  hptime_t t0 ;
  MSRecord *msr = 0  ;
  FILE     *fid ;
  
  if ((fid = fopen (msfil,"r"))==NULL)
    {
      fprintf(stderr,"Warning: (rmseed) Cannot read %s (file rejected)\n",msfil);
      flag = 1 ;
    }
  else
    {
      fclose(fid);
      /* loop over miniseed blocks */
      while ( (retcode = ms_readmsr (&msr, msfil, -1, NULL, NULL, 1, 1, 0)) == MS_NOERROR )
	{ 
	  t0 = msr->starttime+(double)msr->numsamples*(double)HPTMODULUS/msr->samprate;
	  if (t0_in <= t0 && t1_in >= msr->starttime) 
	    {
	      if (fill_sac(o_hdr, o_x, msr)) /* fill o_hdr and o_x */
		{
		  fprintf(stderr," Warning: (rmseed) file %s rejected\n",msfil);
		  break ;
		}
	    }
	  else if ( t1_in < msr->starttime )
	    {
	      if (o_hdr->npts ==  -12345) /* Data missing */
		flag = 1;
	      else                        /* Stop reading */
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
  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0); /* Make sure everything is cleaned up */
  return flag ;
}

/* Interface for trim_sac_files_qrt.c */
int
rmseed_file(char *mseedfil, str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1, \
	  double *o_x, sachdr *o_hdr)
{
  int      flag         ;
  hptime_t t0_in, t1_in ;
  /* Initialize variables */
  t0_in = MS_EPOCH2HPTIME(i_t0) ;
  t1_in = MS_EPOCH2HPTIME(i_t1) ;
  o_hdr->npts = -12345 ;
  /* Loop over the input file */
  flag = read_mseed_file(mseedfil,eq, tm0->tm_year+1900, tm0->tm_yday+1, t0_in, t1_in, o_hdr, o_x) ;
  return flag;
} 


/***********************************************************/
/* Main subroutine                                         */
/* miniseed filename is obtained from get_msname           */

/* Read mseed file, fill o_hdr struct and o_x to have data 
 *  between t0_in and t1_in                                */
int
read_mseed(str_quake_params *eq, int year, int yday, hptime_t t0_in,
	     hptime_t t1_in, sachdr *o_hdr, double *o_x)
{
  int  flag = 1       ;
  int  retcode        ;
  int  totalrecs  = 0 ;
  int  totalsamps = 0 ;
  char msfil[FSIZE]   ;
  hptime_t t0 ;
  MSRecord *msr = 0  ;
  FILE     *fid ;
  
  get_msname(o_hdr, eq, year, yday, msfil); /* Get miniseed filename */
  if ((fid = fopen (msfil,"r"))==NULL)
    {
      fprintf(stderr,"Warning: (rmseed) Cannot read %s (file rejected)\n",msfil);
      flag = 1 ;
    }
  else
    {
      fclose(fid);
      /* loop over miniseed blocks */
      while ( (retcode = ms_readmsr (&msr, msfil, -1, NULL, NULL, 1, 1, 0)) == MS_NOERROR )
	{ 
	  t0 = msr->starttime+(double)msr->numsamples*(double)HPTMODULUS/msr->samprate;
	  if (t0_in <= t0 && t1_in >= msr->starttime) 
	    {
	      if (fill_sac(o_hdr, o_x, msr)) /* fill o_hdr and o_x */
		{
		  fprintf(stderr," Warning: (rmseed) file %s rejected\n",msfil);
		  break ;
		}
	    }
	  else if ( t1_in < msr->starttime )
	    {
	      if (o_hdr->npts ==  -12345) /* Data missing */
		flag = 1;
	      else                        /* Stop reading */
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

  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0); /* Make sure everything is cleaned up */
  if (retcode == MS_ENDOFFILE  && flag)     /* Reading next file */
    {
      fprintf(stderr,"Warning: (rmseed) file %s is not complete, reading next file...\n",msfil);           
      if (yday == 365)
	flag = read_mseed(eq, year+1, 0, t0_in, t1_in, o_hdr, o_x);
      else
	flag = read_mseed(eq, year, yday+1, t0_in, t1_in, o_hdr, o_x);
      if (flag == 0)
	{
	  fprintf(stderr,"             ... No problem encountered while reading next file\n");           
	}
    }

  return flag ;
}

/* Interface for trim_sac_files_qrt.c */
int
rmseed(str_quake_params *eq, struct tm *tm0, int i_t0, int i_t1,
	    double *o_x, sachdr *o_hdr)
{
  int      flag         ;
  hptime_t t0_in, t1_in ;
  /* Initialize variables */
  t0_in = MS_EPOCH2HPTIME(i_t0) ;
  t1_in = MS_EPOCH2HPTIME(i_t1) ;
  o_hdr->npts = -12345 ;
  /* Loop over the input file */
  flag = read_mseed(eq, tm0->tm_year+1900, tm0->tm_yday+1, t0_in, t1_in, o_hdr, o_x) ;
  return flag;
} 


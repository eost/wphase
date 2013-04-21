/****************************************************************
*	W phase package - sac files manipulation subroutines
*                                           
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include "proto_alloc.h"
#include "rwsacs.h"

/***************************************************/
/*                 v=float_alloc(n)                */
/***************************************************/
/* Allocates memory for a sac header               */
/* Default header version number = NVHDR           */
/* Default type of file = ITIME (time series file) */

void 
hdr_alloc(sachdr *hdr)
{
  int i;
  hdr->delta = -12345. ;
  hdr->depmin=-12345. ;
  hdr->depmax=-12345. ;
  hdr->scale=-12345. ;
  hdr->odelta=-12345. ;
  hdr->b=-12345. ;
  hdr->e=-12345. ;
  hdr->o=-12345. ;
  hdr->a=-12345. ;
  hdr->internal1=-12345. ;
  for (i=0 ; i<10; i++)
    hdr->t[i]=-12345. ;
  hdr->f=-12345. ;
  for (i=0 ; i<10; i++)
    hdr->resp[i]=-12345. ;
  hdr->stla=-12345. ;
  hdr->stlo=-12345. ;
  hdr->stel=-12345. ;
  hdr->stdp=-12345. ;
  hdr->evla=-12345. ;
  hdr->evlo=-12345. ;
  hdr->evel=-12345. ;
  hdr->evdp=-12345. ;
  hdr->mag=-12345. ;
  for (i=0 ; i<10; i++)
    hdr->user[i]=-12345. ;
  hdr->dist=-12345. ;
  hdr->az=-12345. ;
  hdr->baz=-12345. ;
  hdr->gcarc=-12345. ;
  hdr->internal2=-12345. ;
  hdr->internal3=-12345. ;
  hdr->depmen=-12345. ;
  hdr->cmpaz=-12345. ;
  hdr->cmpinc=-12345. ;
  hdr->xminimum=-12345. ;
  hdr->xmaximum=-12345. ;
  hdr->yminimum=-12345. ;
  hdr->ymaximum=-12345. ;  

  hdr->nzyear=-12345;  
  hdr->nzjday=-12345;  
  hdr->nzhour=-12345;
  hdr->nzmin=-12345;  
  hdr->nzsec=-12345;  
  hdr->nzmsec=-12345;  
  hdr->nvhdr=NVHDR;  
  hdr->norid=-12345;  
  hdr->nevid=-12345;
  hdr->npts=-12345; 
  hdr->internal4=-12345;  
  hdr->nwfid=-12345;
  hdr->nxsize=-12345;
  hdr->nysize=-12345;
  hdr->iftype=ITIME;
  hdr->idep=-12345;
  hdr->iztype=-12345;
  hdr->iinst=-12345;
  hdr->istreg=-12345;
  hdr->ievreg=-12345;
  hdr->ievtyp=-12345;
  hdr->iqual=-12345;
  hdr->isynth=-12345;
  hdr->imagtyp=-12345;
  hdr->imagsrc=-12345;
  hdr->leven=-12345;
  hdr->lpspol=-12345;
  hdr->lovrok=-12345;
  hdr->lcalda=-12345;

  /* hdr->kt    = char_alloc2(10,9) ; */
  /* hdr->kuser = char_alloc2(3,9) ;  */

  strcpy(hdr->kstnm,     "-12345  ");
  strcpy(hdr->kevnm,     "-12345          ");
  strcpy(hdr->khole,     "-12345  ");
  strcpy(hdr->ko,        "-12345  ");
  strcpy(hdr->ka,        "-12345  ");
  for (i=0 ; i<10 ; i++)
    strcpy(hdr->kt[i],   "-12345  ");
  strcpy(hdr->kf,        "-12345  ");
  for (i=0 ; i<3 ; i++)
    strcpy(hdr->kuser[i],"-12345  ");
  strcpy(hdr->kcmpnm,    "-12345  ");
  strcpy(hdr->knetwk,    "-12345  ");
  strcpy(hdr->kdatrd,    "-12345  ");
  strcpy( hdr->kinst,    "-12345  ");
}

/*********************************************/
/* Allocates memory for a tab of sac headers */
void 
hdr_tab(sachdr **hdr, int n)
{
  int i;
  if (((*hdr) = (sachdr*) malloc (n * sizeof(sachdr))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory");
      exit(1);    
    }
  for(i=0; i<n; i++)
    hdr_alloc(&((*hdr)[i])) ;
}


/**************************************************/
/*          rhdrsac(filename, h, ierror)          */
/**************************************************/
/* Read sac file header                           */
/* Input  : filename                              */
/*          ierror                                */
/*	     >if *ierror=1 (input) and if the     */
/*            file does not exist then an         */
/*	      error occur                         */
/*	     >if *ierror=0 (input) and if the     */
/*	      file does not exist, then ierror    */
/*	      is returned as *ierror = 1 (output) */
/*                                                */
/* Output : h : a struct sachdr                   */
void 
rhdrsac(char *file, sachdr *hdr, int *ierror)
{
  const int sz_of_f=4, sz_of_i=4;
  int i;
  FILE *f;
  
  if ((f=fopen(file,"rb"))==NULL)
    {
      if (*ierror == 1)	{
		fprintf(stderr,"ERROR opening %s to read header\n",file);
		exit(1); }
      *ierror = 1 ;
      return;
      
    }
  fread(&hdr->delta,sz_of_f,1,f);
  fread(&hdr->depmin,sz_of_f,1,f);
  fread(&hdr->depmax,sz_of_f,1,f);
  fread(&hdr->scale,sz_of_f,1,f);
  fread(&hdr->odelta,sz_of_f,1,f);
  fread(&hdr->b,sz_of_f,1,f);
  fread(&hdr->e,sz_of_f,1,f);
  fread(&hdr->o,sz_of_f,1,f);
  fread(&hdr->a,sz_of_f,1,f);
  fread(&hdr->internal1,sz_of_f,1,f);
  for (i=0 ; i<10; i++)
    fread(&hdr->t[i],sz_of_f,1,f);
  fread(&hdr->f,sz_of_f,1,f);
  for (i=0 ; i<10; i++)
    fread(&hdr->resp[i],sz_of_f,1,f);
  fread(&hdr->stla,sz_of_f,1,f);
  fread(&hdr->stlo,sz_of_f,1,f);
  fread(&hdr->stel,sz_of_f,1,f);
  fread(&hdr->stdp,sz_of_f,1,f);
  fread(&hdr->evla,sz_of_f,1,f);
  fread(&hdr->evlo,sz_of_f,1,f);
  fread(&hdr->evel,sz_of_f,1,f);
  fread(&hdr->evdp,sz_of_f,1,f);
  fread(&hdr->mag,sz_of_f,1,f);
  for (i=0 ; i<10; i++)
    fread(&hdr->user[i],sz_of_f,1,f);
  fread(&hdr->dist,sz_of_f,1,f);
  fread(&hdr->az,sz_of_f,1,f);
  fread(&hdr->baz,sz_of_f,1,f);
  fread(&hdr->gcarc,sz_of_f,1,f);
  fread(&hdr->internal2,sz_of_f,1,f);
  fread(&hdr->internal3,sz_of_f,1,f);
  fread(&hdr->depmen,sz_of_f,1,f);
  fread(&hdr->cmpaz,sz_of_f,1,f);
  fread(&hdr->cmpinc,sz_of_f,1,f);
  fread(&hdr->xminimum,sz_of_f,1,f);
  fread(&hdr->xmaximum,sz_of_f,1,f);
  fread(&hdr->yminimum,sz_of_f,1,f);
  fread(&hdr->ymaximum,sz_of_f,1,f);  
  fseek(f,7*sz_of_f,SEEK_CUR);
  fread(&hdr->nzyear,sz_of_i,1,f);  
  fread(&hdr->nzjday,sz_of_i,1,f);  
  fread(&hdr->nzhour,sz_of_i,1,f);
  fread(&hdr->nzmin,sz_of_i,1,f);  
  fread(&hdr->nzsec,sz_of_i,1,f);  
  fread(&hdr->nzmsec,sz_of_i,1,f);  
  fread(&hdr->nvhdr,sz_of_i,1,f);  
  fread(&hdr->norid,sz_of_i,1,f);  
  fread(&hdr->nevid,sz_of_i,1,f);
  fread(&hdr->npts,sz_of_i,1,f); 
  fread(&hdr->internal4,sz_of_i,1,f);  
  fread(&hdr->nwfid,sz_of_i,1,f);
  fread(&hdr->nxsize,sz_of_i,1,f);
  fread(&hdr->nysize,sz_of_i,1,f);
  fseek(f,1*sz_of_i,SEEK_CUR);
  fread(&hdr->iftype,sz_of_i,1,f);
  fread(&hdr->idep,sz_of_i,1,f);
  fread(&hdr->iztype,sz_of_i,1,f);
  fseek(f,1*sz_of_i,SEEK_CUR);
  fread(&hdr->iinst,sz_of_i,1,f);
  fread(&hdr->istreg,sz_of_i,1,f);
  fread(&hdr->ievreg,sz_of_i,1,f);
  fread(&hdr->ievtyp,sz_of_i,1,f);
  fread(&hdr->iqual,sz_of_i,1,f);
  fread(&hdr->isynth,sz_of_i,1,f);
  fread(&hdr->imagtyp,sz_of_i,1,f);
  fread(&hdr->imagsrc,sz_of_i,1,f);
  fseek(f,8*sz_of_i,SEEK_CUR);
  fread(&hdr->leven,sz_of_i,1,f);
  fread(&hdr->lpspol,sz_of_i,1,f);
  fread(&hdr->lovrok,sz_of_i,1,f);
  fread(&hdr->lcalda,sz_of_i,1,f);
  fseek(f,1*sz_of_i,SEEK_CUR);
  fgets(hdr->kstnm,9,f);
  fgets(hdr->kevnm,17,f);
  fgets(hdr->khole,9,f);
  fgets(hdr->ko,9,f);
  fgets(hdr->ka,9,f);
  for (i=0 ; i<10 ; i++)
    fgets(hdr->kt[i],9,f);
  fgets(hdr->kf,9,f);
  for (i=0 ; i<3 ; i++)
    fgets(hdr->kuser[i],9,f);
  fgets(hdr->kcmpnm,9,f);
  fgets(hdr->knetwk,9,f);
  fgets(hdr->kdatrd,9, f);
  fgets(hdr->kinst,9,f);
  /* e */
  hdr->e = hdr->b + ((float)(hdr->npts - 1)) * hdr->delta ; 
  fclose(f);
}



/**************************************************/
/*          rhdrsac(filename, h, data)            */
/**************************************************/
/* Read sac file data                             */ 
/*                                                */
/* Input  : filename                              */
/*          h : the struct sachdr of this         */
/*              file                              */
/*          ierror                                */
/*	     >if *ierror=1 (input) and if the     */
/*            file does not exist then an         */
/*	      error occur                         */
/*	     >if *ierror=0 (input) and if the     */
/*	      file does not exist, then ierror    */
/*	      is returned as *ierror = 1 (output) */
/*                                                */
/* Output : data : pointer to the proto_allocd       */
/*                 data array                     */
void 
rdatsac(char *file, sachdr *hdr, double *data, int *ierror)
{
  const int sz_of_f = 4;
  int i;
  float *d;
  FILE  *f;
  
  if ((f=fopen(file,"rb"))==NULL)
    {
      if (*ierror == 1)	
	{
	  fprintf(stderr,"ERROR opening %s to read data\n",file);
	  exit(1);	
	}
      *ierror = 1 ;
      return ;
    }

  fseek(f,632,SEEK_SET);
  d = float_alloc(hdr->npts) ;

  fread(d,sz_of_f,hdr->npts,f);
  for (i=0 ; i<hdr->npts ; i++)
    data[i] = (double)d[i]; /* conversion from float to double */
  /* e */
  hdr->e = hdr->b + ((float)(hdr->npts - 1)) * hdr->delta ; 
  free((void*)d);
  fclose(f);
}


/************************************************/
/*            whdr(f, h)                        */
/************************************************/
/* Write sac file header                        */
/* Input  : f : File descriptor                 */
/*          h : the struct sachdr               */
int
whdr(FILE *f, sachdr *hdr)
{
  const int sz_of_f = 4, sz_of_i = 4;
  int i,dumi;
  
  float dumf;
  char  dumc[9];
  
  /* Undefined values */
  dumi = -12345;
  dumf = -12345.0;
  strcpy(dumc,"-12345");
  
  if (sizeof(int) != 4) 
    {
      fprintf(stderr,"ERROR : int not in 32bits\n");
      return 1;
    }
  
  i = fseek(f,0,SEEK_SET);
  if (i!=0) 
    return 1;
  fwrite(&hdr->delta,sz_of_f,1,f);
  fwrite(&hdr->depmin,sz_of_f,1,f);
  fwrite(&hdr->depmax,sz_of_f,1,f);
  fwrite(&hdr->scale,sz_of_f,1,f);
  fwrite(&hdr->odelta,sz_of_f,1,f);
  fwrite(&hdr->b,sz_of_f,1,f);
  /* e */
  hdr->e = hdr->b + ((float)(hdr->npts - 1)) * hdr->delta ; 
  fwrite(&hdr->e,sz_of_f,1,f);
  fwrite(&hdr->o,sz_of_f,1,f);
  fwrite(&hdr->a,sz_of_f,1,f);
  fwrite(&hdr->internal1,sz_of_f,1,f);
  for (i=0 ; i<10; i++)
    fwrite(&hdr->t[i],sz_of_f,1,f);
  fwrite(&hdr->f,sz_of_f,1,f);
  for (i=0 ; i<10; i++)
    fwrite(&hdr->resp[i],sz_of_f,1,f);
  fwrite(&hdr->stla,sz_of_f,1,f);
  fwrite(&hdr->stlo,sz_of_f,1,f);
  fwrite(&hdr->stel,sz_of_f,1,f);
  fwrite(&hdr->stdp,sz_of_f,1,f);
  fwrite(&hdr->evla,sz_of_f,1,f);
  fwrite(&hdr->evlo,sz_of_f,1,f);
  fwrite(&hdr->evel,sz_of_f,1,f);
  fwrite(&hdr->evdp,sz_of_f,1,f);
  fwrite(&hdr->mag,sz_of_f,1,f);
  for (i=0 ; i<10; i++)
    fwrite(&hdr->user[i],sz_of_f,1,f);
  fwrite(&hdr->dist,sz_of_f,1,f);
  fwrite(&hdr->az,sz_of_f,1,f);
  fwrite(&hdr->baz,sz_of_f,1,f);
  fwrite(&hdr->gcarc,sz_of_f,1,f);
  fwrite(&hdr->internal2,sz_of_f,1,f);
  fwrite(&hdr->internal3,sz_of_f,1,f);
  fwrite(&hdr->depmen,sz_of_f,1,f);
  fwrite(&hdr->cmpaz,sz_of_f,1,f);
  fwrite(&hdr->cmpinc,sz_of_f,1,f);
  fwrite(&hdr->xminimum,sz_of_f,1,f);
  fwrite(&hdr->xmaximum,sz_of_f,1,f);
  fwrite(&hdr->yminimum,sz_of_f,1,f);
  fwrite(&hdr->ymaximum,sz_of_f,1,f);  
  for (i=0; i<7; i++)
    fwrite(&dumf,sz_of_f,1,f);
  /*fseek(f,7*sz_of_f,SEEK_CUR);*/
  fwrite(&hdr->nzyear,sz_of_i,1,f);  
  fwrite(&hdr->nzjday,sz_of_i,1,f);  
  fwrite(&hdr->nzhour,sz_of_i,1,f);
  fwrite(&hdr->nzmin,sz_of_i,1,f);  
  fwrite(&hdr->nzsec,sz_of_i,1,f);  
  fwrite(&hdr->nzmsec,sz_of_i,1,f);  
  fwrite(&hdr->nvhdr,sz_of_i,1,f);  
  fwrite(&hdr->norid,sz_of_i,1,f);  
  fwrite(&hdr->nevid,sz_of_i,1,f);
  fwrite(&hdr->npts,sz_of_i,1,f); 
  fwrite(&hdr->internal4,sz_of_i,1,f);  
  fwrite(&hdr->nwfid,sz_of_i,1,f);
  fwrite(&hdr->nxsize,sz_of_i,1,f);
  fwrite(&hdr->nysize,sz_of_i,1,f);
  fwrite(&dumi,sz_of_i,1,f);
  fwrite(&hdr->iftype,sz_of_i,1,f);
  fwrite(&hdr->idep,sz_of_i,1,f);
  fwrite(&hdr->iztype,sz_of_i,1,f);
  fwrite(&dumi,sz_of_i,1,f);  
  fwrite(&hdr->iinst,sz_of_i,1,f);
  fwrite(&hdr->istreg,sz_of_i,1,f);
  fwrite(&hdr->ievreg,sz_of_i,1,f);
  fwrite(&hdr->ievtyp,sz_of_i,1,f);
  fwrite(&hdr->iqual,sz_of_i,1,f);
  fwrite(&hdr->isynth,sz_of_i,1,f);
  fwrite(&hdr->imagtyp,sz_of_i,1,f);
  fwrite(&hdr->imagsrc,sz_of_i,1,f);
  for (i=0; i<8; i++)
    fwrite(&dumi,sz_of_i,1,f);
  fwrite(&hdr->leven,sz_of_i,1,f);
  fwrite(&hdr->lpspol,sz_of_i,1,f);
  fwrite(&hdr->lovrok,sz_of_i,1,f);
  fwrite(&hdr->lcalda,sz_of_i,1,f);
  fwrite(&dumi,sz_of_i,1,f);
  fputs(hdr->kstnm,f); 
  for (i=0 ; i<16; i++) 
    fprintf(f,"%c",hdr->kevnm[i]);
  fputs(hdr->khole,f); 
  fputs(hdr->ko,f); 
  fputs(hdr->ka,f); 
  for (i=0 ; i<10 ; i++) 
    fputs(hdr->kt[i],f); 
  fputs(hdr->kf,f); 
    for (i=0 ; i<3 ; i++)
    fputs(hdr->kuser[i],f);
  fputs(hdr->kcmpnm,f);
  fputs(hdr->knetwk,f);
  fputs(hdr->kdatrd, f);
  fputs(hdr->kinst,f);
  i = fseek(f,0,SEEK_SET);
  if (i!=0) 
    return 1;
  return 0;
}


void 
minimax(float *y, int np, float *ymin, float *ymax) 
{
  int i ;
  *ymin = y[0] ;
  *ymax = y[0] ;
  for (i=1 ; i<np ; i++)
    {
      if (y[i] > *ymax)
	*ymax = y[i] ;
      if (y[i] < *ymin)
	*ymin = y[i] ;
    }
}


/***************************************************/
/*           wdat(f, h, data)                      */
/***************************************************/
/* Write data in sac file                          */
/* Input  : f : file descriptor                    */
/*          h : the struct sachdr                  */
/*          data : array of data points            */
int
wdat(FILE *f, sachdr *hdr, double *data)
{
  const int sz_of_f = 4, sz_of_i = 4;
  int    i,npts;
  float *d;

  i = fseek(f,316,SEEK_SET);
  fread(&npts,sz_of_i,1,f);
  if (i!=0)
    return 1;
  i = fseek(f,632,SEEK_SET);
  if (i!=0) 
    return 1;
  d = float_alloc(hdr->npts) ;
  for (i=0 ; i<hdr->npts ; i++)
      d[i] = (float)data[i];
  fwrite(d,sz_of_f,hdr->npts,f);
  /* rewrite depmin,depmax */
  minimax(d,hdr->npts,&hdr->depmin,&hdr->depmax) ;
  i = fseek(f,4,SEEK_SET);
  if (i!=0) 
    return 1;
  fwrite(&hdr->depmin,sz_of_f,1,f);
  fwrite(&hdr->depmax,sz_of_f,1,f);    
  /* rewrite e    */
  hdr->e = hdr->b + ((float)(hdr->npts - 1)) * hdr->delta ; 
  i = fseek(f,24,SEEK_SET);
  if (i!=0) 
    return 1;
  fwrite(&hdr->e,sz_of_f,1,f);
  free((void*)d);
  i = fseek(f,0,SEEK_SET);
  if (i!=0) 
    return 1;
  return 0;
}



/************************************************/
/*            whdrsac(filename, h)              */
/************************************************/
/* Overwrite sac file header                    */
/*   >If "filename" already exists, the header  */
/*    header is replaced but the data points    */
/*    are not affected.                         */
/*   >If "filename" does not exists, a new file */
/*    is created                                */
/*                                              */
/* Input  : filename                            */
/*          h : the struct sachdr               */
/* Output : header in the sac file "filename    */
void 
whdrsac(char *file, sachdr *hdr)
{
  FILE *f;
  
  if ((f=fopen(file,"r+b"))==NULL) 
    if ((f=fopen(file,"wb"))==NULL) /* If file does not exist */
      {
	fprintf(stderr,"ERROR opening %s to write sac header\n",file);
	exit(1); 
      }
  if (whdr(f, hdr))
    {
      fprintf(stderr,"ERROR writing header in file %s\n",file);
      exit(1);
    }
  fclose(f);
}


/***************************************************/
/*           wsac(filename, h, data)               */
/***************************************************/
/* Write sac file                                  */
/* Input  : filename                               */
/*          h : the struct sachdr                  */
/*          data : array of data points            */
void 
wsac(char *file, sachdr *hdr, double *data)
{
  FILE  *f;
  
  if ((f=fopen(file,"wb"))==NULL) 
    {
      fprintf(stderr,"ERROR opening file %s to write data\n",file);
      exit(1); 
    }
  if (whdr(f, hdr))
    {
      fprintf(stderr,"ERROR writing header in file %s\n",file);
      exit(1);
    }
  if (wdat(f,hdr,data))
    {
      fprintf(stderr,"ERROR writing data in file %s\n",file);
      exit(1);
    }
  fclose(f);
}



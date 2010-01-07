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
  int i;
  FILE *f;
  
  	 
  if (sizeof(int) != 4)    {
    fprintf(stderr,"ERROR : int not in 32bits\n");
    exit(1); }
  
  if ((f=fopen(file,"rb"))==NULL)
    {
      if (*ierror == 1)	{
	fprintf(stderr,"ERROR opening %s to read header\n",file);
	exit(1); }
      *ierror = 1 ;
      return;
      
    }
  fread(&hdr->delta,sizeof(float),1,f);
  fread(&hdr->depmin,sizeof(float),1,f);
  fread(&hdr->depmax,sizeof(float),1,f);
  fread(&hdr->scale,sizeof(float),1,f);
  fread(&hdr->odelta,sizeof(float),1,f);
  fread(&hdr->b,sizeof(float),1,f);
  fread(&hdr->e,sizeof(float),1,f);
  fread(&hdr->o,sizeof(float),1,f);
  fread(&hdr->a,sizeof(float),1,f);
  fread(&hdr->internal1,sizeof(float),1,f);
  for (i=0 ; i<10; i++)
    fread(&hdr->t[i],sizeof(float),1,f);
  fread(&hdr->f,sizeof(float),1,f);
  for (i=0 ; i<10; i++)
    fread(&hdr->resp[i],sizeof(float),1,f);
  fread(&hdr->stla,sizeof(float),1,f);
  fread(&hdr->stlo,sizeof(float),1,f);
  fread(&hdr->stel,sizeof(float),1,f);
  fread(&hdr->stdp,sizeof(float),1,f);
  fread(&hdr->evla,sizeof(float),1,f);
  fread(&hdr->evlo,sizeof(float),1,f);
  fread(&hdr->evel,sizeof(float),1,f);
  fread(&hdr->evdp,sizeof(float),1,f);
  fread(&hdr->mag,sizeof(float),1,f);
  for (i=0 ; i<10; i++)
    fread(&hdr->user[i],sizeof(float),1,f);
  fread(&hdr->dist,sizeof(float),1,f);
  fread(&hdr->az,sizeof(float),1,f);
  fread(&hdr->baz,sizeof(float),1,f);
  fread(&hdr->gcarc,sizeof(float),1,f);
  fread(&hdr->internal2,sizeof(float),1,f);
  fread(&hdr->internal3,sizeof(float),1,f);
  fread(&hdr->depmen,sizeof(float),1,f);
  fread(&hdr->cmpaz,sizeof(float),1,f);
  fread(&hdr->cmpinc,sizeof(float),1,f);
  fread(&hdr->xminimum,sizeof(float),1,f);
  fread(&hdr->xmaximum,sizeof(float),1,f);
  fread(&hdr->yminimum,sizeof(float),1,f);
  fread(&hdr->ymaximum,sizeof(float),1,f);  
  fseek(f,7*sizeof(float),SEEK_CUR);
  fread(&hdr->nzyear,sizeof(int),1,f);  
  fread(&hdr->nzjday,sizeof(int),1,f);  
  fread(&hdr->nzhour,sizeof(int),1,f);
  fread(&hdr->nzmin,sizeof(int),1,f);  
  fread(&hdr->nzsec,sizeof(int),1,f);  
  fread(&hdr->nzmsec,sizeof(int),1,f);  
  fread(&hdr->nvhdr,sizeof(int),1,f);  
  fread(&hdr->norid,sizeof(int),1,f);  
  fread(&hdr->nevid,sizeof(int),1,f);
  fread(&hdr->npts,sizeof(int),1,f); 
  fread(&hdr->internal4,sizeof(int),1,f);  
  fread(&hdr->nwfid,sizeof(int),1,f);
  fread(&hdr->nxsize,sizeof(int),1,f);
  fread(&hdr->nysize,sizeof(int),1,f);
  fseek(f,1*sizeof(int),SEEK_CUR);
  fread(&hdr->iftype,sizeof(int),1,f);
  fread(&hdr->idep,sizeof(int),1,f);
  fread(&hdr->iztype,sizeof(int),1,f);
  fseek(f,1*sizeof(int),SEEK_CUR);
  fread(&hdr->iinst,sizeof(int),1,f);
  fread(&hdr->istreg,sizeof(int),1,f);
  fread(&hdr->ievreg,sizeof(int),1,f);
  fread(&hdr->ievtyp,sizeof(int),1,f);
  fread(&hdr->iqual,sizeof(int),1,f);
  fread(&hdr->isynth,sizeof(int),1,f);
  fread(&hdr->imagtyp,sizeof(int),1,f);
  fread(&hdr->imagsrc,sizeof(int),1,f);
  fseek(f,8*sizeof(int),SEEK_CUR);
  fread(&hdr->leven,sizeof(int),1,f);
  fread(&hdr->lpspol,sizeof(int),1,f);
  fread(&hdr->lovrok,sizeof(int),1,f);
  fread(&hdr->lcalda,sizeof(int),1,f);
  fseek(f,1*sizeof(int),SEEK_CUR);
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
  int i;
  float *d;
  FILE  *f;
  
  if ((f=fopen(file,"rb"))==NULL)
    {
      if (*ierror == 1)	{
	fprintf(stderr,"ERROR opening %s to read data\n",file);
	exit(1);	}
      *ierror = 1 ;
      return;
    }

  fseek(f,632,SEEK_SET);
  d = float_alloc(hdr->npts) ;

  fread(d,sizeof(float),hdr->npts,f);
  for (i=0 ; i<hdr->npts ; i++)
    data[i] = (double)d[i]; /* conversion from float to double */
  /* e */
  hdr->e = hdr->b + ((float)(hdr->npts - 1)) * hdr->delta ; 
  free((void*)d);
  fclose(f);
}




/************************************************/
/*            whdrsac(filename, h)              */
/************************************************/
/* Write sac file header                        */
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
  int i,dumi;
  float dumf;
  char  dumc[9];
  FILE *f;
  
  /* Undefined values */
  dumi = -12345;
  dumf = -12345.0;
  strcpy(dumc,"-12345");
  
  if (sizeof(int) != 4) {
    fprintf(stderr,"ERROR : int not in 32bits\n");
    exit(1); }
  
  if ((f=fopen(file,"r+b"))==NULL)
    { /* If file does not exist */
      if ((f=fopen(file,"wb"))==NULL) {
	fprintf(stderr,"ERROR opening %s to write sac header\n",file);
	exit(1); }
    }
  fseek(f,0,SEEK_SET);
  fwrite(&hdr->delta,sizeof(float),1,f);
  fwrite(&hdr->depmin,sizeof(float),1,f);
  fwrite(&hdr->depmax,sizeof(float),1,f);
  fwrite(&hdr->scale,sizeof(float),1,f);
  fwrite(&hdr->odelta,sizeof(float),1,f);
  fwrite(&hdr->b,sizeof(float),1,f);
  /* e */
  hdr->e = hdr->b + ((float)(hdr->npts - 1)) * hdr->delta ; 
  fwrite(&hdr->e,sizeof(float),1,f);
  fwrite(&hdr->o,sizeof(float),1,f);
  fwrite(&hdr->a,sizeof(float),1,f);
  fwrite(&hdr->internal1,sizeof(float),1,f);
  for (i=0 ; i<10; i++)
    fwrite(&hdr->t[i],sizeof(float),1,f);
  fwrite(&hdr->f,sizeof(float),1,f);
  for (i=0 ; i<10; i++)
    fwrite(&hdr->resp[i],sizeof(float),1,f);
  fwrite(&hdr->stla,sizeof(float),1,f);
  fwrite(&hdr->stlo,sizeof(float),1,f);
  fwrite(&hdr->stel,sizeof(float),1,f);
  fwrite(&hdr->stdp,sizeof(float),1,f);
  fwrite(&hdr->evla,sizeof(float),1,f);
  fwrite(&hdr->evlo,sizeof(float),1,f);
  fwrite(&hdr->evel,sizeof(float),1,f);
  fwrite(&hdr->evdp,sizeof(float),1,f);
  fwrite(&hdr->mag,sizeof(float),1,f);
  for (i=0 ; i<10; i++)
    fwrite(&hdr->user[i],sizeof(float),1,f);
  fwrite(&hdr->dist,sizeof(float),1,f);
  fwrite(&hdr->az,sizeof(float),1,f);
  fwrite(&hdr->baz,sizeof(float),1,f);
  fwrite(&hdr->gcarc,sizeof(float),1,f);
  fwrite(&hdr->internal2,sizeof(float),1,f);
  fwrite(&hdr->internal3,sizeof(float),1,f);
  fwrite(&hdr->depmen,sizeof(float),1,f);
  fwrite(&hdr->cmpaz,sizeof(float),1,f);
  fwrite(&hdr->cmpinc,sizeof(float),1,f);
  fwrite(&hdr->xminimum,sizeof(float),1,f);
  fwrite(&hdr->xmaximum,sizeof(float),1,f);
  fwrite(&hdr->yminimum,sizeof(float),1,f);
  fwrite(&hdr->ymaximum,sizeof(float),1,f);  
  for (i=0; i<7; i++)
    fwrite(&dumf,sizeof(float),1,f);
  /*fseek(f,7*sizeof(float),SEEK_CUR);*/
  fwrite(&hdr->nzyear,sizeof(int),1,f);  
  fwrite(&hdr->nzjday,sizeof(int),1,f);  
  fwrite(&hdr->nzhour,sizeof(int),1,f);
  fwrite(&hdr->nzmin,sizeof(int),1,f);  
  fwrite(&hdr->nzsec,sizeof(int),1,f);  
  fwrite(&hdr->nzmsec,sizeof(int),1,f);  
  fwrite(&hdr->nvhdr,sizeof(int),1,f);  
  fwrite(&hdr->norid,sizeof(int),1,f);  
  fwrite(&hdr->nevid,sizeof(int),1,f);
  fwrite(&hdr->npts,sizeof(int),1,f); 
  fwrite(&hdr->internal4,sizeof(int),1,f);  
  fwrite(&hdr->nwfid,sizeof(int),1,f);
  fwrite(&hdr->nxsize,sizeof(int),1,f);
  fwrite(&hdr->nysize,sizeof(int),1,f);
  fwrite(&dumi,sizeof(int),1,f);
  fwrite(&hdr->iftype,sizeof(int),1,f);
  fwrite(&hdr->idep,sizeof(int),1,f);
  fwrite(&hdr->iztype,sizeof(int),1,f);
  fwrite(&dumi,sizeof(int),1,f);  
  fwrite(&hdr->iinst,sizeof(int),1,f);
  fwrite(&hdr->istreg,sizeof(int),1,f);
  fwrite(&hdr->ievreg,sizeof(int),1,f);
  fwrite(&hdr->ievtyp,sizeof(int),1,f);
  fwrite(&hdr->iqual,sizeof(int),1,f);
  fwrite(&hdr->isynth,sizeof(int),1,f);
  fwrite(&hdr->imagtyp,sizeof(int),1,f);
  fwrite(&hdr->imagsrc,sizeof(int),1,f);
  for (i=0; i<8; i++)
    fwrite(&dumi,sizeof(int),1,f);
  fwrite(&hdr->leven,sizeof(int),1,f);
  fwrite(&hdr->lpspol,sizeof(int),1,f);
  fwrite(&hdr->lovrok,sizeof(int),1,f);
  fwrite(&hdr->lcalda,sizeof(int),1,f);
  fwrite(&dumi,sizeof(int),1,f);
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
  fclose(f);
}


/***************************************************/
/*           wdatsac(filename, h, data)            */
/***************************************************/
/* Write data in sac file                          */
/*   > The sac file must already exist otherwise   */
/*     an error signal is sended.                  */
/* Input  : filename                               */
/*          h : the struct sachdr                  */
/*          data : array of data points            */
/* Output : data points in the sac file "filename" */
void 
wdatsac(char *file, sachdr *hdr, double *data)
{
  int    i,npts;
  float *d;
  FILE  *f;
  
    
  if ((f=fopen(file,"r+b"))==NULL) {
    fprintf(stderr,"ERROR opening file %s to write data\n",file);
    exit(1); }

  i = fseek(f,316,SEEK_SET);
  fread(&npts,sizeof(int),1,f);
  if (i!=0 || npts != hdr->npts) {
    fprintf(stderr,"ERROR writing data in file %s\n",file);
    fprintf(stderr,"npts header = %d, npts memory = %d\n",npts,hdr->npts);
    exit(1); }

  i = fseek(f,632,SEEK_SET);
  if (i!=0) {
    fprintf(stderr,"ERROR writing data in file %s\n",file);
    exit(1); }
  d = float_alloc(hdr->npts) ;
  for (i=0 ; i<hdr->npts ; i++)
      d[i] = (float)data[i];

  fwrite(d,sizeof(float),hdr->npts,f);
  free((void*)d);
  fclose(f);
}



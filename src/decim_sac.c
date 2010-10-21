#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include "decimate.h"
#include "rwsacs.h"
#include "proto_alloc.h"

int init_FIR(double *coeffs, int Ncoeffs, FIR_filter *FIR);
int decimate(FIR_filter *FIR, int dec_fac, double *yi, int ni, double *yo, int *no);

void 
error_decim_sac(char *arg)
{
  fprintf(stderr, "Syntax: %s i_sac_file o_sac_file dec_facs... \n", arg);
  exit(1);
}

int 
main(int argc, char *argv[])
{
  int    i,dec_fac,ierror=1;
  char   isacfile[128],osacfile[128];
  FIR_filter FIR2,FIR4,FIR5,*FIR;
  double *y ;
  sachdr hdr ;
  
  if (argc < 4)
    error_decim_sac(argv[0]);

  strcpy(isacfile,argv[1]);
  strcpy(osacfile,argv[2]);

  hdr_alloc(&hdr) ;  
  rhdrsac(isacfile, &hdr, &ierror) ;
  y = double_alloc(hdr.npts) ;
  rdatsac(isacfile, &hdr, y, &ierror) ;
  init_FIR(dec2,Ndec2,&FIR2) ;
  init_FIR(dec4,Ndec4,&FIR4) ;
  init_FIR(dec5,Ndec5,&FIR5) ;

  for(i=3;i<argc;i++)
    {
      if (sscanf(argv[i],"%d",&dec_fac)!=1)
	error_decim_sac(argv[0]);
      if (dec_fac==2)
	FIR = &FIR2 ;
      else if(dec_fac==4)
	FIR = &FIR4 ;
      else if(dec_fac==5)
	FIR = &FIR5 ;
      else
	{
	  fprintf(stderr,"ERROR Incorrect decimation factor ");
	  fprintf(stderr,"(can be 2, 4 or 5)\n");
	  exit(1);
	}
      decimate(FIR,dec_fac,y,hdr.npts,y,&hdr.npts);
      hdr.delta *= (float)dec_fac;
    }
  whdrsac(osacfile,&hdr)  ;
  wdatsac(osacfile,&hdr,y);
  free((void*)FIR2.coeffs);
  free((void*)FIR4.coeffs);
  free((void*)FIR5.coeffs);
  free((void*)y);
  return 0;
}

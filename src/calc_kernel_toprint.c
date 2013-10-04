calc_kernel(str_quake_params *eq,structopt *opt,sachdr *hd_synth,int nsac,char *itype,
			int nd,double *dv,double *tv, double ***G,FILE *o_log)
{
  int   i,j,ns,jsac,nsects,flag,flag2,ngfcomp=6,ierror=1 ;
  long  int nerr ;
  char  ori  ; 
  float *stlats,*stlons,*dists,*azs,*bazs,*xdegs ;
  double **GFs,*S,*TH,*PH,*x_conv   ;
  double gcarc,Ptt,twp_beg,twp_end, *ref_vm=NULL ;
  double *b1,*b2,*a1,*a2,gain,dt=1.; 
  sachdr hdr; 
  make_chan_list(hd_synth,nsac,&stlats,&stlons);
  /* Memory Allocations */
  if (opt->ref_flag)
    {
      ref_vm = double_alloc(NM);
      for(i=0;i<NM;i++) 
		ref_vm[i] = eq->vm[1][i] ;
    }
  dists = float_calloc(nsac) ; 
  azs   = float_calloc(nsac) ;
  bazs  = float_calloc(nsac) ;
  xdegs = float_calloc(nsac) ;
  GFs   = double_alloc2(10,__LEN_SIG__) ;/* GFs: Rrr, Rtt, Rpp, Rrt  */
  S     = double_alloc(__LEN_SIG__)     ;/*    Vertical components   */
  TH    = double_alloc(__LEN_SIG__)     ;/*    Radial components     */
  PH    = double_alloc(__LEN_SIG__)     ;/*    Transverse components */
  x_conv   = double_alloc(__LEN_SIG__)  ;
  hdr_init(&hdr) ; /* SAC header allocation       */
  nsects = (eq->flow > 0.)? eq->filtorder : eq->filtorder/2 ;
  b1 = double_alloc(nsects) ; 
  b2 = double_alloc(nsects) ;
  a1 = double_alloc(nsects) ; 
  a2 = double_alloc(nsects) ;
  /* Distance, azimuth, back-azimuth, etc          */
  distaz(eq->evla, eq->evlo, stlats, stlons, nsac, dists, azs, bazs, xdegs, &nerr) ;
  flag = 0 ;  
  for(i=0;i<ngfcomp;i++) /* Main loop */
    {
      ns   = 0 ; /* Channel counter */
      for(j=0;j<ngfcomp;j++) /* Inititializing the MT components */
		eq->vm[1][j] = 0. ;
      eq->vm[1][i]   = 1. ;
	  for(jsac=0;jsac<nsac;jsac++)
		{
		  gcarc = (double)hd_synth[jsac].gcarc;
		  trav_time(gcarc,tv,dv,nd,&Ptt,&ierror) ;
		  wp_time_window(gcarc,eq->wp_win4,&twp_beg,&twp_end) ; /* Data Time Window  */
		  flag2 = 0;
		  ori = hd_synth[ns].kcmpnm[2];
		  
		  if ( ori == 'Z' )
			  fast_synth_only_Z_sub(azs[jsac],bazs[jsac],xdegs[jsac], tv,dv,nd,eq,&hdr,GFs,S);
		  else if ( ori == 'N' || ori == 'E' || ori == '1' || ori == '2' ) 
			{
			  fast_synth_only_Hs_sub(azs[jsac],bazs[jsac],xdegs[jsac],tv,dv,nd,eq,&hdr,GFs,TH,PH);
			  rotate_traces(TH, PH, hd_synth[jsac].baz-hd_synth[jsac].cmpaz, hdr.npts, S) ; /*Rotating TH, PH to H*/
			}
		  else
			continue;
		  
		  conv_by_stf(eq->ts,eq->hd,itype,&hdr,S,x_conv) ;/* Perform convolution */
		  if (flag == 0) /* Set the butterworth sos (dt must be the same for all stations)   */
			{
			  flag = 1 ; 
			  dt = (double)hdr.delta;
			  if (eq->flow>0.)
				bpbu2sos(eq->flow,eq->fhigh,dt,eq->filtorder,&gain,b1,b2,a1,a2);
			  else
				lpbu2sos(eq->fhigh,dt,eq->filtorder,&gain,b1,b2,a1,a2);		  		      
			}
		  else if (dt != (double)hdr.delta) /* Check the sampling frequency (must be uniform) */
			{
			  fprintf(stderr, "*** ERROR: non uniform samp. period between sac files (%s)\n",hd_synth[jsac].kstnm);
			  fprintf(stderr,"    -> Exiting\n") ;
			  fflush(stderr);
			  exit(1);
			}	  
		  filter_with_sos(gain,b1,b2,a1,a2,nsects,x_conv,hdr.npts) ; /* Apply sos */
		  flag2 = fill_kernel_G(&hdr,&(hd_synth[jsac]),Ptt,twp_beg,twp_end,x_conv,G[jsac][i],opt,o_log);
		  if (flag2)
			{
			  fprintf(stderr,"*** ERROR: Incomplete green function for %s\n",hd_synth[jsac].kstnm) ;
			  fprintf(stderr,"    -> Exiting\n") ;
			  fflush(stderr);
			  exit(1);
			}
		  hd_synth[jsac].b  = hdr.b + hd_synth[jsac].o - hdr.o ;
		  ns++;
		}
	  if (nsac!=ns)
		{	
		  fprintf(stderr,"\n*** ERROR: Kernel G is incomplete (%d vs %d)\n",nsac,ns);  
		  fprintf(stderr,"    -> Exiting\n") ;
		  fflush(stderr);
		  exit(1);	  
		}
	}
  
  /* Free Memory */
  if (opt->ref_flag)
    {
      for(i=0;i<NM;i++) 
		eq->vm[1][i] = ref_vm[i] ;
      free((void*)ref_vm);
    }
  free((void*)stlats) ;
  free((void*)stlons) ;
  free((void*)dists)  ;
  free((void*)xdegs)  ;
  free((void*)bazs)   ;
  free((void*)azs)    ;
  for(i=0;i<10;i++)
    free((void*)GFs[i]) ;
  free((void**)GFs)     ;
  free((void*)x_conv)   ;
  free((void*)S)    ;
  free((void*)TH)   ;
  free((void*)PH)   ;
  free((void*)a1)   ;
  free((void*)a2)   ;
  free((void*)b1)   ;
  free((void*)b2)   ;
}

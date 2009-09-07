#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "read_i_files.h"
#include "proto_alloc.h"
#include "rwtextfiles.h"


void 
add_slash(char *c)
{
  if (c[strlen(c)-1] != '/')
    strcat(c,"/");
}


/**********************************************************************/
/*                         get_cmtf(eq, flag)                         */
/**********************************************************************/
/* Read cmtfile (shared routine).                                     */
/* Input : eq->cmtfile                                                */
/*         flag : if flag = 0 : read  pde                             */
/*                if flag = 1 :     +|evname                          */
/*                                   |t-shift                         */
/*                                   |half-dur                        */
/*                if flag = 2 :     +|ref lat                         */
/*                                   |ref lon                         */
/*                                   |ref dep                         */
/*                                   |ref mom                         */
/* Output : eq : parameters in a structure                            */
/*          eq->vm is a pointer on 2 arrays (allocated before calling */
/*                 this routine if flag=2)                            */
/*          > eq->vm[0] moment tensor elements to be determined later */
/*                      by inverting Wphase.                          */
/*          > eq->vm[1] moment tensor elements of the reference sol.  */
/*                      loaded from cmtfile by this routine if flag=2 */
void 
get_cmtf(eq, flag)
     int   flag            ;
     str_quake_params *eq  ; 
{
  int    nl, nb, nb2, tmp1 ;
  char   *line             ;
  double tmp2              ;
  FILE   *cmtfile          ;
  
  /* Allocate memory    */
  line = char_alloc(LSIZE);
  
  /* Openning file      */
  cmtfile = openfile_rt(eq->cmtfile,&nl);

  /* Checking variables */
  if (nl < 1) {
    fprintf(stderr,"ERROR : incomplete cmtfile : %s\n", eq->cmtfile) ;
    fclose(cmtfile) ;
    exit(1) ; }  
  else if (flag!=0 && flag!=1 && flag!=2){
    fprintf(stderr,"ERROR reading cmtfile : %s (incorrect flag) \n", eq->cmtfile) ;
    fclose(cmtfile) ;
    exit(1) ; }  


  /* Read lines */
  nl = 0;
  while( fgets(line,LSIZE,cmtfile) != NULL )
    {
      nb = nb_blank(line) ;
      if (nb != strlen(line))
	nb2 = nb;
      else
	continue ;
      if (strncmp(&line[nb2],"PDE",3)==0 || strncmp(&line[nb2],"MLI",3)==0 || 
	  strncmp(&line[nb2],"ISC",3)==0 || strncmp(&line[nb2],"REB",3)==0 || 
	  strncmp(&line[nb2],"SWE",3)==0 || strncmp(&line[nb2],"HSW",3)==0 || 
	  strncmp(&line[nb2],"JMA",3)==0)
	{
	  strcpy(eq->pdeline,line);
	  if (strncmp(&line[nb2],"PDEW",4) || strncmp(&line[nb2],"SWEQ",4))
	    nb2+= 4 ;
	  else if (line[nb2+3] == ' ')
	    nb2+= 3 ;
	  else {
	    fprintf(stderr,"ERROR reading PDE in %s (unknown format) \n", eq->cmtfile) ;
	    fclose(cmtfile) ;
	    exit(1)         ; }  
	  tmp1 = sscanf(&line[nb2], "%i %i %i %i %i %lf %lf %lf %lf", 
			&eq->ot_ye, &eq->ot_mo, &eq->ot_dm, 
			&eq->ot_ho, &eq->ot_mi, &tmp2     , 
			&eq->pde_evla , &eq->pde_evlo , &eq->pde_evdp ) ;
	  check_scan(9, tmp1   , eq->cmtfile  , cmtfile);
	  eq->ot_se = (int) tmp2;
	  eq->ot_ms = (int) ((tmp2 - eq->ot_se)*1000 + 0.5);
	  nl++;	
	}
      else if (flag == 1 || flag == 2)
	{
	  if (strncmp(&line[nb2],"event name",10)==0){
	    nb2+= 10+nb_blank(line+10)+1 ;
	    tmp1 = sscanf (&line[nb2], "%s", eq->evid) ;
	    check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
	    nl++;   }
	  else if (strncmp(&line[nb2],"time shift",10)==0){
	    nb2+= 10+nb_blank(line+10)+1 ;
	    tmp1 = sscanf (&line[nb2], "%lf", &eq->ts) ;
	    check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
	    nl++;	}
	  else if (strncmp(&line[nb2],"half duration",13)==0){
	    nb2+= 13+nb_blank(line+13)+1 ;
	    tmp1 = sscanf (&line[nb2], "%lf", &eq->hd) ;
	    check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
	    nl++;	}
	  else if (flag == 2)
	    {
	      if (strncmp(&line[nb2],"latitude",8)==0){
		nb2+= 8+nb_blank(line+8)+1 ;
		tmp1 = sscanf (&line[nb2], "%lf", &eq->evla) ;
		check_scan(1, tmp1, eq->cmtfile, cmtfile)    ;
		nl++;    }     
	      else if (strncmp(&line[nb2],"longitude",9)==0){
		nb2+= 9+nb_blank(line+9)+1 ;
		tmp1 = sscanf (&line[nb2], "%lf", &eq->evlo) ;
		check_scan(1, tmp1, eq->cmtfile, cmtfile)    ;
		nl++;    }     
	      else if (strncmp(&line[nb2],"depth",5)==0){
		nb2+= 5+nb_blank(line+5)+1 ;
		tmp1 = sscanf (&line[nb2], "%lf", &eq->evdp) ;
		check_scan(1, tmp1, eq->cmtfile, cmtfile)    ;
		nl++;    }
	      else if (strncmp(&line[nb2],"Mrr",3)==0){
		nb2    += 3+nb_blank(line+3)+1 ;
		tmp1    = sscanf (&line[nb2], "%lf", &eq->vm[1][0]) ;
		eq->vm[1][0] /= (double)POW;
		check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
		nl++;    }
	      else if (strncmp(&line[nb2],"Mtt",3)==0){
		nb2    += 3+nb_blank(line+3)+1 ;
		tmp1    = sscanf (&line[nb2], "%lf", &eq->vm[1][1]) ;
		eq->vm[1][1] /= (double)POW;
		check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
		nl++;    }
	      else if (strncmp(&line[nb2],"Mpp",3)==0){
		nb2    += 3+nb_blank(line+3)+1 ;
		tmp1    = sscanf (&line[nb2], "%lf", &eq->vm[1][2]) ;
		eq->vm[1][2] /= (double)POW;
		check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
		nl++;    }
	      else if (strncmp(&line[nb2],"Mrt",3)==0){
		nb2    += 3+nb_blank(line+3)+1 ;
		tmp1    = sscanf (&line[nb2], "%lf", &eq->vm[1][3]) ;
		eq->vm[1][3] /= (double)POW;
		check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
		nl++;    }
	      else if (strncmp(&line[nb2],"Mrp",3)==0){
		nb2    += 3+nb_blank(line+3)+1 ;
		tmp1    = sscanf (&line[nb2], "%lf", &eq->vm[1][4]) ;
		eq->vm[1][4] /= (double)POW;	  
		check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
		nl++;    }      
	      else if (strncmp(&line[nb2],"Mtp",3)==0){
		nb2    += 3+nb_blank(line+3)+1 ;
		tmp1    = sscanf (&line[nb2], "%lf", &eq->vm[1][5]) ;
		eq->vm[1][5] /= (double)POW;
		check_scan(1, tmp1, eq->cmtfile, cmtfile) ;
		nl++;    }
	    }
	}
    }

  fclose(cmtfile);

  if (flag == 1){
    strcpy(eq->cmtfile, "NO REF. SOL. USED (PDE only)") ;
    eq->evla = eq->evlo = eq->evdp = 0.0     ;}  

  /* Memory Freeing */
  free((void *)line) ;  
  
  /* Check if nl is correct */
  if ((flag == 1 && nl < 4) || (flag == 2 && nl < 13)){
    fprintf(stderr,"ERROR: reading cmtfile %s (field(s) are missing)\n",eq->cmtfile) ;
    exit(1);}
  else if ((flag == 1 && nl > 4) || (flag == 2 && nl > 13)){
    fprintf(stderr,"ERROR: reading cmtfile %s (field(s) redundancy)\n",eq->cmtfile) ;
    exit(1);}
  else if (!flag && nl < 1) {
    fprintf(stderr,"ERROR: PDE missing in cmtfile : %s \n",eq->cmtfile) ;
    exit(1);}
  else if (!flag && nl > 1) {
    fprintf(stderr,"ERROR: PDE redundancy in cmtfile : %s \n",eq->cmtfile) ;
    exit(1);}
}


void 
decode_wp_win(buffer, wp_win4)
     double *wp_win4 ;     
     char *buffer    ; 
{
  int    np;
  double f1, f2, f3, f4;

  np = sscanf(buffer, "%lf %lf %lf %lf", &f1, &f2, &f3, &f4);
  wp_win4[0] =   0.;
  wp_win4[1] =   0.;
  wp_win4[3] =   0.;
  wp_win4[3] = 180.;
  if(np < 1  || np > 4)
	{
	  fprintf(stderr, "Error reading the Wp time window parameters (i_wpinversion)\n");
	  exit(1);
	}
  else if(np == 1) 
	{
	  wp_win4[0] =   0.;
	  wp_win4[1] =   f1;
	  wp_win4[2] =   0.;
	  wp_win4[3] = 180.;
	}
  else if(np == 2) 
	{
	  wp_win4[0] =   f1;
	  wp_win4[1] =   f2;
	  wp_win4[2] =   0.;
	  wp_win4[3] = 180.;
	}
  else if(np == 3) 
	{
	  wp_win4[0] =   f1;
	  wp_win4[1] =   f2;
	  wp_win4[2] =   f3;
	  wp_win4[3] = 180.;
	}
  else 
	{
	  wp_win4[0] =   f1;
	  wp_win4[1] =   f2;
	  wp_win4[2] =   f3;
	  wp_win4[3] =   f4;
	}
}


/************************************************/
/*       get_i_master(file, keys, n, eq)        */
/************************************************/
/* Read i_masterfile (shared routine).          */
/* Input : >file : i_master filename            */
/*         >keys : list of variables keyword to */
/*                 read.                        */
/*         >nb of keywords                      */
/*                                              */
/* Output : eq : parameters in a structure      */
void 
get_i_master(file, keys, n, eq)
     int   n ;
     char  *file, **keys  ; 
     str_quake_params *eq ;
{
  int  i, nb, nb2, tmp, nl;
  char *line,*buf;
  FILE *i_file ;

  /* Allocate memory          */
  line  = char_alloc(LSIZE) ;
  buf   = char_alloc(LSIZE) ;

  /* Openning file            */
  i_file = openfile_rt(file,&nl) ;

  /* check nb of lines */
  if (n > nl) {
    fprintf(stderr,"ERROR : incomplete i_master file : %s\n", file) ;
    fclose(i_file) ;
    exit(1) ; }

  /* Initialize GFDIR  */
  /* strcpy(eq->gf_dir,"./GF/");  */

  /* Initialize DMIN and DMAX */
  nl = 0;
  for (i=0 ; i<n ; i++)
    {
      fflush(stdout);
      if (strncmp("DMAX",keys[i],strlen(keys[i]))==0)
	nl++; 
      else if (strncmp("DMIN",keys[i],strlen(keys[i]))==0){
	eq->dmin = 0. ;
	nl++; }
    }

  /* Read lines */
  while( fgets(line,LSIZE,i_file) != NULL )
    {
      nb = nb_blank(line) ;
      if (nb != strlen(line))
	nb2 = nb;
      else
	continue ;
      for (i=0 ; i<n ; i++) /* Keyword search */
	if (strncmp(&line[nb2],keys[i],strlen(keys[i]))==0) 
	  {
	    if (strncmp(&line[nb2],"EVNAME",6)==0) {
	      nb2+= 6+nb_blank(&line[nb2+6])+1          ;
	      tmp = sscanf (&line[nb2], "%s", eq->evnm) ;
	      check_scan(1, tmp, file, i_file)          ;
	       nl++;   }
	    else if (strncmp(&line[nb2],"SEED",4)==0) {
	      nb2+= 4+nb_blank(&line[nb2+4])+1          ;
	      tmp = sscanf (&line[nb2], "%s", eq->seed) ;
	      check_scan(1, tmp, file, i_file)          ;
	       nl++;   }
	    else if (strncmp(&line[nb2],"CMTFILE",7)==0){
	      nb2+= 7+nb_blank(&line[nb2+7])+1               ;
	      if (strlen(eq->cmtfile) == 0){
		tmp = sscanf (&line[nb2], "%s", eq->cmtfile) ;
		check_scan(1, tmp, file, i_file)             ;}
	      nl++;   }
	    else if (strncmp(&line[nb2],"DMIN",4)==0){
	      nb2+= 4+nb_blank(&line[nb2+4])+1            ;
	      tmp = sscanf (&line[nb2], "%lf", &eq->dmin) ;
	      check_scan(1, tmp, file, i_file)            ;}
	    else if (strncmp(&line[nb2],"DMAX",4)==0){
	      nb2+= 4+nb_blank(&line[nb2+4])+1              ;
	      tmp = sscanf ((&line[nb2]), "%lf", &eq->dmax) ;
	      check_scan(1, tmp, file, i_file)              ;}
	    else if (strncmp(&line[nb2],"filt_order",10)==0){
	      nb2+= 10+nb_blank(&line[nb2+10])+1                ;
	      tmp = sscanf ((&line[nb2]), "%d", &eq->filtorder) ;
	      check_scan(1, tmp, file, i_file)                  ;
	      nl++ ;  }
	    else if (strncmp(&line[nb2],"filt_pass",9)==0){
	      nb2+= 9+nb_blank(&line[nb2+9])+1                  ;
	      tmp = sscanf ((&line[nb2]), "%d", &eq->filtnpass) ;
	      check_scan(1, tmp, file, i_file)                  ;
	      nl++ ;  }
	    else if (strncmp(&line[nb2],"filt_cf1",8)==0){
	      nb2+= 8+nb_blank(&line[nb2+8])+1              ;
	      tmp = sscanf ((&line[nb2]), "%lf", &eq->flow) ;
	      check_scan(1, tmp, file, i_file)              ;
	      nl++ ;  }
	    else if (strncmp(&line[nb2],"filt_cf2",8)==0){
	      nb2+= 8+nb_blank(&line[nb2+8])+1               ;
	      tmp = sscanf ((&line[nb2]), "%lf", &eq->fhigh) ;
	      check_scan(1, tmp, file, i_file)               ;
	      nl++ ;  }
	    else if (strncmp(&line[nb2],"IDEC_2",6)==0){
	      nb2+= 6+nb_blank(&line[nb2+6])+1           ;
	      tmp = sscanf (&line[nb2], "%lf %lf %lf", 
			    &eq->idtr, &eq->preevent, &eq->fend) ;
	      check_scan(3, tmp, file, i_file)           ;
	      nl++;   }
	    else if (strncmp(&line[nb2],"IDEC_3",6)==0) {
	      nb2+= 6+nb_blank(&line[nb2+6])+1            ;
	      tmp = sscanf (&line[nb2], "%lf %lf %d %lf", 
			      &eq->fl, &eq->fh, 
			    &eq->nf, &eq->tol)            ;
	      check_scan(4, tmp, file, i_file)            ;
	      nl++;   }
	    else if (strncmp(&line[nb2],"GFDIR",5)==0) {
	      nb2+= 5+nb_blank(&line[nb2+5])+1            ;
	      tmp = sscanf (&line[nb2], "%s", eq->gf_dir) ;
	      check_scan(1, tmp, file, i_file)            ;
	      add_slash(eq->gf_dir);
	      nl++;   }
	    else if (strncmp(&line[nb2],"WP_WIN",6)==0) {
	      nb2+= 6+nb_blank(&line[nb2+6])+1       ;
	      decode_wp_win(&line[nb2], eq->wp_win4) ;
	      nl++;   }
	  }
    }

  fclose(i_file)     ;

  /* Check nl */
  if (n > nl)  {
    fprintf(stderr,"ERROR 2 : incomplete i_master file : %s\n", file) ;
    exit(1) ; }
  else if (n < nl) {
    fprintf(stderr,"ERROR : redundancy in i_master file : %s\n", file) ;
    exit(1) ; }
  /* Memory Freeing */
  free((void *)buf)  ;
  free((void *)line) ;
}


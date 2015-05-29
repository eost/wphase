/***************************************************************************
*
*                     W phase source inversion package              
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Universite de Strasbourg / CNRS 
*                                  April 2013
*
*    Neither the name of the California Institute of Technology (Caltech) 
*    nor the names of its contributors may be used to endorse or promote 
*    products derived from this software without specific prior written 
*    permission
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
****************************************************************************/

/*      Recursive resitution and filtering         */

#define _GNU_SOURCE
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Subroutines headers */
#include "proto_alloc.h"
#include "rwsacs.h"     
#include "rwtextfiles.h" 
#include "read_i_files.h"
#include "sort_tree.h"
#include "butterworth.h"

/* internal functions */
void makeid(char *id, sachdr *hdr) ;
int  findid(char *id, char **ids, int n) ;
void trapz(double *x_in, int npts, double dt, double *x_int) ;
void get_id_dec(char *file, int *n, char ***id, double **a1, double **a2, double **a3) ;
void get_filt_params(char *file, str_quake_params *eq) ;


int main(int argc, char *argv[])
{
    int     i, n ,tmp, count, nsects, ierror = 1;
    int     beg, end, MAX, nsac;
    double  *c1, *c2, *c3, *b1,*b2,*a1,*a2, scale, gain ;
    double  *x_in, *x_int, dt=1. ;
    char    **ids, *x_int_fil ;
    char    *decfile, *i_master, *i_sacfile ;
    char    *id, *o_sacfile ;
    FILE    *i_sacf, *o_sacf;
    sachdr  hdr;
    str_quake_params eq ;
    struct tree *datfil ;

    /* Input params */
    if(argc<5)
    {
        fprintf (stderr, "*** ERROR (input parameters : 4 params needed (%d given))\n", argc-1) ;
        fprintf (stderr, "Syntax : %s i_coeffs(input) i_master i_scr_data_list(input) o_scr_data_list(output)\n", argv[0]) ;
        exit (1) ; 
    }
    decfile   = char_alloc(FSIZE) ;
    i_sacfile = char_alloc(FSIZE) ;
    o_sacfile = char_alloc(FSIZE) ;
    i_master  = char_alloc(FSIZE) ;
    strcpy(  decfile, argv[1])    ;
    strcpy( i_master, argv[2])    ;
    strcpy(i_sacfile, argv[3])    ; 
    strcpy(o_sacfile, argv[4])    ; 
    get_filt_params(i_master,&eq)           ;
    get_id_dec(decfile,&n,&ids,&c1,&c2,&c3) ;
    /* Allocate data tabs */
    MAX   = (int) __LEN_SIG__ ;
    x_in  = double_alloc(MAX) ;
    x_int = double_alloc(MAX) ;
    /* Allocate tree and sac header */
    hdr_init(&hdr)      ;
    datfil = alloctree();
    datfil->d = NULL ;
    datfil->g = NULL ;
    /* Allocate sac filenames */
    x_int_fil = char_alloc(FSIZE);
    /* Allocate id */
    id = char_alloc(IDSIZE);
    /* Allocate sos */
    nsects = (eq.flow > 0.)? eq.filtorder+1 : eq.filtorder/2+1;
    b1 = double_alloc(nsects);
    b2 = double_alloc(nsects);
    a1 = double_alloc(nsects);
    a2 = double_alloc(nsects);
    /* Open Data File Lists */
    i_sacf = openfile_rt(i_sacfile,&nsac) ;
    o_sacf = openfile_wt(o_sacfile)       ;
    i = 0 ; count = 0 ;
    datfil->occur = 1 ;
    datfil->npts = -12345 ;
    while( (tmp=fscanf(i_sacf,"%s %s %s %s %s %lf %lf %lf %lf %lf\n",
                     datfil->file,datfil->sta,datfil->net,datfil->cmp,
                         datfil->locid,&datfil->stla,&datfil->stlo,&datfil->cmpaz,
                     &datfil->az,&datfil->xdeg)) != EOF )
    {
        if (tmp == 0)
            break;
        else
            check_scan(10,tmp,i_sacfile,i_sacf) ;
        /* Read input sac file */
        rhdrsac(datfil->file,&hdr,&ierror) ;
        if (hdr.npts > MAX) 
            hdr.npts = MAX ;      
        rdatsac(datfil->file,&hdr,x_in,&ierror) ;
        if (hdr.npts<2)
        {
            fprintf(stderr,"Warning (rec_dec_filt) : %s (rejected) npts < 2\n",datfil->file);
            continue;
        }
        /* Set the butterworth sos (samp. rate must be the same for all stations)*/
        if (count==0 && i==0)     
        {
            dt = (double)hdr.delta;
            if (eq.flow>0.)
                bpbu2sos(eq.flow,eq.fhigh,dt,eq.filtorder,&gain,b1+1,b2+1,a1+1,a2+1);
            else
                lpbu2sos(eq.fhigh,dt,eq.filtorder,&gain,b1+1,b2+1,a1+1,a2+1);
        }
        else if (dt!=(double)hdr.delta)
        {
            fprintf(stderr, "ERROR (rec_dec_filt): non uniform samp. period between sac files, file : %s\n",datfil->file);
            fclose(i_sacf);
            exit(1);
        }
        /* Set filenames */
        strcpy(x_int_fil,datfil->file) ;
        strcat(x_int_fil,".id") ;
        strcat(datfil->file,".dec.bp.int") ;

        /* Base line operations */
        if (eq.idtr==1)  /*   detrend and taper, eq.preevent = duration of taper at the beginning */
        {              /*                      eq.fend = duration of taper at the end           */
            dtrd( x_in,hdr.npts) ;                       /* dtrend */
            beg = (int)(eq.preevent/hdr.delta+0.5);
            end = (int)(eq.fend/hdr.delta+0.5);
            taper(x_in,hdr.npts,eq.preevent,eq.fend) ;   /* taper  */
        }
        else if (eq.idtr==2)   /*   shift the base line, eq.preevent=duration over which the  */
        {                    /*                      base line is defined, eq.fend = dummy  */
            beg = 0;
            end = (int)(eq.preevent/hdr.delta)-1;
            rmean(x_in,hdr.npts,&beg,&end);
        }
        /* Trapezoidal numerical integration */
        trapz(x_in,hdr.npts,(double)hdr.delta,x_int);
        /* Write header values and data in the sac input file */
        wsac(x_int_fil,&hdr,x_int) ;
        makeid(id,&hdr);
        i = findid(id,ids,n);
        if (i == -1)
            continue ;
        /* Set the first line of the sos (Deconvolution) */
        b1[0] = c2[i]/c3[i] ;
        b2[0] = c1[i]/c3[i] ;
        a1[0] = -1.0 ;
        a2[0] =  0.0 ;
        scale = gain*c3[i] ;
        /* Applying sos */
        filter_with_sos(scale,b1,b2,a1,a2,nsects,x_in,hdr.npts) ; /* Apply sos */
        /* Trapezoidal numerical integration */
        trapz(x_in,hdr.npts,(double)hdr.delta,x_int) ;
        /* Integration again */
        trapz(x_int,hdr.npts,(double)hdr.delta,x_in) ;
        /* Write output sac  */
        wsac(datfil->file,&hdr,x_in) ;
        savetree(datfil,o_sacf,&eq.dmin,&eq.dmax) ;
        count++;
    }
    fclose(i_sacf) ;
    fclose(o_sacf) ;
    printf("%d channels remain after deconvolution and filtering (%d rejected)\n",count,nsac-count) ;

    /* Memory Freeing */
    free((void *)decfile);
    free((void *)i_master);
    free((void *)i_sacfile);
    free((void *)o_sacfile);
    free((void *)x_in);
    free((void *)x_int);
    free((void *)x_int_fil);
    free((void *)a1);
    free((void *)a2);
    free((void *)b1);
    free((void *)b2);
    free((void *)id);
    for(i=0; i<n; i++)
        free((void *)ids[i]);
    free((void**)ids);
    free((void *)c1);
    free((void *)c2);
    free((void *)c3);
    freetree(datfil);
    return 0;
}


/************************************************/
/*           get_filt_params(file, eq)          */
/************************************************/
/*  > Read input file for recursive filtering   */
void get_filt_params(char *file, str_quake_params *eq)
{
    int  i = 0  ;
    char **keys ;

    keys = char_alloc2(6, 16) ;
    strcpy(keys[i++],"filt_order") ;
    strcpy(keys[i++],"filt_cf1")   ;
    strcpy(keys[i++],"filt_cf2")   ;
    strcpy(keys[i++],"IDEC_2")     ;
    strcpy(keys[i++],"DMIN")     ;
    strcpy(keys[i++],"DMAX")     ;
    get_i_master(file,keys,6,eq) ;
    for(i=0 ; i<6 ; i++)
        free((void*)keys[i]) ;
    free((void**) keys )   ;
}

void makeid(char *id, sachdr *hdr)
{
    int tmp;
    strcpy(id, hdr->knetwk) ;
    id[nbchar(hdr->knetwk)] = '\0';
    strcat(id,"_") ;
    strncat(id, hdr->kstnm,nbchar(hdr->kstnm)) ;
    strcat(id,"_") ;
    tmp = nbchar(hdr->khole) ;
    if ( tmp == 0 )
        strcat(id, "--") ;
    else
        strncat(id,hdr->khole,tmp);
    strcat(id,"_") ;
    strncat(id, hdr->kcmpnm, nbchar(hdr->kcmpnm)) ;
    strcat(id,"");
}

int findid(char *id, char **ids, int n)
{
    int i;
    /* Find id */
    for (i=0; i<n; i++)
    {
        if (strcmp(id,ids[i])==0)
            return i;
    }
    fprintf (stderr, "WARNING (rec_dec_filt): %s not found in id list (rejected)\n",id) ;
    return -1;
}


void trapz(double *x_in, int npts, double dt, double *x_int)
{
    int i;
  
    //x_int = (double *) malloc (npts * sizeof(double));
    x_int[0] = 0.0 ;
    for (i=1 ; i < npts ; i++)
        x_int[i] = x_int[i-1] + (x_in[i-1] + x_in[i])*dt/2.0 ;
}


void get_id_dec(char *file, int *n, char ***id, double **a1, double **a2, double **a3)
{
    int   i, tmp ;
    FILE *decin ;  

    /* Open Data File List */
    decin = openfile_rt(file,n);  
    /* Allocates memory */
    (*id) = char_alloc2((*n),IDSIZE);
    (*a1) = double_alloc((*n));
    (*a2) = double_alloc((*n));
    (*a3) = double_alloc((*n));
  
    for(i=0 ; i<(*n) ; i++)
    {
        tmp = fscanf (decin, "%s %lf %lf %lf",(*id)[i],(*a1)+i,(*a2)+i,(*a3)+i);
        check_scan(4,tmp,file,decin);
    }
    fclose(decin);
}

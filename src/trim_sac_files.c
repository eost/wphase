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

/*      Trim sac files       */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>
#include <locale.h>

/* Subroutines headers */
#include "proto_alloc.h"    /* proto_alloc.c    */
#include "rwsacs.h"      /* rwsacs.c      */
#include "rwtextfiles.h" /* rwtextfiles.c */
#include "sort_tree.h"   /* sort_tree.h   */
#include "travel_times.h"
#include "read_i_files.h"

#ifndef SAMPLEPERIOD
#define SAMPLEPERIOD 1.
#endif /* not SAMPLEPERIOD */


/* External routines */
void distaz(double    cmt_lat,  double    cmt_lon,  float*    stlats,  float*    stlons, 
            int       nstat,    float*    dists,    float*    azs,     float*    bazs,   
            float*    xdegs,    long int* nerr);  


/* Internal routines */
void get_params(int argc, char **argv, int *un, char **i_sacs, char **o_sacs, 
                str_quake_params *eq) ;
void date2epoch(int year, int mm, int dd, int hour, int min, int sec, int msec, double *epoch);
void delta_t(int y1, int j1, int h1, int m1, int s1, int ms1,
             int y0, int j0, int h0, int m0, int s0, int ms0, double *tdiff) ;
void ymd2jul(int yyyy,int mm,int dd, int *jul) ;
          
int main(int argc, char *argv[])
{
    int    nc, nh=NDEPTHS, nd=NDISTAS ; 
    int    tmp, ierror = 1, tterr = 0 ; 
    int    ot_jday, sampstart, flag   ;
    long   int nerr                   ;
    float  dist, az, baz, xdeg        ;
    double xdegd, P_tt, tdiff, otime  ;
    double *x_in, *x_out, *dv, *tv    ;
    char   *i_sacs, *fil,*o_sacs      ;
    FILE   *i_sacf, *o_sacf   ;
    str_quake_params   eq     ;
    sachdr hdr                ;
    struct tree   *root, *mod ;
  
    /* OPTIONS */
    int un ;
    /* Input params */
    get_params(argc,argv,&un,&i_sacs,&o_sacs,&eq);
    /* Allocate memory */
    fil    = char_alloc(FSIZE) ;
    x_in   = double_alloc((int)__LEN_SIG__);
    dv     = double_alloc(nd);
    tv     = double_alloc(nd);
    hdr_init(&hdr)    ;  
    root = alloctree();
    mod  = alloctree();
    /* Set travel time table */
    trav_time_init(nh,nd,eq.pde_evdp,dv,tv,&tterr) ;
    /* Main loop */
    i_sacf = openfile_rt(i_sacs,&nc) ;
    nc = 0   ;
    flag = 0 ;
    while( (tmp=fscanf (i_sacf, "%s", fil)) != EOF )
    {
        check_scan(1,tmp,i_sacs,i_sacf);
        /* Read sac file */
        rhdrsac(fil,&hdr,&ierror);
        if (hdr.npts > (int)__LEN_SIG__)  
            hdr.npts = (int)__LEN_SIG__ ;
        rdatsac(fil,&hdr,x_in,&ierror) ;
        /* Set epicentral distances, azimuth, backazimuth */
        dist = 0. ;
        az   = 0. ;
        baz  = 0. ;
        xdeg = 0. ;
        distaz(eq.pde_evla,eq.pde_evlo,&hdr.stla,&hdr.stlo,1,&dist,&az,&baz,&xdeg,&nerr) ;
        /* Set travel time */
        xdegd = (double)xdeg ;
        trav_time(xdegd,tv,dv,nd,&P_tt,&tterr) ;
        /* Set event origin time in sac header variable 'o' (relative to the reference time) */
        ymd2jul(eq.ot_ye,eq.ot_mo,eq.ot_dm,&ot_jday) ;
        delta_t(eq.ot_ye,ot_jday,eq.ot_ho,eq.ot_mi,eq.ot_se,eq.ot_ms, hdr.nzyear,hdr.nzjday,hdr.nzhour,hdr.nzmin,hdr.nzsec,hdr.nzmsec,&tdiff) ;
        otime = tdiff ;
        hdr.t[0] = (float)(P_tt+tdiff) ; /* P arrival         */

        /* Windowing -- Screening by distance */
        tdiff += P_tt-(double)hdr.b-eq.preevent ;   /* time for the 1st sample */
        sampstart = (int)(tdiff/((double)hdr.delta) + 0.5)                   ;
        if ((int)(hdr.delta*1000.+0.5) != (int)((float)SAMPLEPERIOD*1000.+0.5)) /* Sampling period check    */
        {
            fprintf(stderr, "WARNING: incorrect samp. period between sac files\n") ;
            fprintf(stderr, "     ...file : %s with dt = %e is rejected\n", fil, hdr.delta)  ;
            flag++   ;
            continue ;
        }
        else if (sampstart >= 0 && sampstart < hdr.npts)
        {
            /* Set new sac header variables */
            hdr.delta = (float) SAMPLEPERIOD  ;
            hdr.o     = (float) otime         ;   /* event origin time */
            hdr.npts  = hdr.npts - sampstart ;                 /* nb of samples                         */
            hdr.b     = hdr.b + ((float)sampstart)*hdr.delta ; /* shift of the first sample (corrected) */
            hdr.e     = hdr.b + (hdr.npts-1) * hdr.delta    ; 
            hdr.dist  = dist ; /* epicentral distance (km)                      */
            hdr.gcarc = xdeg ; /* station to event great circle arc length(deg) */
            hdr.az    = az   ; /* event to station azimuth (deg)                */
            hdr.baz   = baz  ; /* station to event azimuth (deg)                */
            hdr.evla  = (float) eq.pde_evla ;
            hdr.evlo  = (float) eq.pde_evlo ;
            hdr.evdp  = (float) eq.pde_evdp ;
            /* Write trimmed sac file */
            if (hdr.npts > (int)__LEN_SIG__)
                fprintf(stderr,"Warning : traces cut to %d samples.\n",(int)__LEN_SIG__);
            x_out = &x_in[sampstart] ; /* shift of the first sample */      
            wsac(fil,&hdr,x_out)     ;
            /* Sort sac files */
            if (nc==0)
                splithdr(fil,&hdr,root);
            else
            {
                splithdr(fil,&hdr,mod) ;
                build (root,mod,un)     ;
            }
            nc++ ;
        }
        else 
        {
            fprintf(stderr,"Warning: (trim_sac_files) rejected incomplete file : %s (%d -- %d)\n", fil, sampstart, hdr.npts);
        }
    }
    fclose(i_sacf);
    if (flag > 1)
    {
        fprintf(stderr, "WARNING: %d files have been rejected because of incorrect sampling period\n",flag);
        fprintf(stderr, "    ...: if to much files are rejected, please screen or decimate data files manually\n");
    }
    o_sacf = openfile_wt(o_sacs);
    savetree(root, o_sacf, &eq.dmin, &eq.dmax);
    fclose(o_sacf);
    /* Memory Freeing */
    free((void*)fil)    ;
    free((void*)i_sacs) ;
    free((void*)o_sacs) ;
    free((void*)x_in)   ;
    free((void*) dv)    ;
    free((void*) tv)    ;
    freetree(root);
    freetree(mod);
    return 0;
}


void get_params(int argc, char **argv, int *un, char **i_sacs, char **o_sacs, 
           str_quake_params *eq)
{
    int  i              ;
    char **keys, *i_tmp ;
  
    if ( argc < 4 ) 
    {
        fprintf (stderr, "*** ERROR (minimum of 3 params needed (%d given)) ***",argc-1)        ;
        fprintf (stderr, "Syntax : %s i_master(in) i_sac_list(in) o_sac_list(out) [-u(allow only one network per channel)  -a(all channels)]\n", argv[0]) ;
        exit(1) ; 
    }

    *un = 0 ; 
    if (!strncmp(argv[argc-1],"-u",2))
        *un = 1  ;
    else if (!strncmp(argv[argc-1],"-a",2))
        *un = -1 ;

    i_tmp     = char_alloc(FSIZE) ;  
    (*i_sacs) = char_alloc(FSIZE) ;
    (*o_sacs) = char_alloc(FSIZE) ;
    strcpy(   i_tmp , argv[1]) ;
    strcpy((*i_sacs), argv[2]) ;
    strcpy((*o_sacs), argv[3]) ;
    eq->cmtfile[0] = '\0';

    i = 0;
    keys = char_alloc2(4, 16)   ;
    strcpy(keys[i++],"CMTFILE") ;
    strcpy(keys[i++],"IDEC_2")  ;
    strcpy(keys[i++],"DMIN")    ;
    strcpy(keys[i++],"DMAX")    ;

    get_i_master(i_tmp, keys, 4, eq) ;  
    get_cmtf( eq, 0) ;
  
    /* Memory Freeing */
    for(i=0 ; i<4 ; i++)
        free((void*)keys[i]) ;
    free((void**) keys )   ;
    free((void*) i_tmp )   ;
}

void ymd2jul(int yyyy,int mm,int dd, int *jul)
{
    int k ;
    int ndays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    if ( (yyyy%4)   == 0 ) ndays[1] ++ ;
    if ( (yyyy%100) == 0 ) ndays[1] -- ;
    if ( (yyyy%400) == 0 ) ndays[1] ++ ;

    (*jul) = dd ;
    for (k=0; k<mm-1; k++)
        (*jul) += ndays[k] ;
}


void date2epoch(int year, int mm, int dd, int hour, int min, int sec, int msec, double *epoch)
{
    struct tm date;
    time_t     tmp;
    extern long timezone;
    date.tm_sec   = sec  ;
    date.tm_min   = min  ;
    date.tm_hour  = hour ;
    date.tm_mday  = dd   ;
    date.tm_mon   = mm-1 ;
    date.tm_year  = year - 1900 ;
    date.tm_isdst = 0    ;
    tzset();
    tmp           = mktime(&date) ;
    *epoch        = (double)tmp + (double)(msec)/1000. - timezone ;
}

void delta_t(int y1, int j1, int h1, int m1, int s1, int ms1, int y0, int j0, int h0, int m0, int s0, int ms0, double *tdiff)
{
    double t1, t0 ;
    date2epoch(y1,1,1,h1,m1,s1,ms1,&t1) ;
    date2epoch(y0,1,1,h0,m0,s0,ms0,&t0) ;
    *tdiff = t1 - t0 + (j1 - j0)*86400   ;
}



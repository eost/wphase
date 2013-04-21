/***************************************************************************
*
*	              W phase source inversion package 	            
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Université de Strasbourg / CNRS 
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

/* Header of rwsacs.c subroutines*/

#ifndef ITIME
#define ITIME 1
#endif /* not ITIME */

#ifndef IRLIM
#define IRLIM 2
#endif /* not IRLIM */

#ifndef IAMPH
#define IAMPH 3
#endif /* not IAMPH */

#ifndef IXY
#define IXY   4
#endif /* not IXY */

#ifndef IXYZ
#define IXYZ 51
#endif /* not IXYZ */

#ifndef NVHDR
#define NVHDR 6
#endif /* not NVHDR */

typedef struct
   {
     float delta, depmin, depmax, scale, odelta;
     float b, e, o, a, depmen, cmpaz, cmpinc;
     float stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag;
     float dist, az, baz, gcarc;
     float t[10], resp[10], user[10], f;
     float xminimum, xmaximum, yminimum, ymaximum;
     float internal1,internal2, internal3;
     int   nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec, npts;
     int   nvhdr, norid, nevid, nwfid, nxsize, nysize;
     int   iftype, idep, iztype, iinst, istreg, ievreg,ievtyp;
     int   iqual, isynth, imagtyp, imagsrc;
     int   internal4;
     int leven, lpspol, lovrok, lcalda;
     char  kevnm[17], kstnm[9], knetwk[9], khole[9];
     char  ko[9], ka[9], kt[10][9], kf[9], kuser[3][9], kcmpnm[9];
     char  kdatrd[9], kinst[9];
   } sachdr ;


/*************************************/
/* Allocates memory for a sac header */
void hdr_alloc(sachdr *hdr);

/*********************************************/
/* Allocates memory for a tab of sac headers */
void hdr_tab(sachdr **hdr, int n);

/**************************************************/
/*          rhdrsac(filename, h, ierror)          */
/**************************************************/
/* Read sac file header                           */
/* Input  : filename                              */
/*          ierror                                */
/*	     >if *ierror=1 (input) and if the     */
/*            file does not exist then an         */
/*	      error occur and the program quit    */
/*	     >if *ierror=0 (input) and if the     */
/*	      file does not exist, then ierror    */
/*	      is returned as *ierror = 1 (output) */
/*                                                */
/* Output : h : a struct sachdr                   */
void rhdrsac(char *file, sachdr *hdr, int *ierror);


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
/*	      error occur and the program quit    */
/*	     >if *ierror=0 (input) and if the     */
/*	      file does not exist, then ierror    */
/*	      is returned as *ierror = 1 (output) */
/*                                                */
/* Output : data : pointer to the proto_allocd    */
/*                 data array                     */
void rdatsac(char *file, sachdr *hdr, double *data, int *ierror);

/************************************************/
/*            whdr(f, h)                        */
/************************************************/
/* Write sac file header                        */
/* Input  : f : File descriptor                 */
/*          h : the struct sachdr               */
int whdr(FILE *fd, sachdr *hdr);

/***************************************************/
/*         minimax(y, np, ymin, ymax)              */
/***************************************************/
/* Calculate ymin,ymax from the np first samples   */
/* in y                                            */
void minimax(float *y, int np, float *ymin, float *ymax) ;

/***************************************************/
/*              wdat(f, h, data)                   */
/***************************************************/
/* Write data in sac file                          */
/* Input  : f : file descriptor                    */
/*          h : the struct sachdr                  */
/*          data : array of data points            */
int wdat(FILE *fd, sachdr *hdr, double *data);

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
void whdrsac(char *file, sachdr *hdr);


/***************************************************/
/*             wsac(filename, h, data)             */
/***************************************************/
/* Write data in sac file                          */
/*   > The sac file must already exist otherwise   */
/*     an error signal is sended.                  */
/* Input  : filename                               */
/*          h : the struct sachdr                  */
/*          data : array of data points            */
/* Output : data points in the sac file "filename" */
void wsac(char *file, sachdr *hdr, double *data);



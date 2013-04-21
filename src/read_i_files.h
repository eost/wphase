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

/* Header of read_i_files.c subroutines */

#ifndef POW
#define POW 1.e28
#endif /* not POW */


typedef struct
{
  char   pdeline[__LSIZE__], cmtfile[__FSIZE__], seed[__FSIZE__] ; 
  char   gf_dir[__FSIZE__], evnm[__FSIZE__], evid[__LSIZE__] ;
  int    ot_ye, ot_mo, ot_dm, ot_ho, ot_mi, ot_se, ot_ms ;
  int    filtorder, nf, filtnpass, idtr ;
  double pde_evla, pde_evlo, pde_evdp ;
  double evla, evlo, evdp, dmin, dmax ;
  double fl, fh, tol,  flow, fhigh    ;
  double preevent, fend, ts, hd ;
  double *wp_win4, **vm   ;
} str_quake_params ;



void add_slash(char *c);

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
int get_cmtf(str_quake_params *eq, int flag) ;


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
void get_i_master(char *file, char **keys, int n, str_quake_params *eq) ;
void decode_wp_win(char *buffer, double *wp_win4);
void wp_time_window(double gcarc, double *wp_win4, double *twp_beg, double *twp_end) ;

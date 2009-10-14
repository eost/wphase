/* Header of read_i_files.c subroutines */

#ifndef POW
#define POW 1.e28
#endif /* not POW */


typedef struct
{
  char   pdeline[__LSIZE__], cmtfile[__FSIZE__], seed[__FSIZE__] ; 
  char   gf_dir[__FSIZE__], evnm[__FSIZE__], evid[32] ;
  int    ot_ye, ot_mo, ot_dm, ot_ho, ot_mi, ot_se, ot_ms ;
  int    filtorder, nf, filtnpass     ;
  double pde_evla, pde_evlo, pde_evdp ;
  double evla, evlo, evdp, dmin, dmax ;
  double fl, fh, tol,  flow, fhigh    ;
  double preevent, idtr, fend, ts, hd ;
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
void wp_time_window(double *gcarc, double *wp_win4, double *twp_beg, double *twp_end) ;

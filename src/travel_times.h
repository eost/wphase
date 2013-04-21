/****************************************************************
*	W phase package - travel_times header
*                                           
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#define  NDEPTHS 15

#define  NDISTAS 118

void make_table(char* name, char* table) ;

void trav_time_init(int nh, int nd, double h, double *dv, double *tv, int *ierror) ;

void trav_time(double d, double *tv, double *dv, int nd, double *P_tt, int *ierror);

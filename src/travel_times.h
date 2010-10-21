/****************************************************************
*	W phase package - travel_times header
*                                           
*       History
*             2010  Original Coding
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#ifndef SAFETY_DELAY
#define SAFETY_DELAY __SAFETY_DELAY__
#endif /* not SAFETY_DELAY */


#define  NDEPTHS 15

#define  NDISTAS 118


void make_table(char* name, char* table) ;

void trav_time_init(int *nh, int *nd, double *h, double *dv, double *tv, int *ierror) ;

void trav_time(double *d, double *tv, double *dv, int *nd, double *P_tt, int *ierror);

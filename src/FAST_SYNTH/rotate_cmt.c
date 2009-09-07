#include <math.h>
#define  DEG2RAD   M_PI/180.

/* Rotate moment tensor M_cmt to N_cmt according to az */
void 
rotate_cmt(double *M_cmt, double *N_cmt, double az)
{
  double co, si, co2, si2;
  
  co  = cos(   az*DEG2RAD);
  si  = sin(   az*DEG2RAD);
  co2 = cos(2.*az*DEG2RAD);
  si2 = sin(2.*az*DEG2RAD);

  N_cmt[0] = M_cmt[0];                                       //RR
  N_cmt[1] = M_cmt[1]*co*co + M_cmt[2]*si*si - M_cmt[5]*si2; //TT
  N_cmt[2] = M_cmt[2]*co*co + M_cmt[1]*si*si + M_cmt[5]*si2; //PP
 
  N_cmt[3] =  M_cmt[3]*co - M_cmt[4]*si;                     //RT
  N_cmt[4] =  M_cmt[3]*si + M_cmt[4]*co;                     //RP
  N_cmt[5] = (M_cmt[1]-M_cmt[2])*si2/2. + M_cmt[5]*co2;      //TP
}

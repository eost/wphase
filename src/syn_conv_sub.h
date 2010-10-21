/****************************************************************
*	W phase package - syn_conv_sub headers
*                                           
*       History
*             2010  Original Coding
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

void make_stf(int lh, char *type, double *h) ;

void runave(double *x, double *y, int npts, int lh, char *type) ;

void conv_by_stf(double *delay, double *half_width, char *itype, sachdr *hdr, double *x_in, double *x_conv) ;



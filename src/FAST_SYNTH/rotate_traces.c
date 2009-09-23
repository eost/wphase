#include <math.h>
#include <stdlib.h>

void 
rotate_traces(double *T, double *P, float baz, int npts, double *N, double *E)
{
  double co, si;
  int j;

  co = cos(M_PI*(double)baz/180.);
  si = sin(M_PI*(double)baz/180.);

  for(j=0; j<npts; j++)
	{
	  N[j] = -co*T[j] - si*P[j];
	  E[j] = -si*T[j] + co*P[j];
	}
}

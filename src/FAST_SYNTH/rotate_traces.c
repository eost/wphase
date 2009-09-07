#include <math.h>
#include <stdlib.h>

void 
rotate_traces(float *T, float *P, float baz, int npts, float *N, float *E)
{
  double co, si;
  int j;

  co = cos(M_PI*baz/180.);
  si = sin(M_PI*baz/180.);

  for(j=0; j<npts; j++)
	{
	  N[j] = -co*T[j] - si*P[j];
	  E[j] = -si*T[j] + co*P[j];
	}
}

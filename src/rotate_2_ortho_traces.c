/***************************************************************************
*
*	              W phase source inversion package 	            
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

#include <math.h>
#include <stdlib.h>

void 
rotate_2_ortho_traces(double *T, double *P, float baz, int npts, double *N, double *E)
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

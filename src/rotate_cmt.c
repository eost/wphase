/***************************************************************************
*
*                     W phase source inversion package              
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
#define  DEG2RAD   M_PI/180.

/* Rotate moment tensor M_cmt to N_cmt according to az */
void rotate_cmt(double *M_cmt, double *N_cmt, double az)
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

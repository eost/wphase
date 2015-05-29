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

/*     Convert date to julian days     */

int yyyymmdd2jjj(int year, int month, int day)
{
  int ndays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int k, jul;
      
  if ( year%4   == 0 ) 
      ndays[1] += 1;
  if ( year%100 == 0 ) 
      ndays[1] -= 1;
  if ( year%400 == 0 ) 
      ndays[1] += 1;
      
  jul = day;
  for(k=0; k<month-1; k++)
      jul += ndays[k];
  
  return(jul);
}

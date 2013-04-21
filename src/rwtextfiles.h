/***************************************************************************
*
*	              W phase source inversion package 	            
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Université de Strasbourg / CNRS 
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

/*********************************************************/
/* Return nb of blank in line before the first character */
int  nb_blank(char *line);


/************************************************************/
/* Return nb of character in "in" without the ending blanks */
int nbchar(char *in);


/*******************************/
/* Count nb of lines in a file */
int count_lines(FILE *f);


/********************************************/
/* Check if flag = nv                       */
void check_scan(int nv, int flag, char *file, FILE *P_ttfile);


/********************************************/
/* Open a file "filename" to read text.     */
/* Return the file stream "fstream" and the */
/* number of lines "nl"                     */
FILE *openfile_rt(char *filename, int *nl);

/*****************************************/
/* Open "filename" to write a text file  */
FILE *openfile_wt(char *filename);

/********************************************************/
/* Set standard GF filename                             */
char *get_gf_filename(char *dir, char *stnm, char *netwk, char chan, 
		      char *ext) ;



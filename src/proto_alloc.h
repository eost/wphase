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

#ifndef FSIZE
#define FSIZE  (int)__FSIZE__
#endif /* not FSIZE */

#ifndef LSIZE
#define LSIZE  (int)__LSIZE__
#endif /* not FSIZE */

#ifndef IDSIZE
#define IDSIZE  (int)__IDSIZE__
#endif /* not IDSIZE */


/********************************************/
/*Allocates memory for a tab of n chars     */
char *char_alloc(int n);

char *char_calloc(int n);


/**************************************************/
/*Allocates memory for a tab of n pointer to char */
char **char_alloc2p(int n);

/********************************************/
/*Allocates memory for a tab of n x m chars */
char **char_alloc2(int n, int m); 

char **char_calloc2(int n, int m); 

/**********************************************************/
/* Allocates memory for a tab of n float  precision reals */
float *float_alloc(int n);

float *float_calloc(int n);

/*********************************************************/
/*Allocates memory for a tab of n double precision reals */
double *double_alloc(int n);

double *double_calloc(int n);

/****************************************/
/* Allocate a tab of pointer to pointer */
double **double_alloc2p(int n);

double ***double_alloc3p(int n);

/*********************************************************/
/*Allocates memory for a tab of n double precision reals */
double **double_alloc2(int n, int m);

double **double_calloc2(int n, int m) ;

/******************************************/
/*Allocates memory for a tab of n x m int */
int **int_alloc2(int n, int m);

/**************************************/
/*Allocates memory for a tab of n int */
int *int_alloc(int n);

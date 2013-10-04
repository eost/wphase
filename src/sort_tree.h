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

/**********************************************/
/*  "Tree" sorting subroutines of sac files   */
/* in order of increasing epicentral distance */
/**********************************************/
/* Header : sort_tree.h                       */
/* Source Code : sort_tree.c                  */
/**********************************************/

struct tree
{
  char   file[FSIZE], sta[9], net[9];
  char   cmp[9], locid[9]; 
  double cmpaz, az, xdeg,stla,stlo,stel;
  int    occur, npts;
  struct tree   *g, *d;
};




/**************************************************/
/*             new = alloctree()                  */
/**************************************************/
/* Default size memory allocation for struct tree */
struct tree *alloctree();


/*********************************************************/
/*               splithdr(ch, hdr, new)                  */
/*********************************************************/
/* Fill the struct tree new (already proto_allocd) from sac */
/* header hdr and a filename "ch"                        */
void splithdr(char *ch, sachdr *hdr, struct tree *new);


/***************************************************************/
/*                   newnode = addfile(oldnode)                */
/***************************************************************/
/* Copy oldnode into newnode (a struct tree newly proto_allocd    */
/* and returned in output of addfile.                          */
struct tree *addfile(struct tree *mod);


/***************************************************************/
/*                   build (root, model)                       */
/***************************************************************/
/* Recursive Sorting of sac files in order of increasing       */
/* epicentral distances.                                       */
/* sac file names are stored properly as each node of a tree   */
/* Input -> root  : pointer to the curent node in the tree     */
/*       -> model : pointer to node corresponding to the new   */
/*                  sac file to be added                       */
void build (struct tree * root, struct tree *mod, int flag);


/*****************************************************************/
/*          savetree(tree, stream, DMIN, DMAX)                   */
/*****************************************************************/
/* Write tree on the stream in a way the sac files names apprear */
/* in order of increasing distances between DMIN and DMAX (deg)  */
void savetree(struct tree *root, FILE *o_sacf, double *DMIN, double *DMAX);


/***********************************************************/
/*                     disptree(tree)                      */
/***********************************************************/
/* Display the whole tree in order of increasing distances */
void disptree (struct tree *root);


/***********************************************************/
/*                     freetree(tree)                      */
/***********************************************************/
/* freeing memory for the whole tree                       */
void freetree (struct tree *root);

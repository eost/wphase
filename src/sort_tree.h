/****************************************************************
*	W phase package - "Tree" sorting header
*                                           
*       History
*             2010  Original Coding
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

struct tree
{
  char   file[FSIZE], sta[9], net[9];
  char   cmp[9], locid[9]; 
  double az, xdeg,stla,stlo,stel;
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

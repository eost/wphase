/****************************************************************
*	W phase package - "Tree" sorting subroutines
*                                           
*       History
*             2010  Original Coding
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rwtextfiles.h" /* rwtextfiles.c */
#include "proto_alloc.h"    /* proto_alloc.c    */
#include "rwsacs.h"      /* rwsacs.c      */
#include "sort_tree.h"  


/**************************************************/
/*             new = alloctree()                  */
/**************************************************/
/* Default size memory allocation for struct tree */
struct tree *alloctree()
{
  struct tree *new;
  if (( new=(struct tree*)malloc(sizeof(struct tree)) )==NULL)
    {
      fprintf (stderr,"Error alloc struct tree\n");
      exit (1);
    }
  return new;
}


/*********************************************************/
/*               splithdr(ch, hdr, new)                  */
/*********************************************************/
/* Fill the struct tree new (already proto_allocd) from sac */
/* header hdr and a filename "ch"                        */
void 
splithdr(char *ch, sachdr *hdr, struct tree *new)
{
  int nb ;

  /* Set filename */
  strcpy (new->file,ch) ;

  /* Set staname */
  nb = nbchar(hdr->kstnm) ;
  strncpy(new->sta,hdr->kstnm,nb) ;
  new->sta[nb]='\0' ;

  /* Set network */
  nb = nbchar(hdr->knetwk) ;
  strncpy(new->net,hdr->knetwk,nb) ;
  new->net[nb]='\0' ;

  /* Set component */
  nb = nbchar(hdr->kcmpnm) ;
  strncpy(new->cmp,hdr->kcmpnm,nb) ;
  new->cmp[nb]='\0' ;

  /* Set locid */
  nb = nbchar(hdr->khole) ;
  if ( nb == 0 )
    strcpy(new->locid, "--") ;
  else
    {
      strncpy(new->locid,hdr->khole,nb) ;
      new->locid[nb] = '\0';
    }

  /* Set occurence and npts */
  new->occur = 1 ;
  new->npts = hdr->npts ;

  /* Set azimuth and great circle arc lengths */
  new->az   = (double)hdr->az ;
  new->xdeg = (double)hdr->gcarc ;
  new->stla = (double)hdr->stla ;
  new->stlo = (double)hdr->stlo ;
  new->stel = (double)hdr->stel ;  
  /* Pointers to branches*/
  new->d    = NULL ;
  new->g    = NULL ;
}

/***************************************************************/
/*                   newnode = addfile(oldnode)                */
/***************************************************************/
/* Copy oldnode into newnode (a struct tree newly proto_allocd    */
/* and returned in output of addfile.                          */
struct tree *
addfile(struct tree *mod)
{
  struct tree *new;
  
  if (( new=(struct tree*)malloc(sizeof(struct tree)) )==NULL)
    {
      fprintf (stderr,"Error sorting files : alloc file %s\n",new->file);
      exit (1);
    }

  /* Set tree */
  strcpy (new->file,mod->file) ;
  strcpy(new->sta,mod->sta) ;
  strcpy(new->net,mod->net) ;
  strcpy(new->cmp,mod->cmp) ;
  new->occur = mod->occur ;
  strcpy(new->locid,mod->locid) ;
  new->az   = mod->az ;
  new->xdeg = mod->xdeg ;
  new->stla = mod->stla ;
  new->stlo = mod->stlo ;
  new->stel = mod->stel ;
  new->npts = mod->npts ;  
  new->d    = NULL ;
  new->g    = NULL ;

  return new;
}


/***************************************************************/
/*                   build (root, model)                       */
/***************************************************************/
/* Recursive Sorting of sac files in order of increasing       */
/* epicentral distances.                                       */
/* sac file names are stored properly as each node of a tree   */
/* Input -> root  : pointer to the curent node in the tree     */
/*       -> model : pointer to node corresponding to the new   */
/*                  sac file to be added                       */
void 
build (struct tree * root, struct tree *mod, int flag)
{
  int    i,j,k;
  struct tree   *beg = root;
  struct tree   *new = mod;

  
  if (flag <= 0)
    i = strcmp(new->net,beg->net) ;
  else
    i = 0 ;
  j   = strcmp(new->sta,beg->sta)   ;
  k   = strcmp(new->cmp,beg->cmp)   ;

  if ((i == 0) && (j == 0) && (k == 0)) /* same network, station and channel */
    {
      if (strcmp(new->locid,beg->locid) == 0)/* same locid => priority to largest npts */
	{
	  beg->occur++ ;
	  if (beg->npts < new->npts) 
	    { 
	      beg->npts = new->npts       ; 
	      strcpy(beg->file,new->file) ;
	      if (flag > 0)
		strcpy(beg->net,new->net) ;
	    }
	}
      else if ( flag >= 0 )
	{
	  if (strcmp(beg->locid,"--") == 0) /* Priority to empty locid */
	    return ;
	  else if (strcmp(beg->locid,"00") == 0) /* then 00 */
	    {
	      if (strcmp(new->locid,"--") == 0)
		{
		  if ((new->xdeg) != (beg->xdeg))
		    fprintf(stderr,"Warning (sort tree): station position in files : %s - %s\n",new->file,beg->file);
		  strcpy(beg->file,new->file)   ;
		  strcpy(beg->locid,new->locid) ;
		  beg->npts = new->npts         ;
		  beg->occur = 1                ;
		  if (flag > 0)
		    strcpy(beg->net,new->net) ;
		}
	    }
	  else if (strcmp(beg->locid,"10") == 0) /* then 10 */
	    { 
	      if ( (strcmp(new->locid,"--") == 0) || (strcmp(new->locid,"00") == 0) )
		{
		  if ((new->xdeg) != (beg->xdeg))
		    fprintf(stderr,"Warning (sort tree): station position in files : %s - %s\n",new->file,beg->file);
		  beg->occur = 1;
		  strcpy(beg->file,new->file)   ;
		  strcpy(beg->locid,new->locid) ;
		  if (flag > 0)
		    strcpy(beg->net,new->net)   ;
		}
	    }
	  else                                  /* then largest npts */
	    {
	      i = strcmp(new->locid,"--") ; 
	      j = strcmp(new->locid,"00") ;
	      k = strcmp(new->locid,"10") ;
	      if ( (i==0) || (j==0) || (k==0) )
		{
		  if ((new->xdeg) != (beg->xdeg))
		    fprintf(stderr,"Warning (sort tree): station position in files : %s - %s\n",new->file,beg->file);
		  beg->occur = 1               ;
		  strcpy(beg->file,new->file)  ;
		  strcpy(beg->locid,new->locid);
		  if ( flag > 0)
		    strcpy(beg->net,new->net)  ;
		}
	      else if (beg->npts < new->npts) 
		{ 
		  beg->occur = 1                ;
		  beg->npts  = new->npts        ; 	      
		  strcpy(beg->file,new->file)   ;
		  strcpy(beg->locid,new->locid) ;
		  if ( flag > 0)
		    strcpy(beg->net,new->net)   ;
		}
	    }
	}
      else
	{
	  if (beg->g == NULL)
	    beg->g = addfile(new) ;
	  else
	    build(beg->g, new, flag);
	}
    }
  else if ((i == 0) && (j == 0) && (k != 0)) /* same station but different channel */
    {
      if (beg->d == NULL)
	  beg->d = addfile(new) ;
      else
	  build(beg->d, new, flag) ;
    }
  else if ((new->xdeg) <= (beg->xdeg))
    {
      if (beg->g == NULL)
	  beg->g = addfile(new) ;
      else
	  build(beg->g, new, flag);
    }
  else if ((new->xdeg) > (beg->xdeg)) 
    {
      if (beg->d == NULL)
	  beg->d = addfile(new) ;
      else
	  build(beg->d, new, flag) ;
    }
}



/*****************************************************************/
/*          savetree(tree, stream, DMIN, DMAX)                   */
/*****************************************************************/
/* Write tree on the stream in a way the sac files names apprear */
/* in order of increasing distances between DMIN and DMAX (deg)  */
void 
savetree(struct tree *root, FILE *o_sacf, double *DMIN, double *DMAX)
{
  struct tree *beg=root; 
  
  if (beg)
    {
      savetree (beg->g, o_sacf, DMIN, DMAX);
      if ( (beg->xdeg >= *DMIN) && (beg->xdeg <= *DMAX) )
	{	         //file sta  net  cmp  stla   stlo stel     az     xdeg
	  fprintf(o_sacf,"%-65s %-9s %-9s %-9s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
		  beg->file,beg->sta,beg->net,beg->cmp,beg->stla,beg->stlo,beg->stel,
		  beg->az, beg->xdeg);
	  if (beg->occur > 1)
	    {
	      fprintf(stderr,"Warning (sort tree): multiple data files for station %s and channel %s\n", beg->sta,beg->cmp) ;
	      fprintf(stderr,"    ... only %s is considered.\n",beg->file) ;
	    }
	}
      savetree (beg->d, o_sacf, DMIN, DMAX);
    } 
}



/***********************************************************/
/*                     disptree(tree)                      */
/***********************************************************/
/* Display the whole tree in order of increasing distances */
void 
disptree (struct tree *root)
{
  struct tree *beg=root;
  if (beg)
    {
      disptree (beg->g);
      printf ("%45s %8s %8s %8s %8s %12.4f %12.4f %4d\n", beg->file, beg->sta, beg->net, beg->cmp, beg->locid, beg->az, beg->xdeg, beg->occur);
      //printf ("%45s %12.4f\n", beg->file, beg->xdeg);
      disptree (beg->d);
    }
}


/***********************************************************/
/*                     freetree(tree)                      */
/***********************************************************/
/* freeing memory for the whole tree                       */
void 
freetree (struct tree *root)
{
  if (root)
    {
      freetree (root->g);
      freetree (root->d);
      free((void *)root);
    }
}

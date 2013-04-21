/****************************************************************
*	W phase package - Memory allocation subroutines
*                                           
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "proto_alloc.h"



/******************************************/
/*             v=char_alloc(n)            */
/******************************************/
/*Allocates memory for a tab of n chars   */
char *
char_alloc(int n)
{
  char *v;
  if ((v = (char * ) malloc  ( n * sizeof(char))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}


/******************************************/
/*            v=char_calloc(n)            */
/******************************************/
/*Allocates memory for a tab of n chars   */
char *
char_calloc(int n)
{
  char *v;
  if ((v = (char * ) calloc  ( n , sizeof(char))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}


/**************************************************/
/*               v=char_alloc2p(n)                */
/**************************************************/
/*Allocates memory for a tab of n pointer to char */
char **
char_alloc2p(int n)
{
  char **v;
  if ((v = (char ** ) malloc  ( n * sizeof(char*))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}


/********************************************/
/*            v=char_alloc2(n)              */
/********************************************/
/*Allocates memory for a tab of n x m chars */
char **
char_alloc2(int n, int m)
{
  char **v;
  int  i;
  v = char_alloc2p(n);
  for (i=0 ; i<n ; i++)
    v[i] = char_alloc(m);
 return v;
}


/********************************************/
/*           v=char_calloc2(n)              */
/********************************************/
/*Allocates memory for a tab of n x m chars */
char **
char_calloc2(int n, int m)
{
  char **v;
  int  i;
  if ((v = (char ** ) calloc  ( n , sizeof(char*))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
  for (i=0 ; i<n ; i++)
    if ((v[i] = (char *) calloc (m , sizeof(char))) == NULL)
      {
	fprintf(stderr,"FATAL ERROR: Out of memory\n");
	exit(1);
      }
 return v;
}



/**********************************************************/
/*                    v=float_alloc(n)                    */
/**********************************************************/
/* Allocates memory for a tab of n float  precision reals */
float *
float_alloc(int n)
{
  float *v;
  if ((v = (float * ) malloc  ( n * sizeof(float))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
  return v;
}

/**********************************************************/
/*                    v=float_calloc(n)                   */
/**********************************************************/
/* Allocates memory for a tab of n float  precision reals */
float *
float_calloc(int n)
{
  float *v;
  if ((v = (float * ) calloc  ( n , sizeof(float))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
  return v;
}

/*********************************************************/
/*                   v=double_alloc(n)                   */
/*********************************************************/
/*Allocates memory for a tab of n double precision reals */
double *
double_alloc(int n)
{
  double *v;
  if ((v = (double * ) malloc  ( n * sizeof(double))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}

/*********************************************************/
/*                  v=double_calloc(n)                   */
/*********************************************************/
/*Allocates memory for a tab of n double precision reals */
double *
double_calloc(int n)
{
  double *v;
  if ((v = (double * ) calloc  ( n , sizeof(double))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}


/***************************************/
/*         v=double_alloc2p(n)         */
/***************************************/
/* Allocate a tab of pointer to pointer*/
double **
double_alloc2p(int n)
{
  double **v;
  if ((v = (double ** ) malloc  ( n * sizeof(double*))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}


/***************************************/
/*         v=double_alloc3p(n)         */
/***************************************/
/* Allocate a tab of pointer to pointer*/
double ***
double_alloc3p(int n)
{
  double ***v;
  if ((v = (double *** ) malloc  ( n * sizeof(double**))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}



/**********************************************/
/*            v=double_alloc2(n,m)            */
/**********************************************/
/*Allocates memory for a tab of n x m doubles */
double **
double_alloc2(int n, int m)
{
  double **v;
  int  i;
  v = double_alloc2p(n) ;
  for (i=0 ; i<n ; i++)
    v[i] = double_alloc(m);
  return v;
}


/**********************************************/
/*           v=double_calloc2(n,m)            */
/**********************************************/
/*Allocates memory for a tab of n x m doubles */
double **
double_calloc2(int n, int m)
{
  double **v;
  int  i;
  if ((v = (double **) calloc  ( n , sizeof(double *))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
  for (i=0 ; i<n ; i++)
    if ((v[i] = (double *) calloc (m , sizeof(double))) == NULL)
      {
	fprintf(stderr,"FATAL ERROR: Out of memory\n");
	exit(1);
      }
 return v;
}

int **
int_alloc2(int n, int m)
{
  int **v, i;
  if ((v = (int **) malloc  ( n * sizeof(int *))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
  for(i=0; i<n; i++)
    if ((v[i] = (int *) malloc  ( m * sizeof(int))) == NULL)
      {
	fprintf(stderr,"FATAL ERROR: Out of memory\n");
	exit(1);
      }

 return v;
}


int *
int_alloc(int n)
{
  int *v;
  if ((v = (int * ) malloc  ( n * sizeof(int))) == NULL)
    {
      fprintf(stderr,"FATAL ERROR: Out of memory\n");
      exit(1);
    }
 return v;
}


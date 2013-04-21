/****************************************************************
*	W phase package - memory allocation headers
*                                           
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

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

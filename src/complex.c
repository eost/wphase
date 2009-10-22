#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"


/***********************************************************/
/*              v=complex_alloc(n)                         */
/***********************************************************/
/*Allocates memory for a tab of n double precision complex */
complex *
complex_alloc(int n)
{
  complex *v;
  if ((v = (complex * ) malloc  ( n * sizeof(complex))) == NULL)
    fprintf(stderr,"FATAL ERROR: Out of memory");
  return v;
}


complex 
cmplx (double u, double v)
{
  complex w;
  
  w.real = u ;
  w.imag = v ;
  
  return (w) ;
}

double 
abs_c (complex u)
{
  double m;
  
  m = sqrt (pow(u.real,2.) + pow(u.imag,2.)) ;
  return (m) ;
}

complex add_c (complex u, complex v)
{
  complex w;

  w.real = u.real + v.real ;
  w.imag = u.imag + v.imag ;

   return (w) ;
}

/*************************************************************************/
/*               sub_c                                                   */
/*************************************************************************/
/*
   Routine to subtract two complex numbers

       w = sub_c(u,v)
       w = u - v
*/
complex 
sub_c(complex u,complex v)
{
  complex w ;
  
  w.real = u.real - v.real ;
  w.imag = u.imag - v.imag ;
  
   return (w) ;
}

/***********************************************************************/
/*             mul_c                                                   */
/***********************************************************************/
/*
       Routine to multiply two complex numbers

               w = mul_c(u,v)
               w = u * v

*/
complex 
mul_c (complex u, complex v)
{
  complex w ;

  w.real = u.real*v.real - u.imag*v.imag ;
  w.imag = u.real*v.imag + u.imag*v.real ;
  
  return (w) ;
}

/***********************************************************************/
/*            cmul_c (a,u)                                             */
/***********************************************************************/
/*
   Routine to multiply a real number times a complex number

       w = cmul_c (a,u)

       a - real number
       u - complex number
*/
complex     
cmul_c (double a, complex u)
{
   complex     w ;

   w.real = a * u.real ;
   w.imag = a * u.imag ;

   return (w) ;
}

/******************************************************************/
/*                    div_c                                       */
/******************************************************************/
/*
   Routine to divide two complex numbers

       w = div_c(u,v)
       w = u/v

*/
complex 
div_c (complex u, complex v, int *ierror)
{
  complex     w ;
  
  /*   check for divide by 0    */
  if (v.real != 0 || v.imag != 0)
    {
      w.real = ((u.real * v.real) + (u.imag * v.imag)) /
	((v.real * v.real) + (v.imag * v.imag)) ;
      w.imag = ((u.imag * v.real) - (u.real * v.imag)) /
	((v.real * v.real) + (v.imag * v.imag)) ;
      
      return (w) ;
    }
  else
    {
      fprintf (stderr, "ERROR: complex division by 0 in div_c\n") ;
      if (*ierror == 1)
	exit (1) ;
      *ierror = 1 ;
      return (w) ;
    }
}


/***************************************************************/
/*                 conj_c                                      */
/***************************************************************/
/*
       Routine to calculate the complex conjugate

               w = conjugate(u)

*/
complex 
conj_c (complex u)
{
       complex         w ;

       w.real = u.real ;
       w.imag = -u.imag ;

       return (w) ;
}

complex 
sqrt_c (complex u)
{
       complex w ;
       double  arg, norm ; 
       arg    = atan2(u.imag, u.real)/2. ;
       norm   = sqrt(sqrt(u.real*u.real+u.imag*u.imag)) ;
       w.real = norm*cos(arg) ;
       w.imag = norm*sin(arg) ;
       return (w) ;
}



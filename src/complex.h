typedef struct {
  float re;
  float im;
} complexf;


typedef struct
   {
       double real ;
       double imag ;
   } complex ;


/***********************************************************/
/*              v=complex_alloc(n)                         */
/***********************************************************/
/*Allocates memory for a tab of n double precision complex */
complex *complex_alloc(int n);

complex cmplx (double u, double v);

double abs_c (complex u);

complex add_c (complex u, complex v);


/*************************************************************************/
/*               sub_c                                                   */
/*************************************************************************/
/*
   Routine to subtract two complex numbers

       w = sub_c(u,v)
       w = u - v
*/
complex sub_c(complex u,complex v);


/***********************************************************************/
/*             mul_c                                                   */
/***********************************************************************/
/*
       Routine to multiply two complex numbers

               w = mul_c(u,v)
               w = u * v

*/
complex mul_c (complex u, complex v);

/***********************************************************************/
/*            cmul_c (a,u)                                             */
/***********************************************************************/
/*
   Routine to multiply a real number times a complex number

       w = cmul_c (a,u)

       a - real number
       u - complex number
*/
complex     cmul_c (double a, complex u);


/******************************************************************/
/*                    div_c                                       */
/******************************************************************/
/*
   Routine to divide two complex numbers

       w = div_c(u,v)
       w = u/v

*/
complex div_c (complex u, complex v, int *ierror);


/***************************************************************/
/*                 conj_c                                      */
/***************************************************************/
/*
       Routine to calculate the complex conjugate

               w = conjugate(u)

*/
complex conj_c (complex u);


complex sqrt_c (complex u);



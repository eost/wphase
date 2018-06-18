/***************************************************************************
*
*                     W phase source inversion package              
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_alloc.h"
#define NDEPTHS_LOC 200

/* External subroutine */
void get_depths(char *path, double *depths, int *nd);

/******************************************************************/
/*     get_prefix(dep, xdeg, bsegm1, bsegm2, bdepth, bdist )      */
/******************************************************************/
/* Select the best depth and dist for Green's function extraction */
/* Input : dep  : focal depths (from pde)                         */
/*         xdeg : epicentral distance                             */
/* Ouput : bsegm1 : path to best depths GF                        */
/*         bsegm2 : sac file containing the best dist GF          */
/*         bdepths: best focal depths                             */
/*         bdist  : best focal distance                           */
void get_prefix(double cmt_dep, double xdeg, char *best_segm1, char *best_segm2, 
                                double *best_depth, double *best_dist )
{
    char   gf_path[FSIZE];
    double dd,dd_best,depths[NDEPTHS_LOC];
    int    jd, nd=88, idist, idistmax;

    /* ############################################################### */
    /* This section should be updated if the distance grid is modified */
    idistmax = 899;
    #ifdef __GFS_01D__
    idist = (int)floor(10.0*xdeg+0.5);
    #else
    #ifdef __GFS_0005D__
    idist = (int)floor((1000.0*xdeg+2.5)/5.)*5;
    idistmax = 30000;
    #else
    idist = 2*(int)floor(5.*xdeg)+1;
    #endif
    #endif /* not FSIZE */  
  
    if (idist < 1) 
    {
        fprintf(stderr, "Warning distance too short!; using 0.1 deg\n");
        idist = 1;
    }
    if (idist > idistmax) 
    {
        fprintf(stderr, "Warning distance too big!;  using 89.9 deg\n");
        idist = idistmax;
    } 
    /* ############################################################### */
  
    if(getenv("GF_PATH"))
        strncpy(gf_path, getenv("GF_PATH"), FSIZE);
    else
    {
        fprintf(stderr, "Error: GF_PATH environment variable not defined\n");
        exit(1);
    }

    get_depths(gf_path, depths, &nd);
  
    /* Finding the best depth */
    *best_depth = depths[0];
    dd_best = fabs(*best_depth-cmt_dep);
    for (jd=1; jd < nd; jd++)
    {
        dd = fabs(depths[jd]-cmt_dep);
        if(dd<dd_best || (dd == dd_best && depths[jd]>*best_depth)) 
        {
            *best_depth = depths[jd];
            dd_best     = dd;
        }
    }
    strncpy(best_segm1, gf_path, FSIZE);
    sprintf(best_segm2, "/H%05.1f",  *best_depth);
    strcat (best_segm1, best_segm2);
    sprintf(best_segm2, "GF.%04d.SY.LH", idist); 
    *best_dist = (float)idist/10.;
    //fprintf(stdout, "Using gf's %s/*/%sZ.SAC\n", best_segm1,best_segm2);
    return;
}

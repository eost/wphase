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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include "decimate_2sps.h"
#include "rwsacs.h"
#include "proto_alloc.h"
#define  EPS 1.e-2


int init_FIR(double *coeffs, int Ncoeffs, FIR_filter *FIR);
int decimate(FIR_filter *FIR, int dec_fac, double *yi, int ni, double *yo, int *no);

void error_decim_sac(char *arg)
{
    fprintf(stderr, "Syntax: %s  i_sac_file_name o_dir [o_sac_file_name]\n", arg);
    exit(1);
}

char* try_update_name_LH(char *name)
{
    // Heuristically try to update the file name to say 'LHZ'
    int k, n;
    n = strlen(name);
    for (k=0; k<n-5; k++)
        if (name[k] == '.' && name[k+4] == '.')
            if (name[k+1] == 'H' || name[k+1] == 'B')
                if (name[k+2] == 'H' || name[k+2] == 'L')
                    if (name[k+3] == 'Z' || name[k+3] == 'N' || name[k+3] == 'E' || name[k+3] == '1' || name[k+3] == '2')
                        name[k+1] = 'L';

    return name;
}

int main(int argc, char *argv[])
{
    int        ierror=1, n, p, p0;
    char       i_sac_filename[256], o_sac_filename[256], o_dir[256];
    FIR_filter FIR2,FIR3,FIR4,FIR5,*FIR;
    double    *y ;
    sachdr     hdr ;
    Cascade   *cascades;

    cascades = (Cascade *)calloc(NCASC, sizeof(Cascade));
    init_casc(&cascades[ 0], N_casc_len, facs_120);
    init_casc(&cascades[ 1], N_casc_len, facs_100);
    init_casc(&cascades[ 2], N_casc_len, facs_090);
    init_casc(&cascades[ 3], N_casc_len, facs_080);
    init_casc(&cascades[ 4], N_casc_len, facs_060);
    init_casc(&cascades[ 5], N_casc_len, facs_050);
    init_casc(&cascades[ 6], N_casc_len, facs_040);
    init_casc(&cascades[ 7], N_casc_len, facs_030);
    init_casc(&cascades[ 8], N_casc_len, facs_020);
    init_casc(&cascades[ 9], N_casc_len, facs_010);
    init_casc(&cascades[10], N_casc_len, facs_002);

    // Parses arguments
    if ( argc < 3 || argc > 4 )
        error_decim_sac(argv[0]);

    strcpy(i_sac_filename, argv[1]);
    strcpy(o_dir, argv[2]);
    n = strlen(o_dir);
    if ( o_dir[n-1] != '/' ) strcat(o_dir, "/");
        strcpy(o_sac_filename, o_dir);
    if ( argc == 4 ) 
        strcat(o_sac_filename, argv[3]);
    else
    {
        strcat(o_sac_filename, basename(i_sac_filename));
        //try_update_name_LH(o_sac_filename);
        n = strlen(o_sac_filename);
        if (strncmp(o_sac_filename+n-5, ".SAC", 4)) 
            strcpy(o_sac_filename+n-3, "1sps.SAC");
    }
    printf("%s\t\t%s\n",   i_sac_filename, o_sac_filename);

    // Inits filters
    init_FIR(dec2,Ndec2,&FIR2);
    init_FIR(dec3,Ndec3,&FIR3);
    init_FIR(dec4,Ndec4,&FIR4);
    init_FIR(dec5,Ndec5,&FIR5);

    // Decimates to 1sps
    hdr_init(&hdr);

    rhdrsac(i_sac_filename, &hdr, &ierror);
    y = double_alloc(hdr.npts);
    rdatsac(i_sac_filename, &hdr, y, &ierror);
    p0 = -1;
    for (p=0; p<NCASC; p++)
    {
        //fprintf( stderr, "%4d %12.6f %12.6f\n", p, 1./hdr.delta, (float)(cascades[p].ratio));
        if ( fabs( cascades[p].ratio * hdr.delta * 2. -  1.)  < EPS )
        {
            p0 = p;
            break;
        }
    }
    if (p0 == -1)
    {
        fprintf(stderr, "Not available decimation filter for file %s (dt = %12.9f); skipping it.\n", i_sac_filename, hdr.delta);
        return 1;
    }
    for (p=0; p<N_casc_len; p++)
    {
        if      (cascades[p0].facs[p] == 1) continue;
        else if (cascades[p0].facs[p] == 2) FIR = &FIR2;
        else if (cascades[p0].facs[p] == 3) FIR = &FIR3;
        else if (cascades[p0].facs[p] == 4) FIR = &FIR4;
        else if (cascades[p0].facs[p] == 5) FIR = &FIR5;
        else
        {
            fprintf(stderr,"ERROR Incorrect decimation factor ");
            fprintf(stderr,"Only 2, 3, 4 or 5 are available.\n");
            continue;
        }
        decimate(FIR, cascades[p0].facs[p], y,hdr.npts,y,&hdr.npts);
        hdr.delta *= (float)cascades[p0].facs[p];
    }
    //hdr.kcmpnm[0] = 'L';
    wsac(o_sac_filename, &hdr,y);
    free((void *)y);  
    return 0;
}

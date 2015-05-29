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

/*      Read/write text files subroutines      */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "proto_alloc.h"
#include "rwtextfiles.h"


int nb_blank(char *line)
{
    int i,N; 
    N = strlen(line);
  
    for (i=0; i<=N; i++)
    {
        if (line[i] != ' ')
            return i ;
    }
    return N;
}


int nbchar(char *in)
{
    int i,N;
    N = strlen(in);

    for (i=1; i<=N; i++)
    {
        if (in[N-i] != ' ')
            return (N-i+1);
    }
    return 0;
}



/******************************************/
/*          count_lines(stream)           */
/******************************************/
/* Count nb of lines in a file            */
int count_lines(FILE *f)
{
    char *line;
    int  nl = 0;
    line = char_alloc(LSIZE);

    while( fgets(line,LSIZE,f) != NULL )
        ++nl;
  
    rewind(f);
    free(line);
    return nl;
}


/********************************************/
/*  check_scan(nv, flag, filename, file)    */
/********************************************/
/* Check if flag = nv                       */
/* if flag != nv then an error message      */
/* is written on stderr and the file stream */
/* is "file" is closed                      */
void check_scan(int nv, int flag, char *file, FILE *fstream)
{
    if (flag != nv)
    { 
        fprintf(stderr, "ERROR: reading %s\n",file); 
        fclose(fstream);
        fflush(stderr);
        exit(1); 
    } 
}



/********************************************/
/*    fstream = openfile_rt(filename,nl)    */
/********************************************/
/* Open a file "filename" to read text.     */
/* Return the file stream "fstream" and the */
/* number of lines nl                       */
FILE *openfile_rt(char *filename, int *nl)
{
    int nc;
    FILE *stream;
    if ((stream = fopen (filename,"rt"))==NULL)
    {
        fprintf (stderr, "ERROR (read) : opening file: %s \n", filename) ;
        exit (1) ;
    }
  
    nc = count_lines(stream);
    if(nc < 1) /* error if less than 1 line */
    {
        fprintf(stderr,"No sacfiles found in file: %s\n",filename) ;
        fclose(stream);
        exit(1);
    }
    *nl = nc ;
    return stream;
}


/*****************************************/
/*     openfile_wt(filename)             */
/*****************************************/
/* Open "filename" to write a text file  */
FILE *openfile_wt(char *filename)
{
    FILE *stream;
 
    if ((stream = fopen (filename,"wt"))==NULL)
    {
        fprintf (stderr, "ERROR (write) : opening file: %s \n", filename) ;
        exit (1) ;
    }
    return stream ;
}
  

/*********************************************************/
/*   C = get_gf_filename(dir, stnm, net, cmpnm, loc, ext)*/
/*********************************************************/
/* Set standard GF filename                              */
/* input params are :  -> dir  : directory               */
/*                     -> stnm : station name            */
/*                     -> net  : network                 */
/*                     -> chan : channel                 */
/*                     -> ext  : extention               */
char *get_gf_filename(char *dir, char *stnm, char *netwk, char *cmpnm, char *loc, char *ext)
{
    int n,m ;
    char *sac_filename ;
    sac_filename = char_alloc(FSIZE);
  
    n = strlen(dir);
    m = nbchar(stnm);
    if (n!=0) 
    {
        if (dir[n-1] != '/')
        {
            strcpy(sac_filename, dir);
            strcat(sac_filename, "/");
        }
        else
        strcpy(sac_filename, dir);
        strncat(sac_filename, stnm,m);  
    }
    else
    {
        strcpy(sac_filename, dir);
        strncat(sac_filename, stnm,m);
    }
    // If loc is empty  (more precisely "  "; it is set to "--")
    if ( strncmp(loc, "  ", 2 ) == 0 )  
        strncpy(loc, "--", 2);
    strcat(sac_filename, "."); n = nbchar(netwk); strncat(sac_filename, netwk,n);
    strcat(sac_filename, "."); n = nbchar(cmpnm); strncat(sac_filename, cmpnm,n);
    strcat(sac_filename, "."); n = nbchar(loc);   strncat(sac_filename, loc,2);
    strcat(sac_filename, "."); n = nbchar(ext);
    if (ext[0] == '.')
        m=1;
    else 
        m=0;
    strncat(sac_filename, &ext[m],n);
    strcat(sac_filename, "");
    return sac_filename ;
}




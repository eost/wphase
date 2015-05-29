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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include <dirent.h>

/********************************************************/
/*             get_depths(path, depths, nd)             */
/********************************************************/
/* Set the available source depths of Green's function  */
/* Input : path to Green's function directory           */
/*                                                      */
/* Output : depths                                      */
/*          nd : nb of depths                           */
/*                                                      */
void get_depths(char *path, double *depths, int *nd) 
{
    struct dirent *entry;
    struct stat file_status;
    char  c_depth[6], name[256], dirname[256];
    int   n,fd;
    DIR   *dir;
  
    strcpy(dirname,path);  
    if (dirname[strlen(dirname)-1] != '/')
        strcat(dirname,"/");

    dir = opendir(path);
    if (!dir) 
    {
        perror(path);
        exit(1);
    }
    n  = 0 ;
    c_depth[5] = '\0';
    fd = dirfd(dir);
    flock(fd, LOCK_EX); // LOCK ON
    while ( (entry = readdir(dir)) != NULL) 
    {
        strcpy(name,dirname);
        strcat(name, entry->d_name);
        if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0))
            if (stat(name, &file_status) == 0)
                if (S_ISDIR(file_status.st_mode) && (entry->d_name[0] == 'H')) 
                {
                    strncpy(c_depth,entry->d_name+1, 5);
                    sscanf(c_depth, "%lf", &depths[n]);
                    n++;
                }
    }
    flock(fd, LOCK_UN); // LOCK OFF
    if (closedir(dir) == -1) 
    {
        perror(path);
        exit(1);
    }
    if(n < 1)
    {
        printf("No H-directories found in %s\n", path);
        exit(1);
    }
    *nd = n;
}


/***************************************************************************
*
*	              W phase source inversion package 	            
*                               -------------
*
*        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*                      
* (c) California Institute of Technology and Université de Strasbourg / CNRS 
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

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<dirent.h>


int rm_dir_rec(const char *path)
{
  DIR *d          = opendir(path);
  size_t path_len = strlen(path);
  int r = -1;
  if (d)
    {
      struct dirent *p;
      r = 0;
      while (!r && (p=readdir(d)))
	{
          int r2 = -1;
          char *buf;
          size_t len;

          /* Skip the names "." and ".." as we don't want to recurse on them. */
          if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
	    continue;
	  
          len = path_len + strlen(p->d_name) + 2; 
          buf = malloc(len);
	  
          if (buf)
	    {
	      struct stat statbuf;
	      snprintf(buf, len, "%s/%s", path, p->d_name);
	      if (!stat(buf, &statbuf))
		{
		  if (S_ISDIR(statbuf.st_mode))
		    r2 = rm_dir_rec(buf);
		  else
		    r2 = unlink(buf);
		}
	      free(buf);
	    }
          r = r2;
	}
      closedir(d);
    }
  if (!r)
    r = rmdir(path);
  return r;
}


void crea_dir(const char *path)
{
  int flag = 0;
  struct stat st;  
  if (stat(path,&st) == 0)
    flag = rm_dir_rec(path);
  if(mkdir(path,0777)==-1 || flag != 0)
    {
      fprintf(stderr,"Error creating directory %s",path);
      exit(1);
    }
}

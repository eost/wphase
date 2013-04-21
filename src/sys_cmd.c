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

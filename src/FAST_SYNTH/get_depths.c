#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
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
void 
get_depths(char *path, float *depths, int *nd) 
{
  struct dirent *entry;
  struct stat file_status;
  char  c_depth[6], name[256], wd[256];
  int   n;
  DIR *dir;
  
  getcwd(wd, 256); // Current working directory

  dir = opendir(path);
		
  if (!dir) {
    perror(path);
    exit(1);
  }
  chdir(path);
  n = 0;
  while ( (entry = readdir(dir)) != NULL) {
    strncpy(name, entry->d_name, 256);
    if ((strcmp(name, ".") != 0) && (strcmp(name, "..") != 0))
      if (stat(name, &file_status) == 0)
	if (S_ISDIR(file_status.st_mode) && (name[0] == 'H')) n++;
  }
  if (closedir(dir) == -1) {
    perror(path);
    exit(1);
  }
  if(n < 1){
    printf("No H-directories found in %s\n", path);
    exit(1);
  }
  
  n = 0;
  dir   = opendir(path);
  c_depth[5] = '\0';
  while ( (entry = readdir(dir)) != NULL) {
    strncpy(name, entry->d_name, 256);
    if ((strcmp(name, ".") != 0) && (strcmp(name, "..") != 0)) {
      if (stat(name, &file_status) == 0) {
	if (S_ISDIR(file_status.st_mode) && (name[0] == 'H')) {
	  strncpy(c_depth, name+1, 5);
	  sscanf(c_depth, "%f", &depths[n]);
	  //printf("%s %f\n", c_depth, depths[n]);
	  n++;
	}
      }
    }
  }
  if (closedir(dir) == -1) {
    perror(path);
    exit(1);
  }
  *nd = n;
  chdir(wd);
}

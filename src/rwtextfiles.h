/****************************************************************
*	W phase package - read/write text file subroutine header
*                                           
*       History
*             2010  Original Coding
*
*       Zacharie Duputel, Luis Rivera and Hiroo Kanamori
*
*****************************************************************/

/*********************************************************/
/* Return nb of blank in line before the first character */
int  nb_blank(char *line);


/************************************************************/
/* Return nb of character in "in" without the ending blanks */
int nbchar(char *in);


/*******************************/
/* Count nb of lines in a file */
int count_lines(FILE *f);


/********************************************/
/* Check if flag = nv                       */
void check_scan(int nv, int flag, char *file, FILE *P_ttfile);


/********************************************/
/* Open a file "filename" to read text.     */
/* Return the file stream "fstream" and the */
/* number of lines "nl"                     */
FILE *openfile_rt(char *filename, int *nl);

/*****************************************/
/* Open "filename" to write a text file  */
FILE *openfile_wt(char *filename);

/********************************************************/
/* Set standard GF filename                             */
char *get_gf_filename(char *dir, char *stnm, char *netwk, char *chan, 
		      char *ext) ;



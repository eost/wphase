/************************************************************
*       List the SAC header fields
*       Usage
*               saclst header_variable_names ... f sac_files ...
*       History
*             1999  Original Coding             Lupei Zhi
*             2003  Further Updates             Qinya Liu
*       28 05 2006  Addition to SAC Codebase    Brian Savage
*                   Combined saclst.c and sacio.c
*       License
*               Distributed under the same License as the SAC Source
*               and binaries, used here under permission by Lupei Zhu
*
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sac.h"

#define OUTPUT_DATE      201
#define OUTPUT_TIME      202
#define OUTPUT_MONTH     204
#define OUTPUT_DAY       205

#define OUTPUT_DEFAULT   301
#define OUTPUT_ALL       302
#define OUTPUT_FULL      303
#define OUTPUT_DEFAULT1  304

#define OUTPUT_PICKS     401

int sac_header_position(const char *s);
void kidate(int year, int jday, int *month, int *day);

void output_header_list() 
{
    fprintf(stderr,"Usage: saclst header_values f file_lists\n");
    fprintf(stderr,"   ex. saclst delta npts kstnm f sacfile1 sacfile2\n");
    fprintf(stderr, "    All Values are case insensitive, except F\n");
    fprintf(stderr, "    If header_values = default  - All Defined Values, 2 Columns\n");
    fprintf(stderr, "    If header_values = default1 - All Defined Values, 1 Column\n");
    fprintf(stderr, "    If header_values = all      - All Values\n");
    fprintf(stderr, "    If header_values = Full     - All Values Formatted (capital F)\n");
    fprintf(stderr, "Available SAC Header Values\n");
    fprintf(stderr, "    \t\tTime-series Values\n");
    fprintf(stderr, "\tb e o a F ko ka kf\n");
    fprintf(stderr, "\tnpts delta depmin depmax depmen scale nvhdr\n");
    fprintf(stderr, "    \t\tStation and Event Values\n");
    fprintf(stderr, "\tkstnm stlo stla stel stdp\n");
    fprintf(stderr, "\tkevnm evlo evla evel evdp\n");
    fprintf(stderr, "\tdist az baz gcarc khole\n");
    fprintf(stderr, "\tkcmpnm knetwk kdatrd kinst cmpaz cmpinc\n");
    fprintf(stderr, "\tiftype idep iztype iinst istreg ievreg ievtyp iqual isynth\n");
    fprintf(stderr, "    \t\tTiming Values\n");
    fprintf(stderr, "\tkzdate kztime odelta\n");
    fprintf(stderr, "\tnzyear nzjday nzmonth nzday nzhour nzmin nzsec nzmsec\n");
    fprintf(stderr, "    \t\tPicks, Response, and User Values\n");
    fprintf(stderr, "\tt0    t1    t2    t3    t4    t5    t6    t7    t8    t9\n");
    fprintf(stderr, "\tkt0   kt1   kt2   kt3   kt4   kt5   kt6   kt7   kt8   kt9\n");
    fprintf(stderr, "\tresp0 resp1 resp2 resp3 resp4 resp5 resp6 resp7 resp8 resp9\n");
    fprintf(stderr, "\tuser0 user1 user2 user3 user4 user5 user6 user7 user8 user9\n");
    fprintf(stderr, "\tkuser0 kuser1 kuser2\n");
  
    exit(-1);
}


int sac_header_value_type(int j) 
{
    if(j >= SAC_HEADER_FLOAT_MIN   && j <= SAC_HEADER_FLOAT_MAX)   return SAC_HEADER_FLOAT_TYPE;
    if(j >= SAC_HEADER_INT_MIN     && j <= SAC_HEADER_INT_MAX)     return SAC_HEADER_INT_TYPE;
    if(j >= SAC_HEADER_ENUM_MIN    && j <= SAC_HEADER_ENUM_MAX)    return SAC_HEADER_ENUM_TYPE;
    if(j >= SAC_HEADER_LOGICAL_MIN && j <= SAC_HEADER_LOGICAL_MAX) return SAC_HEADER_LOGICAL_TYPE;
    if(j == SAC_HEADER_CHAR_DOUBLE)                                return SAC_HEADER_CHAR16_TYPE;
    if(j == SAC_HEADER_CHAR_DOUBLE_END)                            return SAC_HEADER_UNDEFINED_TYPE;
    if(j >= SAC_HEADER_CHAR_MIN    && j <= SAC_HEADER_CHAR_MAX)    return SAC_HEADER_CHAR8_TYPE;
    return SAC_HEADER_UNDEFINED_TYPE;
}  

char *sac_header_value_char(SACHEAD *hd, int j) 
{
    char *c;
    char *cpt;
    float *fpt = (float *) hd;
    int k      = SAC_HEADER_CHAR_MIN + (j - SAC_HEADER_CHAR_MIN) * 2;
    int klen   = (j == SAC_HEADER_CHAR_DOUBLE) ? 16 : 8;

    fpt += k;
    cpt = (char *)fpt;
    c = (char *)malloc(sizeof(char) * (klen + 1));
    c = strncpy(c, cpt, klen);
    c[klen] = '\0';
    return c;
}

float sac_header_value_float(SACHEAD *hd, int j) 
{
    float f;
    float *fpt = (float *) hd;

    fpt += j;
    f = *fpt;
    return f;
}

int sac_header_value_int(SACHEAD *hd, int j) 
{
    int i;
    int *ipt = (int *) hd;

    ipt += j;
    i = *ipt;
    return i;
}

char *sac_header_value_enum(SACHEAD *hd, int j) 
{
    char *s  = NULL;
    int ipt = sac_header_value_int(hd, j);
    if(ipt == SAC_HEADER_INT_UNDEFINED)
        s = strdup(SAC_HEADER_UNDEFINED);
    else if(ipt >= 0 && ipt < SacHeaderEnumsLength)
        s = strdup(SacHeaderEnums[ipt]);
    return s;
}

char * sac_header_value_string(SACHEAD *hd, int j) 
{
    int ipt;
    float fpt;
    char *s = NULL;
  
    switch(sac_header_value_type(j)) 
    {
        case SAC_HEADER_FLOAT_TYPE:
            fpt = sac_header_value_float(hd, j);
            s = (char *)malloc(sizeof(char) * 512);
            s = strncpy(s, "", 0);
            sprintf(s, "%11.5E", fpt);
        break;
        case SAC_HEADER_INT_TYPE:
            ipt = sac_header_value_int(hd, j);
            s = (char *)malloc(sizeof(char) * 512);
            s = strncpy(s, "", 0);
            sprintf(s, "%11d", ipt);
        break;
        case SAC_HEADER_ENUM_TYPE:
            s = sac_header_value_enum(hd, j);
        break;
        case SAC_HEADER_LOGICAL_TYPE:
            s = (char *)malloc(sizeof(char) * 512);
            s = strncpy(s, "", 0);
            sprintf(s, "%s", (sac_header_value_int(hd, j)) ? "TRUE" : "FALSE");
        break;
        case SAC_HEADER_CHAR8_TYPE:
            /* empty */
        case SAC_HEADER_CHAR16_TYPE:
            s = sac_header_value_char(hd, j);
        break;
        default:
            /* empty */
        break;
    }
    return s;
}

int main(int argc, char **argv)
{
    int      i,j,ls[60],nl,k;
    float    fpt;
    int      ipt;
    char    *cpt;
    int      month, day;
    SACHEAD  hd;
    int      def, all, full;
    char    *header_name;
    char    *kstring;
    int      newline;

    if(argc < 2) 
    {
        fprintf(stderr,"Usage: saclst header_lists f file_lists\n");
        fprintf(stderr,"   ex. saclst delta npts kstnm f sacfile1 sacfile2\n");
        fprintf(stderr,"       saclst help - outputs a list of possible values\n"); 
        return -1;
    }
    def = all = full = 0;
    nl=0; argv++; argc--;

    while ( *argv[0] != 'f' ) 
    {
        if(strcasecmp("help", argv[0]) == 0)
            output_header_list();
    
        ls[nl] = sac_header_position(argv[0]);
        if(ls[nl] == OUTPUT_DEFAULT)
            def = 1;
        if(ls[nl] == OUTPUT_ALL)
            all = 1;
        if(ls[nl] == OUTPUT_FULL)
            full = 1;
        if(ls[nl] == OUTPUT_DEFAULT1) 
            def = 2;

        if (ls[nl] == -999) 
        {
            fprintf(stderr, "saclst: Error in header_list  %s\n", argv[0]);
            return -1;
        }
        nl++; argv++; argc--;
    }
    newline = 0;

    if(def || all || full) 
    {
        for(i = 1; i < argc; i++) 
        {
            if( read_sachead(argv[i], &hd) != -1 ) 
            {
                newline = 0;
                printf("File: %s\n", argv[i]);
                for(j = 0; j <= SAC_HEADER_FIELDS; j++) 
                {
                    header_name = SacHeaderName[j];

                    if(full) 
                    {
                        if(j == SAC_HEADER_FLOAT_MIN) 
                            printf(" REAL        INDEX  NAME        Int Value  Real Value\n");
                        if(j == SAC_HEADER_INT_MIN) 
                            printf(" INTEGER     INDEX  NAME        Int Value  Int Value\n");
                        if(j == SAC_HEADER_ENUM_MIN) 
                            printf(" ENUMERATED  INDEX  NAME        Int Value  Enu Value\n");
                        if(j == SAC_HEADER_LOGICAL_MIN) 
                            printf(" LOGICAL     INDEX  NAME        Int Value  Log Value\n");
                        if(j == SAC_HEADER_CHAR_MIN) 
                            printf(" CHARACTER   INDEX  NAME        Int Value  Char Value\n");
                    }

                    switch(sac_header_value_type(j)) 
                    {
                        case SAC_HEADER_FLOAT_TYPE:
                            fpt = sac_header_value_float(&hd, j);
                            if( def && fpt != SAC_HEADER_FLOAT_UNDEFINED ) 
                            {
                                printf("      %-8s  %16.5E", header_name, fpt);
                                newline++;
                            }
                            else if(all)
                                printf("               %3d  %-8s  %16.5E\n", j, header_name, fpt);
                            else if(full)
                                printf("               %3d  %-8s  %16.5E %s\n", j, header_name, fpt, sac_header_value_string(&hd, j));
                        break;
                        case SAC_HEADER_INT_TYPE:
                            /* empty */
                        case SAC_HEADER_ENUM_TYPE:
                            /* empty */
                        case SAC_HEADER_LOGICAL_TYPE:
                            ipt = sac_header_value_int(&hd, j);
                            if( def && ipt != SAC_HEADER_INT_UNDEFINED ) 
                            {
                                printf("      %-8s  %16d", header_name, ipt);
                                newline++;
                            }
                            else if(all)
                                printf("               %3d  %-8s  %16d\n", j, header_name, ipt);
                            else if(full)
                                printf("               %3d  %-8s  %16d %s\n", j, header_name, ipt, sac_header_value_string(&hd, j));
                        break;
                        case SAC_HEADER_CHAR8_TYPE:
                            /* empty */
                        case SAC_HEADER_CHAR16_TYPE:
                            kstring = SAC_HEADER_CHAR_UNDEFINED;
                            cpt = sac_header_value_char(&hd, j);
                            if( def && strncmp(cpt, kstring, 8) != 0 ) 
                            {
                                printf("      %-8s  %16s", header_name, cpt);
                                newline++;
                            }
                            else if(all)
                                printf("               %3d  %-8s  %16s\n", j, header_name, cpt);
                            else if(full)
                                printf("               %3d  %-8s  %16s %s\n", j, header_name, cpt, sac_header_value_string(&hd, j));
                        break;
                    }

                    if(newline == 2 || (def == 2 && newline > 0)) 
                    { 
                        printf("\n");
                        newline = 0;
                    }
                } /* endfor */

                if(def && newline > 0) 
                    printf("\n");
            } /* endif */
        } /* endfor */
        return 0;
    } /* endif */

    for (i=1; i<argc; i++) 
    {    
        if ( read_sachead(argv[i], &hd) != -1) 
        {
            printf("%s ", argv[i]);
            for (j=0; j<nl; j++) 
            {
                switch(sac_header_value_type(ls[j])) 
                {
                    case SAC_HEADER_FLOAT_TYPE:
                        printf("%12.6g",sac_header_value_float(&hd, ls[j]));
                    break;
                    case SAC_HEADER_INT_TYPE:
                        /* empty */
                    case SAC_HEADER_ENUM_TYPE:
                        /* empty */
                    case SAC_HEADER_LOGICAL_TYPE:
                        printf("%10d", sac_header_value_int(&hd, ls[j]));
                    break;
                    case SAC_HEADER_CHAR8_TYPE:
                        /* empty */
                    case SAC_HEADER_CHAR16_TYPE:
                        printf("   %s", sac_header_value_char(&hd, ls[j]));
                    break;
                    default:
                        kidate(hd.nzyear, hd.nzjday, &month, &day);
                        if(ls[j] == OUTPUT_DAY)
                            printf("%10d",day);
                        else if(ls[j] == OUTPUT_MONTH)
                            printf("%10d",month);
                        else if(ls[j] == OUTPUT_DATE)
                            printf(" %5d/%.02d/%.02d", hd.nzyear, month, day);
                        else if(ls[j] == OUTPUT_TIME)
                            printf(" %.02d:%.02d:%06.3f", hd.nzhour, hd.nzmin, hd.nzsec + 0.001 * hd.nzmsec);
                        else if(ls[j] == OUTPUT_PICKS) 
                        {
                            for(k = 0; k < 10; k++) 
                                printf("%12.6g", sac_header_value_float(&hd, SAC_HEADER_TMARK_POSITION + k));
                        }
                    break;
                }
            } /* endfor*/
            printf("\n");
        } /* endif */
    } /* endfor */
    return 0;
}

int sac_header_position(const char *s) 
{
    int i;
    for(i = 0; i <= SAC_HEADER_FIELDS; i++) 
    {
        if(i == 9 || i == 20) 
        {
            if(strcmp(s, SacHeaderName[i]) == 0) 
                return i;
        } 
        else if(strcasecmp(s, SacHeaderName[i]) == 0)
            return i;
    }
    if      (strcasecmp(s,"kzdate")==0)    return(OUTPUT_DATE);
    else if (strcasecmp(s,"kztime")==0)    return(OUTPUT_TIME);
    else if (strcasecmp(s,"nzday")==0)     return(OUTPUT_DAY);
    else if (strcasecmp(s,"nzmonth")==0)   return(OUTPUT_MONTH);
    else if (strcasecmp(s,"default")==0)   return(OUTPUT_DEFAULT);
    else if (strcasecmp(s,"default1")==0)  return(OUTPUT_DEFAULT1);
    else if (strcasecmp(s,"all")==0)       return(OUTPUT_ALL);
    else if (strcasecmp(s,"full")==0)      return(OUTPUT_FULL);
    else if (strcasecmp(s,"picks")==0)     return(OUTPUT_PICKS);

    return -999;
}  



/* sacio.c */
/*******************************************************************
*                       sacio.c
*       swab4           reverse byte order for integer/float
*       sac_byte_order  determine sac byte order
*       read_sachead    read SAC header
*
*********************************************************************/


/*****************************************************

  swab4

  Description:  reverse byte order for float/integer

  Author:       Lupei Zhu

  Arguments:    char *pt        pointer to byte array
                int    n        number of bytes

  Return:       none

  Modify history:
        12/03/96        Lupei Zhu       Initial coding

************************************************************/

void swab4( char *pt, int n )
{
    int i;
    char temp;
    for(i=0;i<n;i+=4) 
    {
        temp = pt[i+3];
        pt[i+3] = pt[i];
        pt[i] = temp;
        temp = pt[i+2];
        pt[i+2] = pt[i+1];
        pt[i+1] = temp;
    }
}


/***********************************************************
   
   sac_byte_order (hd)

   Description:  Determine if the header variable (internal4)
                 which defines the header version is between 0 and 6
   do_swap  = 1 Try and byte swap
            = 0 Do not byte swap
   Returns:
              0 - No swapping needs to be done
              1 - Swapping is needed
*/
int sac_byte_order( SACHEAD *hd ) 
{
    if( ( hd->nvhdr > 0  ) && ( hd->nvhdr <= 6 ) ) 
        return(0);
    return(1);
}

/***********************************************************

  read_sachead

  Description:  read binary SAC header from file.

  Author:       Lupei Zhu

  Arguments:    const char *name        file name
                SACHEAD *hd             SAC header to be filled

  Return:       0 if success, -1 if failed

  Modify history:
        05/29/97        Lupei Zhu       Initial coding
************************************************************/

int read_sachead(const char *name, SACHEAD *hd ) 
{
    FILE *strm;
  
    if ((strm = fopen(name, "rb")) == NULL) 
    {
        fprintf(stderr, "Unable to open %s\n",name);
        return -1;
    }
  
    if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) 
    {
        fprintf(stderr, "saclst: Error in reading SAC header %s\n",name);
        fclose(strm);
        return -1;
    }
    if(sac_byte_order(hd)) 
    {
        swab4((char *)hd, SAC_HEADER_SIZE_NUMBERS);
        if(sac_byte_order(hd)) 
        {
            fprintf(stderr, "saclst: Error determining SAC header: %s\n", name);
            fclose(strm);
            hd = NULL;
            return(-1);
        }
    }
  
    fclose(strm);
    return 0;
}

void kidate(int year, int jday, int *month, int *day) 
{
    int im;
    static long ndays[12]={31,28,31,30,31,30,31,31,30,31,30,31};

    if(year == -12345 || jday == -12345) 
    {
        *month = -12345;
        *day   = -12345;
        return;
    }
    if( (year / 4)*4 == year) 
        ndays[1] = 29;
    else
        ndays[1] = 28;
  
    *day = jday;
    for( *month = 1; *month <= 12; (*month)++) 
    {
        im = *month - 1;
        if(*day <= ndays[im]) 
            return;
        *day = *day - ndays[im];
    }
    *month = -12345;
    *day   = -12345;
}

/*************************************************************************
  Name:		sac.h

  Purpose:	structure for header of a SAC (Seismic Analysis Code)
  		data file, and prototype for basic SAC I/O

  Notes:	Key to comment flags describing each field:
  	Column 1:
  		R	required by SAC
    	  (blank)	optional
  	Column 2:
  		A = settable from a priori knowledge
  		D = available in data
  		F = available in or derivable from SEED fixed data header
  		T = available in SEED header tables
    	  (blank) = not directly available from SEED data, header
			tables, or elsewhere

  Problems:	none known

  References:	O'Neill, D. (1987).  IRIS Interim Data Distribution Format
                (SAC ASCII), Version 1.0 (12 November 1987).  Incorporated
		Research Institutions for Seismology, 1616 North Fort Myer
		Drive, Suite 1440, Arlington, Virginia 22209.  11 pp.
		Tull, J. (1987).  SAC User's Manual, Version 10.2, October 7,
		1987.  Lawrence Livermore National Laboratory, L-205,
		Livermore, California 94550.  ??? pp.

  Language:	C, hopefully ANSI standard

  Author:	Dennis O'Neill

  Revisions:	07/15/88  Dennis O'Neill  Initial preliminary release 0.9
                11/21/88  Dennis O'Neill  Production release 1.0
		01/27/91  Lorraine Hwang  Header number is now version 6
		07/06/93  Xiaoming Ding   structure name sac -> sac_head
		                          typedef structure to be SACHEAD
		12/06/96  Lupei Zhu	  prototype sacio functions
**************************************************************************/


#ifndef _sachead_h
#define _sachead_h

/* True/false definitions */
#ifndef TRUE
#define FALSE	0
#define TRUE	1
#endif

#define SAC_HEADER_FIELDS          133
#define SAC_HEADER_SIZE_NUMBERS    440
#define SAC_HEADER_SIZE            632 
#define SAC_HEADER_TMARK_POSITION  10
#define SAC_HEADER_USERN_POSITION  40

#define SAC_HEADER_FLOAT_MIN       0
#define SAC_HEADER_FLOAT_MAX       69
#define SAC_HEADER_INT_MIN         70
#define SAC_HEADER_INT_MAX         84
#define SAC_HEADER_ENUM_MIN        85
#define SAC_HEADER_ENUM_MAX        104
#define SAC_HEADER_LOGICAL_MIN     105
#define SAC_HEADER_LOGICAL_MAX     109
#define SAC_HEADER_CHAR_MIN        110
#define SAC_HEADER_CHAR_MAX        133
#define SAC_HEADER_CHAR_DOUBLE     111
#define SAC_HEADER_CHAR_DOUBLE_END 112

char *SacHeaderNameNull = "";

#define SAC_HEADER_FLOAT_UNDEFINED (-12345.0)
#define SAC_HEADER_INT_UNDEFINED   (-12345)
#define SAC_HEADER_CHAR_UNDEFINED  ("-12345  ")
#define SAC_HEADER_UNDEFINED       ("UNDEFINED")

typedef struct sac_head
{
  float	delta;			/* RF time increment, sec    */
  float	depmin;			/*    minimum amplitude      */
  float	depmax;			/*    maximum amplitude      */
  float	scale;			/*    amplitude scale factor */
  float	odelta;			/*    observed time inc      */
  float	b;			/* RD initial time - wrt nz* */
  float	e;			/* RD end time               */
  float	o;			/*    event start            */
  float	a;			/*    1st arrival time       */
  float	fmt;		        /*    internal use           */
  float	t0;			/*    user-defined time pick */
  float	t1;			/*    user-defined time pick */
  float	t2;			/*    user-defined time pick */
  float	t3;			/*    user-defined time pick */
  float	t4;			/*    user-defined time pick */
  float	t5;			/*    user-defined time pick */
  float	t6;			/*    user-defined time pick */
  float	t7;			/*    user-defined time pick */
  float	t8;			/*    user-defined time pick */
  float	t9;			/*    user-defined time pick */
  float	f;			/*    event end, sec > 0     */
  float	resp0;			/*    instrument respnse parm*/
  float	resp1;			/*    instrument respnse parm*/
  float	resp2;			/*    instrument respnse parm*/
  float	resp3;			/*    instrument respnse parm*/
  float	resp4;			/*    instrument respnse parm*/
  float	resp5;			/*    instrument respnse parm*/
  float	resp6;			/*    instrument respnse parm*/
  float	resp7;			/*    instrument respnse parm*/
  float	resp8;			/*    instrument respnse parm*/
  float	resp9;			/*    instrument respnse parm*/
  float	stla;			/*  T station latititude     */
  float	stlo;			/*  T station longitude      */
  float	stel;			/*  T station elevation, m   */
  float	stdp;			/*  T station depth, m       */
  float	evla;			/*    event latitude         */
  float	evlo;			/*    event longitude        */
  float	evel;			/*    event elevation        */
  float	evdp;			/*    event depth            */
  float	mag;    		/*    magnitude value        */
  float	user0;			/*    available to user      */
  float	user1;			/*    available to user      */
  float	user2;			/*    available to user      */
  float	user3;			/*    available to user      */
  float	user4;			/*    available to user      */
  float	user5;			/*    available to user      */
  float	user6;			/*    available to user      */
  float	user7;			/*    available to user      */
  float	user8;			/*    available to user      */
  float	user9;			/*    available to user      */
  float	dist;			/*    stn-event distance, km */
  float	az;			/*    event-stn azimuth      */
  float	baz;			/*    stn-event azimuth      */
  float	gcarc;			/*    stn-event dist, degrees*/
  float	sb;     		/*    saved b value          */
  float	sdelta; 		/*    saved delta value      */
  float	depmen;			/*    mean value, amplitude  */
  float	cmpaz;			/*  T component azimuth      */
  float	cmpinc;			/*  T component inclination  */
  float	xminimum;		/*    XYZ X minimum value    */
  float	xmaximum;		/*    XYZ X maximum value    */
  float	yminimum;		/*    XYZ Y minimum value    */
  float	ymaximum;		/*    XYZ Y maximum value    */
  float	unused6;		/*    reserved for future use*/
  float	unused7;		/*    reserved for future use*/
  float	unused8;		/*    reserved for future use*/
  float	unused9;		/*    reserved for future use*/
  float	unused10;		/*    reserved for future use*/
  float	unused11;		/*    reserved for future use*/
  float	unused12;		/*    reserved for future use*/
  int	nzyear;			/*  F zero time of file, yr  */
  int	nzjday;			/*  F zero time of file, day */
  int	nzhour;			/*  F zero time of file, hr  */
  int	nzmin;			/*  F zero time of file, min */
  int	nzsec;			/*  F zero time of file, sec */
  int	nzmsec;			/*  F zero time of file, msec*/
  int	nvhdr;  		/*  R header version number  */
  int	norid;  		/*    Origin ID              */
  int	nevid;  		/*    Event ID               */
  int	npts;			/* RF number of samples      */
  int	nsnpts; 		/*    saved npts             */
  int	nwfid; 		        /*    Waveform ID            */
  int	nxsize;	                /*    XYZ X size             */
  int	nysize;   		/*    XYZ Y size             */
  int	unused15;		/*    reserved for future use*/
  int	iftype;			/* RA type of file          */
  int	idep;			/*    type of amplitude      */
  int	iztype;			/*    zero time equivalence  */
  int	unused16;		/*    reserved for future use*/
  int	iinst;			/*    recording instrument   */
  int	istreg;			/*    stn geographic region  */
  int	ievreg;			/*    event geographic region*/
  int	ievtyp;			/*    event type             */
  int	iqual;			/*    quality of data        */
  int	isynth;			/*    synthetic data flag    */
  int	imagtyp;		/*    magnitude type         */
  int	imagsrc;		/*    magnitude source       */
  int	unused19;		/*    reserved for future use*/
  int	unused20;		/*    reserved for future use*/
  int	unused21;		/*    reserved for future use*/
  int	unused22;		/*    reserved for future use*/
  int	unused23;		/*    reserved for future use*/
  int	unused24;		/*    reserved for future use*/
  int	unused25;		/*    reserved for future use*/
  int	unused26;		/*    reserved for future use*/
  int	leven;			/* RA data-evenly-spaced flag*/
  int	lpspol;			/*    station polarity flag  */
  int	lovrok;			/*    overwrite permission   */
  int	lcalda;			/*    calc distance, azimuth */
  int	unused27;		/*    reserved for future use*/
  char	kstnm[8];		/*  F station name           */
  char	kevnm[16];		/*    event name             */
  char	khole[8];		/*    man-made event name    */
  char	ko[8];			/*    event origin time id   */
  char	ka[8];			/*    1st arrival time ident */
  char	kt0[8];			/*    time pick 0 ident      */
  char	kt1[8];			/*    time pick 1 ident      */
  char	kt2[8];			/*    time pick 2 ident      */
  char	kt3[8];			/*    time pick 3 ident      */
  char	kt4[8];			/*    time pick 4 ident      */
  char	kt5[8];			/*    time pick 5 ident      */
  char	kt6[8];			/*    time pick 6 ident      */
  char	kt7[8];			/*    time pick 7 ident      */
  char	kt8[8];			/*    time pick 8 ident      */
  char	kt9[8];			/*    time pick 9 ident      */
  char	kf[8];			/*    end of event ident     */
  char	kuser0[8];		/*    available to user      */
  char	kuser1[8];		/*    available to user      */
  char	kuser2[8];		/*    available to user      */
  char	kcmpnm[8];		/*  F component name         */
  char	knetwk[8];		/*    network name           */
  char	kdatrd[8];		/*    date data read         */
  char	kinst[8];		/*    instrument name        */
} SACHEAD;

/* prototype for SACIO functions */
int	read_sachead(const char *, SACHEAD *);
void	swab4(char *, int);

char *SacHeaderName[] = {
  "delta",		/* RF time increment, sec    */
  "depmin",		/*    minimum amplitude      */
  "depmax",		/*    maximum amplitude      */
  "scale",		/*    amplitude scale factor */
  "odelta",    		/*    observed time inc      */
  "b",			/* RD initial time - wrt nz* */
  "e",			/* RD end time               */
  "o",    		/*    event start            */
  "a",    		/*    1st arrival time       */
  "Fmt",  		/*    internal use           */
  
  "t0",   		/*    user-defined time pick */
  "t1",   		/*    user-defined time pick */
  "t2",   		/*    user-defined time pick */
  "t3",   		/*    user-defined time pick */
  "t4",   		/*    user-defined time pick */
  "t5",   		/*    user-defined time pick */
  "t6",   		/*    user-defined time pick */
  "t7",   		/*    user-defined time pick */
  "t8",   		/*    user-defined time pick */
  "t9",   		/*    user-defined time pick */
  
  "F",    		/*    event end, sec > 0     */
  "resp0",		/*    instrument respnse parm*/
  "resp1",		/*    instrument respnse parm*/
  "resp2",		/*    instrument respnse parm*/
  "resp3",		/*    instrument respnse parm*/
  "resp4",		/*    instrument respnse parm*/
  "resp5",		/*    instrument respnse parm*/
  "resp6",		/*    instrument respnse parm*/
  "resp7",		/*    instrument respnse parm*/
  "resp8",		/*    instrument respnse parm*/
  
  "resp9",		/*    instrument respnse parm*/
  "stla",		/*  T station latititude     */
  "stlo",		/*  T station longitude      */
  "stel",		/*  T station elevation, m   */
  "stdp",		/*  T station depth, m       */
  "evla",		/*    event latitude         */
  "evlo",		/*    event longitude        */
  "evel",		/*    event elevation        */
  "evdp",		/*    event depth            */
  "mag",	        /*    reserved for future use*/
 
  "user0",		/*    available to user      */
  "user1",		/*    available to user      */
  "user2",		/*    available to user      */
  "user3",		/*    available to user      */
  "user4",		/*    available to user      */
  "user5",		/*    available to user      */
  "user6",		/*    available to user      */
  "user7",		/*    available to user      */
  "user8",		/*    available to user      */
  "user9",		/*    available to user      */
  
  "dist",		/*    stn-event distance, km */
  "az",			/*    event-stn azimuth      */
  "baz",		/*    stn-event azimuth      */
  "gcarc",		/*    stn-event dist, degrees*/
  "sb",   		/*    internal use           */
  "sdelta", 		/*    internal use           */
  "depmen",		/*    mean value, amplitude  */
  "cmpaz",		/*  T component azimuth      */
  "cmpinc",		/*  T component inclination  */
  "xminimum",		/*    reserved for future use*/
  
  "xmaximum",		/*    reserved for future use*/
  "yminimum",		/*    reserved for future use*/
  "ymaximum",		/*    reserved for future use*/
  "unused6",		/*    reserved for future use*/
  "unused7",		/*    reserved for future use*/
  "unused8",		/*    reserved for future use*/
  "unused9",		/*    reserved for future use*/
  "unused10",		/*    reserved for future use*/
  "unused11",		/*    reserved for future use*/
  "unused12",		/*    reserved for future use*/
  
  /* ints */
  "nzyear",		/*  F zero time of file, yr  */
  "nzjday",		/*  F zero time of file, day */
  "nzhour",		/*  F zero time of file, hr  */
  "nzmin",		/*  F zero time of file, min */
  "nzsec",		/*  F zero time of file, sec */
  "nzmsec",		/*  F zero time of file, msec*/
  "nvhdr",	        /*  R header version number  */
  "norid",		/*    internal use           */
  "nevid",		/*    internal use           */
  "npts",		/* RF number of samples      */
  
  "nsnpts",		/*    internal use           */
  "nwfid",		/*    internal use           */
  "xsize",		/*    reserved for future use*/
  "ysize",		/*    reserved for future use*/
  "unused15",		/*    reserved for future use*/
  "iftype",		/* RA type of file           */
  "idep",		/*    type of amplitude      */
  "iztype",		/*    zero time equivalence  */
  "unused16",		/*    reserved for future use*/
  "iinst",		/*    recording instrument   */
  "istreg",		/*    stn geographic region  */
  "ievreg",		/*    event geographic region*/
  "ievtyp",		/*    event type             */
  "iqual",		/*    quality of data        */
  "isynth",		/*    synthetic data flag    */
  "imagtyp",      	/*    reserved for future use*/
  "imagsrc",      	/*    reserved for future use*/
  "unused19",		/*    reserved for future use*/
  "unused20",		/*    reserved for future use*/
  "unused21",		/*    reserved for future use*/
  "unused22",		/*    reserved for future use*/
  "unused23",		/*    reserved for future use*/
  "unused24",		/*    reserved for future use*/
  "unused25",		/*    reserved for future use*/
  "unused26",		/*    reserved for future use*/
  "leven",		/* RA data-evenly-spaced flag*/
  "lpspol",		/*    station polarity flag  */
  "lovrok",		/*    overwrite permission   */
  "lcalda",		/*    calc distance, azimuth */
  "unused27",		/*    reserved for future use*/
  "kstnm",		/*  F station name           */
  "kevnm",		/*    event name             */
  "kevnm empty",        /*                           */
  "khole",		/*    man-made event name    */
  "ko",			/*    event origin time id   */
  "ka",			/*    1st arrival time ident */
  "kt0",		/*    time pick 0 ident      */
  "kt1",		/*    time pick 1 ident      */
  "kt2",		/*    time pick 2 ident      */
  "kt3",		/*    time pick 3 ident      */
  "kt4",		/*    time pick 4 ident      */
  "kt5",		/*    time pick 5 ident      */
  "kt6",		/*    time pick 6 ident      */
  "kt7",		/*    time pick 7 ident      */
  "kt8",		/*    time pick 8 ident      */
  "kt9",		/*    time pick 9 ident      */
  "kf",			/*    end of event ident     */
  "kuser0",		/*    available to user      */
  "kuser1",		/*    available to user      */
  "kuser2",		/*    available to user      */
  "kcmpnm",		/*  F component name         */
  "knetwk",		/*    network name           */
  "kdatrd",		/*    date data read         */
  "kinst"               /*    instrument name        */
};

enum SAC_HEADER_TYPES {
  SAC_HEADER_UNDEFINED_TYPE = 1, 
  SAC_HEADER_FLOAT_TYPE,   SAC_HEADER_INT_TYPE,   SAC_HEADER_ENUM_TYPE,
  SAC_HEADER_LOGICAL_TYPE, SAC_HEADER_CHAR8_TYPE, SAC_HEADER_CHAR16_TYPE
};
  
enum SAC_HEADER_ORDER {
  SAC_HEADER_DELTA = 0, 
  SAC_HEADER_DEPMIN,   SAC_HEADER_DEPMAX,    SAC_HEADER_SCALE,    SAC_HEADER_ODELTA,   SAC_HEADER_B,
  SAC_HEADER_E,        SAC_HEADER_O,         SAC_HEADER_A,        SAC_HEADER_FMT,      SAC_HEADER_T0,       
  SAC_HEADER_T1,       SAC_HEADER_T2,        SAC_HEADER_T3,       SAC_HEADER_T4,       SAC_HEADER_T5,       
  SAC_HEADER_T6,       SAC_HEADER_T7,        SAC_HEADER_T8,       SAC_HEADER_T9,       SAC_HEADER_F,        
  SAC_HEADER_RESP0,    SAC_HEADER_RESP1,     SAC_HEADER_RESP2,    SAC_HEADER_RESP3,    SAC_HEADER_RESP4,    
  SAC_HEADER_RESP5,    SAC_HEADER_RESP6,     SAC_HEADER_RESP7,    SAC_HEADER_RESP8,    SAC_HEADER_RESP9,    
  SAC_HEADER_STLA,     SAC_HEADER_STLO,      SAC_HEADER_STEL,     SAC_HEADER_STDP,     SAC_HEADER_EVLA,     
  SAC_HEADER_EVLO,     SAC_HEADER_EVEL,      SAC_HEADER_EVDP,     SAC_HEADER_MAG,      SAC_HEADER_USER0,    
  SAC_HEADER_USER1,    SAC_HEADER_USER2,     SAC_HEADER_USER3,    SAC_HEADER_USER4,    SAC_HEADER_USER5,    
  SAC_HEADER_USER6,    SAC_HEADER_USER7,     SAC_HEADER_USER8,    SAC_HEADER_USER9,    SAC_HEADER_DIST,     
  SAC_HEADER_AZ,       SAC_HEADER_BAZ,       SAC_HEADER_GCARC,    SAC_HEADER_SB,       SAC_HEADER_SDELTA,  
  SAC_HEADER_DEPMEN,   SAC_HEADER_CMPAZ,     SAC_HEADER_CMPINC,   SAC_HEADER_XMINIMUM, SAC_HEADER_XMAXIMUM, 
  SAC_HEADER_YMINIMUM, SAC_HEADER_YMAXIMUM,  SAC_HEADER_UNUSED6,  SAC_HEADER_UNUSED7,  SAC_HEADER_UNUSED8,  
  SAC_HEADER_UNUSED9,  SAC_HEADER_UNUSED10,  SAC_HEADER_UNUSED11, SAC_HEADER_UNUSED12, SAC_HEADER_NZYEAR,
  SAC_HEADER_NZJDAY,   SAC_HEADER_NZHOUR,    SAC_HEADER_NZMIN,    SAC_HEADER_NZSEC,    SAC_HEADER_NSMSEC,
  SAC_HEADER_NVHDR,    SAC_HEADER_NORID,     SAC_HEADER_NEVID,    SAC_HEADER_NPTS,     SAC_HEADER_NSNPTS,
  SAC_HEADER_NWFID,    SAC_HEADER_NXSIZE,    SAC_HEADER_NYSIZE,   SAC_HEADER_UNUSED15, SAC_HEADER_IFTYPE,
  SAC_HEADER_IDEP,     SAC_HEADER_IZTYPE,    SAC_HEADER_UNUSED16, SAC_HEADER_IINST,    SAC_HEADER_ISTREG, 
  SAC_HEADER_IEVREG,   SAC_HEADER_IEVTYP,    SAC_HEADER_IQUAL,    SAC_HEADER_ISYNTH,   SAC_HEADER_IMAGTYPE, 
  SAC_HEADER_IMAGSRC,  SAC_HEADER_UNUSED19,  SAC_HEADER_UNUSED20, SAC_HEADER_UNUSED21, SAC_HEADER_UNUSED22, 
  SAC_HEADER_UNUSED23, SAC_HEADER_UNUSED24,  SAC_HEADER_UNUSED25, SAC_HEADER_UNUSED26, SAC_HEADER_LEVEN,
  SAC_HEADER_LPSPOL,   SAC_HEADER_LOVROK,    SAC_HEADER_LCALDA,   SAC_HEADER_UNUSED27, SAC_HEADER_KSTNM,
  SAC_HEADER_KEVNM,    SAC_HEADER_KEVNM_END, SAC_HEADER_KHOLE,    SAC_HEADER_KO,       SAC_HEADER_KA,
  SAC_HEADER_KT0,      SAC_HEADER_KT1,       SAC_HEADER_KT2,      SAC_HEADER_KT3,      SAC_HEADER_KT4,
  SAC_HEADER_KT5,      SAC_HEADER_KT6,       SAC_HEADER_KT7,      SAC_HEADER_KT8,      SAC_HEADER_KT9,
  SAC_HEADER_KF,       SAC_HEADER_KUSER0,    SAC_HEADER_KUSER1,   SAC_HEADER_KUSER2,   SAC_HEADER_KCMPNM,
  SAC_HEADER_KNETWK,   SAC_HEADER_KDATRD,    SAC_HEADER_KINST
};

enum SAC_HEADER_ENUMS {
  /* enumerated header values */
  IREAL    = 0,  /* Undocumented                */
  ITIME    = 1,  /* Time series file            */
  IRLIM    = 2,  /* Spectral file-real/imag     */
  IAMPH    = 3,  /* Spectral file-ampl/phase    */
  IXY      = 4,  /* General x vs y file         */
  IUNKN    = 5,  /* Unknown                     */
  IDISP    = 6,  /* Displacement (NM)           */
  IVEL     = 7,  /* Velocity (NM/SEC)           */
  IACC     = 8,  /* Acceleration (CM/SEC/SEC)   */
  IB       = 9,  /* Begin time                  */
  IDAY     = 10,  /* GMT day                     */
  IO       = 11,  /* Event origin time           */
  IA       = 12,  /* First arrival time          */
  IT0      = 13,  /* User defined time pick 0    */
  IT1      = 14,  /* User defined time pick 1    */
  IT2      = 15,  /* User defined time pick 2    */
  IT3      = 16,  /* User defined time pick 3    */
  IT4      = 17,  /* User defined time pick 4    */
  IT5      = 18,  /* User defined time pick 5    */
  IT6      = 19,  /* User defined time pick 6    */
  IT7      = 20,  /* User defined time pick 7    */
  IT8      = 21,  /* User defined time pick 8    */
  IT9      = 22,  /* User defined time pick 9    */
  IRADNV   = 23,  /* Radial (NTS)                */
  ITANNV   = 24,  /* Tangential (NTS)            */
  IRADEV   = 25,  /* Radial (EVENT)              */
  ITANEV   = 26,  /* Tangential (EVENT)          */
  INORTH   = 27,  /* North positive              */
  IEAST    = 28,  /* East positive               */
  IHORZA   = 29,  /* Horizontal (ARB)            */
  IDOWN    = 30,  /* Down positive               */
  IUP      = 31,  /* Up positive                 */
  ILLLBB   = 32,  /* LLL broadband               */
  IWWSN1   = 33,  /* WWSN 15-100                 */
  IWWSN2   = 34,  /* WWSN 30-100                 */
  IHGLP    = 35,  /* High-gain long-period       */
  ISRO     = 36,  /* SRO                         */
  INUCL    = 37,  /* Nuclear event               */
  IPREN    = 38,  /* Nuclear pre-shot event      */
  IPOSTN   = 39,  /* Nuclear post-shot event     */
  IQUAKE   = 40,  /* Earthquake                  */
  IPREQ    = 41,  /* Foreshock                   */
  IPOSTQ   = 42,  /* Aftershock                  */
  ICHEM    = 43,  /* Chemical explosion          */
  IOTHER   = 44,  /* Other                       */
  IGOOD    = 45,  /* Good                        */
  IGLCH    = 46,  /* Gliches                     */
  IDROP    = 47,  /* Dropouts                    */
  ILOWSN   = 48,  /* Low signal to noise ratio   */
  IRLDTA   = 49,  /* Real data                   */
  IVOLTS   = 50,  /* Velocity (volts)            */
  IXYZ     = 51,  /* General XYZ (3-D) file      */
  /* These 18 added to describe magnitude type and source maf 970205 */
  IMB      = 52,  /* Bodywave Magnitude          */
  IMS      = 53,  /* Surface Magnitude           */
  IML      = 54,  /* Local Magnitude             */
  IMW      = 55,  /* Moment Magnitude            */
  IMD      = 56,  /* Duration Magnitude          */
  IMX      = 57,  /* User Defined Magnitude      */
  INEIC    = 58,  /* INEIC                       */
  IPDEQ    = 59,  /* IPDEQ                       */
  IPDEW    = 60,  /* IPDEW                       */
  IPDE     = 61,  /* IPDE                        */
  IISC     = 62,  /* IISC                        */
  IREB     = 63,  /* IREB                        */
  IUSGS    = 64,  /* IUSGS                       */
  IBRK     = 65,  /* IBRK                        */
  ICALTECH = 66,  /* ICALTECH                    */
  ILLNL    = 67,  /* ILLNL                       */
  IEVLOC   = 68,  /* IEVLOC                      */
  IJSOP    = 69,  /* IJSOP                       */
  IUSER    = 70,  /* IUSER                       */
  IUNKNOWN = 71,  /* IUNKNOWN                    */
  /* These  17 added for ievtyp. maf 970325 */
  IQB	   = 72,  /* Quarry or mine blast confirmed by quarry */
  IQB1	   = 73,  /* Quarry or mine blast with designed shot information-ripple fired */
  IQB2	   = 74,  /* Quarry or mine blast with observed shot information-ripple fired */
  IQBX     = 75,  /* Quarry or mine blast - single shot */
  IQMT     = 76,  /* Quarry or mining-induced events: tremors and rockbursts */
  IEQ      = 77,  /* Earthquake                  */
  IEQ1     = 78,  /* Earthquakes in a swarm or aftershock sequence */
  IEQ2     = 79,  /* Felt earthquake             */
  IME      = 80,  /* Marine explosion            */
  IEX	   = 81,  /* Other explosion             */
  INU	   = 82,  /* Nuclear explosion           */
  INC	   = 83,  /* Nuclear cavity collapse     */
  IO_	   = 84,  /* Other source of known origin */
  IL	   = 85,  /* Local event of unknown origin */
  IR	   = 86,  /* Regional event of unknown origin */
  IT	   = 87,  /* Teleseismic event of unknown origin */
  IU	   = 88,  /* Undetermined or conflicting information  */
  /* These 9 added for ievtype to keep up with database. maf 000530 */
  IEQ3     = 89,  /* Damaging Earthquake         */
  IEQ0     = 90,  /* Probable earthquake         */
  IEX0     = 91,  /* Probable explosion          */
  IQC      = 92,  /* Mine collapse               */
  IQB0     = 93,  /* Probable Mine Blast         */
  IGEY     = 94,  /* Geyser                      */
  ILIT     = 95,  /* Light                       */
  IMET     = 96,  /* Meteroic event              */
  IODOR    = 97   /* Odors                       */
};

char *SacHeaderEnums[] = {
  "IREAL",     /* 0    To be consistent with defines above */
  /* iftype */
  "ITIME",        /* 1    Time series file            */
  "IRLIM",        /* 2    Spectral file-real/imag     */
  "IAMPH",        /* 3    Spectral file-ampl/phase    */
  "IXY",          /* 4    General x vs y file         */
  "IUNKN",        /* 5    Unknown                     */

  /* idep */
  "IDISP",        /* 6    Displacement (NM)           */
  "IVEL",         /* 7    Velocity (NM/SEC)           */
  "IACC",         /* 8    Acceleration (CM/SEC/SEC)   */

  /* iztype */
  "IB",           /* 9    Begin time                  */
  "IDAY",         /* 10   GMT day                     */
  "IO",           /* 11   Event origin time           */
  "IA",           /* 12   First arrival time          */
  "IT0",          /* 13   User defined time pick 0    */
  "IT1",          /* 14   User defined time pick 1    */
  "IT2",          /* 15   User defined time pick 2    */
  "IT3",          /* 16   User defined time pick 3    */
  "IT4",          /* 17   User defined time pick 4    */
  "IT5",          /* 18   User defined time pick 5    */
  "IT6",          /* 19   User defined time pick 6    */
  "IT7",          /* 20   User defined time pick 7    */
  "IT8",          /* 21   User defined time pick 8    */
  "IT9",          /* 22   User defined time pick 9    */
  
  /* iinst */
  "IRADNV",	/* 23   Radial (NTS)                */
  "ITANNV",	/* 24   Tangential (NTS)            */
  "IRADEV",	/* 25   Radial (EVENT)              */
  "ITANEV",	/* 26   Tangential (EVENT)          */
  "INORTH",	/* 27   North positive              */
  "IEAST",	/* 28   East positive               */
  "IHORZA",	/* 29   Horizontal (ARB)            */
  "IDOWN",	/* 30   Down positive               */
  "IUP",	/* 31   Up positive                 */
  "ILLLBB",	/* 32   LLL broadband               */
  "IWWSN1",	/* 33   WWSN 15-100                 */
  "IWWSN2",	/* 34   WWSN 30-100                 */
  "IHGLP",	/* 35   High-gain long-period       */
  "ISRO",	/* 36   SRO                         */

  /* ievtyp */
  "INUCL",	/* 37   Nuclear event               */
  "IPREN",	/* 38   Nuclear pre-shot event      */
  "IPOSTN",	/* 39   Nuclear post-shot event     */
  "IQUAKE",	/* 40   Earthquake                  */
  "IPREQ",	/* 41   Foreshock                   */
  "IPOSTQ",	/* 42   Aftershock                  */
  "ICHEM",	/* 43   Chemical explosion          */
  "IOTHER",	/* 44   Other                       */

  /* iqual */
  "IGOOD",	/* 45   Good                        */
  "IGLCH",	/* 46   Gliches                     */
  "IDROP",	/* 47   Dropouts                    */
  "ILOWSN",	/* 48   Low signal to noise ratio   */

  /* isynth */
  "IRLDTA",	/* 49   Real data                   */
  "IVOLTS",	/* 50   Velocity (volts)            */
  "IXYZ",	/* 51   General XYZ (3-D) file      */

  /* These 18 added to describe magnitude type and source maf 970205 */
  "IMB",        /* 52   Bodywave Magnitude */  
  "IMS",        /* 53   Surface Magnitude */   
  "IML",        /* 54   Local Magnitude  */    
  "IMW",        /* 55   Moment Magnitude */
  "IMD",        /* 56   Duration Magnitude */
  "IMX",        /* 57   User Defined Magnitude */
  "INEIC",	/* 58   INEIC */
  "IPDEQ",	/* 59   IPDEQ */
  "IPDEW",	/* 60   IPDEW */
  "IPDE",       /* 61   IPDE */
  "IISC",       /* 62   IISC */
  "IREB",       /* 63   IREB */
  "IUSGS",	/* 64   IUSGS */
  "IBRK",       /* 65   IBRK */
  "ICALTECH",	/* 66   ICALTECH */
  "ILLNL",	/* 67   ILLNL */
  "IEVLOC",	/* 68   IEVLOC */
  "IJSOP",	/* 69   IJSOP */
  "IUSER",	/* 70   IUSER */
  "IUNKNOWN",	/* 71   IUNKNOWN */
  
  /*   These 17 added for ievtyp. maf 970325 */
  "IQB",        /* 72   Quarry or mine blast confirmed by quarry */
  "IQB1",       /* 73   Quarry or mine blast with designed shot information-ripple fired*/
  "IQB2",       /* 74   Quarry or mine blast with observed shot information-ripple fired*/
  "IQBX",       /* 75   Quarry or mine blast - single shot */
  "IQMT",       /* 76   Quarry or mining-induced events: tremors and rockbursts */
  "IEQ",        /* 77   Earthquake */
  "IEQ1",       /* 78   Earthquakes in a swarm or aftershock sequence */
  "IEQ2",       /* 79   Felt earthquake */
  "IME",       	/* 80   Marine explosion */
  "IEX",        /* 81   Other explosion */
  "INU",        /* 82   Nuclear explosion */
  "INC",        /* 83   Nuclear cavity collapse */
  "IO_",        /* 84   Other source of known origin */
  "IL",	        /* 85   Local event of unknown origin */
  "IR",	        /* 86   Regional event of unknown origin */
  "IT",	        /* 87   Teleseismic event of unknown origin */
  "IU",	        /* 88   Undetermined or conflicting information  */
  
  /*   These 9 added for ievtype to keep up with database. maf 000530 */
  "IEQ3",       /* 89   Damaging Earthquake */
  "IEQ0",       /* 90   Probable earthquake */
  "IEX0",       /* 91   Probable explosion */
  "IQC",        /* 92   Mine collapse */
  "IQB0",       /* 93   Probable Mine Blast */
  "IGEY",       /* 94   Geyser */
  "ILIT",       /* 95   Light */
  "IMET",       /* 96   Meteroic event */
  "IODOR"       /* 97   Odors */
};

const int SacHeaderEnumsLength = sizeof(SacHeaderEnums) / sizeof(char *);

char *SacHeaderEnumsDescription[] = {
  "Undocumented", 
  /* iftype */
  "Time Series File",   "Spectral File-Real/Imag",  "Spectral File-Ampl/Phase",
  "General X vs Y file",  "Unknown", 

  /* idep */
  "Displacement (nm)",  "Velocity (nm/sec)",  "Acceleration (cm/sec/sec)",

  /* iztype */
  "Begin Time",  "GMT Day",   "Event Origin Time",  "First Arrival Time",
  "User Defined Time Pick 0",  "User Defined Time Pick 1",
  "User Defined Time Pick 2",  "User Defined Time Pick 3",
  "User Defined Time Pick 4",  "User Defined Time Pick 5",
  "User Defined Time Pick 6",  "User Defined Time Pick 7",
  "User Defined Time Pick 8",  "User Defined Time Pick 9",
  
  /* iinst */
  "Radial (NTS)",  "Tangential (NTS)",  "Radial (Event)",  "Tangential (Event)",
  "North Positive",  "East Positive",  "Horizontal (ARB)",  "Down Positive",
  "Up Positive",  "LLL Broadband",  "WWSN 15-100",  "WWSN 30-100",
  "High Gain Long Period",  "SRO",

  /* ievtyp */
  "Nuclear Event",  "Nuclear Pre-Shot Event",  "Nuclear Post-Shot Event",
  "Earthquake",  "Foreshock",  "Aftershock",  "Chemical Explosion",
  "Other",

  /* iqual */
  "Good",  "Glitches",  "Dropouts",  "Low Signal to Noise Ratio",

  /* isynth */
  "Real Data",  "Velocity (Volts)",  "General XYZ (3-D) file",

  /* These 18 added to describe magnitude type and source maf 970205 */
  "Body Wave Magnitude (mb)",  "Surface Wave Magnitude (Ms)",
  "Local Magnitude (ML)",  "Moment Magnitude (Mw)",
  "Duration Magnitude (Md)",  "User Defined Magnitude",
  "NEIC",  "PDEQ",  "PDEW",  "PDE",  "ISC",  "REB",  "USGS",  "Berkeley",
  "Caltech",  "LLNL",  "EVLOC",  "JSOP",  "User",  "Unknown",
  
  /*   These 17 added for ievtyp. maf 970325 */
  "Quarry/Mine Blast, Confirmed by Quarry", 
  "Quarry/Mine Blast with Shot Information, Ripple Fired", 
  "Quarry/Mine Blast with Observed Shot Information, Ripple Fired",
  "Quarry/Mine Blast, Single Shot",
  "Quarry or Mining Induced Events, Tremors and Rockbursts",
  "Earthquake",  "Earthquake, Swarm or Aftershock Sequence",
  "Earthquake, Felt",  "Marine Explosion",
  "Other Explosion",  "Nuclear Explosion",
  "Nuclear Cavity Collapse",   "Other Source, Unknown Origin",
  "Local Event, Unknown Origin",  "Regional Event, Unknown Origin",
  "Teleseismic Event, Unknown Origin",	
  "Undetermined or Conflicting Information",
  
  /*   These 9 added for ievtype to keep up with database. maf 000530 */
  "Damaging Earthquake",  "Probable Earthquake",  "Probable Explosion",
  "Mine Collapse",      "Probable Mine Blast",   "Geyser",    
  "Light",       "Meteroic Event",   "Odors" 
};

#endif

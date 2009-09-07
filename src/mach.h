/* ../../inc/mach.h */

#ifndef _MACH_H
#define _MACH_H

#ifdef PI
#undef PI
#endif
#define	FMAXFLOAT	3.40282e38
#define	FMINFLOAT	(-FMAXFLOAT)
#define	KDIRDL	'/'
#define	KSUBDL	'/'
#define	KWCONC	"[]"
#define	KWMULT	'*'
#define	KWSNGL	'?'
#define	MCMSG	1001
#define	MCPFN	__FSIZE__	/* max length of filename maybe. maf 960619 */
#define	MCPW	8
#define	MDFL	1000		/* max number of files in memory? */
#define MAXCHARS ( MCPFN * MDFL ) /* max number of chars in list of files in mem. maf 970806 */
/* #define MFILELIST MDFL*MCPFN/2 */
#define	MLARGE	2147483647
#define	MODEFILECASE	0
#define	MUNINP	stdin
#define	MUNOUT	stdout
#define munout  MUNOUT
#define	PI	((float)(3.141592654))
#define	RNDOFF	1.0e-06
#define	TODEG	57.29577950
#define	TORAD	(1./TODEG)
#define	VLARGE	3.40282e38
#define	VSMALL	1.0e-30


#define lge2(a,b,n) (memcmp(a,b,n) >= 0)
#define lgt2(a,b,n) (memcmp(a,b,n) >  0)
#define lle2(a,b,n) (memcmp(a,b,n) <= 0)
#define llt2(a,b,n) (memcmp(a,b,n) <  0)

#endif /* _MACH_H */

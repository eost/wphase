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

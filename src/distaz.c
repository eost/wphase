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

/*      Distance, azimuth, backazimuth calculation        */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "complex.h"
#include "mach.h"

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

double powi(b,x)
    double b;
    long x;
{
    double temp;
    long i;

    if ( b == 0.0 ) 
        return( (double) 0.0 );
    if ( x == 0L ) 
        return( (double) 1.0 ) ;

    if ( x > 0L ) 
    {
        temp = b;
        for ( i = x-1; i > 0L; i-- ) temp *= b;
        return temp;
    }

    if ( x < 0L ) 
    {
        temp = 1.0 / b;
        for ( i = x+1; i < 0L; i++ ) temp *= (1.0/b);
        return temp;
    }

  return 0 ;
}


void /*FUNCTION*/ distaz(double the, double phe, float *ths, float *phs, int ns, float *dist, float *az, float *baz, float *xdeg, int *nerr)
{
    long laz, lbaz, ldist, lxdeg;
    long int idx;
    double a, a1, a12, a12bot, a12top, al, b, b0, b1, c, c0, c1, c2, 
          c4, cosa12, costhi, costhk, d, d1, dl, du, e, e1, e1p1, e2, e3, 
          ec2, el, f, f1, g, g1, h, h1, onemec2, p1, p2, pdist, pherad, 
          phsrad, sc, sd, sina12, sinthi, sinthk, sqrte1p1, ss, t1, t2, 
          tanthi, tanthk, temp, therad, thg, thsrad, u1, u1bot, u2, u2bot, 
          u2top, v1, v2, x2, y2, z1, z2;//eps
    static float rad = 6378.160;
    static float fl = 0.00335293;
    static float twopideg = 360.;
    static float c00 = 1.;
    static float c01 = 0.25;
    static float c02 = -4.6875e-02;
    static float c03 = 1.953125e-02;
    static float c21 = -0.125;
    static float c22 = 3.125e-02;
    static float c23 = -1.46484375e-02;
    static float c42 = -3.90625e-03;
    static float c43 = 2.9296875e-03;
    static float degtokm = 111.3199;

    float *const Az = &az[0] - 1;
    float *const Baz = &baz[0] - 1;
    float *const Dist = &dist[0] - 1;
    float *const Phs = &phs[0] - 1;
    float *const Ths = &ths[0] - 1;
    float *const Xdeg = &xdeg[0] - 1;

    /*=====================================================================
     * PURPOSE:  To compute the distance and azimuth between locations.
     *=====================================================================
     * INPUT ARGUMENTS:
     *    THE:     Event latitude in decimal degrees, North positive. [r]
     *    PHE:     Event longitude, East positive. [r]
     *    THS:     Array of station latitudes. [r]
     *    PHS:     Array of station longitudes. [r]
     *    NS:      Length of THS and PHS. [i]
     *=====================================================================
     * OUTPUT ARGUMENTS:
     *    DIST:    Array of epicentral distances in km. [r]
     *    AZ:      Array of azimuths in degrees. [r]
     *    BAZ:     Array of back azimuths. [r]
     *    XDEG:    Array of great circle arc lengths. [r]
     *    NERR:    Error flag:
     *             =    0   No error.
     *             = 0904   Calculation failed internal consistency checks.
     *=====================================================================
     * MODULE/LEVEL:  DFM/4
     *=====================================================================
     * GLOBAL INPUT:
     *    MACH:
     *=====================================================================
     * SUBROUTINES CALLED:
     *    SACLIB:  SETMSG, APCMSG
     *=====================================================================
     * LOCAL VARIABLES:
     *=====================================================================
     * KNOWN ERRORS:
     * - Problem with equation for distance. See discussion below.
     *===================================================================== */
    /* PROCEDURE: */
    /* - Calculations are based upon the reference spheroid of 1968 and
     *   are defined by the major radius (RAD) and the flattening (FL). */
    /* - Initialize. */
    *nerr = 0;
    ec2 = 2.*fl - fl*fl;
    onemec2 = 1. - ec2;
    //eps = 1. + ec2/onemec2;

    /* - Check which output items are required. */

    laz = TRUE;
    if( Az[1] < 0. )
        laz = FALSE;
    lbaz = TRUE;
    if( Baz[1] < 0. )
        lbaz = FALSE;
    ldist = TRUE;
    if( Dist[1] < 0. )
        ldist = FALSE;
    lxdeg = TRUE;
    if( Xdeg[1] < 0. )
        lxdeg = FALSE;

    /* - Convert event location to radians.
     *   (Equations are unstable for latidudes of exactly 0 degrees.) */

    temp = the;
    if( temp == 0. )
        temp = 1.0e-08;
    therad = TORAD*temp;
    pherad = TORAD*phe;

    /* - Must convert from geographic to geocentric coordinates in order
     *   to use the spherical trig equations.  This requires a latitude
     *   correction given by: 1-EC2=1-2*FL+FL*FL */

    if ( the == 90 || the == -90 )  /* special attention at the poles */
        thg = the*TORAD ;           /* ... to avoid division by zero. */
    else
        thg = atan( onemec2*tan( therad ) );

    d = sin( pherad );
    e = -cos( pherad );
    f = -cos( thg );
    c = sin( thg );
    a = f*e;
    b = -f*d;
    g = -c*e;
    h = c*d;

    /* - Loop on stations: */

    for( idx = 1; idx <= ns; idx++ )
    {
        /* -- Convert to radians. */
        temp = Ths[idx];
        if( temp == 0. )
            temp = 1.0e-08;
        thsrad = TORAD*temp;
        phsrad = TORAD*Phs[idx];

        /* -- Calculate some trig constants. */
        if ( Ths[idx] == 90 || Ths[idx] == -90 )
            thg = Ths[idx] * TORAD ;
        else
            thg = atan( onemec2*tan( thsrad ) );
        d1 = sin( phsrad );
        e1 = -cos( phsrad );
        f1 = -cos( thg );
        c1 = sin( thg );
        a1 = f1*e1;
        b1 = -f1*d1;
        g1 = -c1*e1;
        h1 = c1*d1;
        sc = a*a1 + b*b1 + c*c1;

        /* - Spherical trig relationships used to compute angles. */
        if( lxdeg )
        {
            sd = 0.5*sqrt( (powi(a - a1,2) + powi(b - b1,2) + powi(c - 
                 c1,2))*(powi(a + a1,2) + powi(b + b1,2) + powi(c + c1,2)) );
            Xdeg[idx] = atan2( sd, sc )*TODEG;
            if( Xdeg[idx] < 0. )
                Xdeg[idx] = Xdeg[idx] + twopideg;
        }
        if( laz )
        {
            ss = powi(a1 - d,2) + powi(b1 - e,2) + powi(c1,2) - 2.;
            sc = powi(a1 - g,2) + powi(b1 - h,2) + powi(c1 - f,2) - 2.;
            Az[idx] = atan2( ss, sc )*TODEG;
            if( Az[idx] < 0. )
                Az[idx] = Az[idx] + twopideg;
        }
        if( lbaz )
        {
            ss = powi(a - d1,2) + powi(b - e1,2) + powi(c,2) - 2.;
            sc = powi(a - g1,2) + powi(b - h1,2) + powi(c - f1,2) - 2.;
            Baz[idx] = atan2( ss, sc )*TODEG;
            if( Baz[idx] < 0. )
                Baz[idx] = Baz[idx] + twopideg;
        }

        /* - Now compute the distance between the two points using Rudoe's
         *   formula given in GEODESY, section 2.15(b).
         *   (There is some numerical problem with the following formulae.
         *   If the station is in the southern hemisphere and the event in
         *   in the northern, these equations give the longer, not the
         *   shorter distance between the two locations.  Since the
         *   equations are fairly messy, the simplist solution is to reverse
         *   the meanings of the two locations for this case.) */
        if( ldist )
        {
            if( thsrad > 0. )
            {
                t1 = thsrad;
                p1 = phsrad;
                t2 = therad;
                p2 = pherad;

                /* special attention at the poles to avoid atan2 troubles 
                   and division by zero. */
                if( the == 90.0 ) 
                {
                    costhk = 0.0 ;
                    sinthk = 1.0 ;
                    tanthk = FLT_MAX ;
                }
                else if ( the == -90.0 ) 
                {
                    costhk = 0.0 ;
                    sinthk = -1.0 ;
                    tanthk = -FLT_MAX ;
                }
                else
                {
                    costhk = cos( t2 );
                    sinthk = sin( t2 );
                    tanthk = sinthk/costhk;
                }

                /* special attention at the poles continued. */
                if ( Ths[idx] == 90.0 ) 
                {
                    costhi = 0.0 ;
                    sinthi = 1.0 ;
                    tanthi = FLT_MAX ;
                }
                else if ( Ths[idx] == -90.0 ) 
                {
                    costhi = 0.0 ;
                    sinthi = -1.0 ;
                    tanthi = -FLT_MAX ;
                }
                else 
                {
                    costhi = cos( t1 );
                    sinthi = sin( t1 );
                    tanthi = sinthi/costhi;
                } 
            }
            else
            {
                t1 = therad;
                p1 = pherad;
                t2 = thsrad;
                p2 = phsrad;

                /* more special attention at the poles */
                if ( Ths[idx] == 90.0 ) 
                {
                    costhk = 0.0 ;
                    sinthk = 1.0 ;
                    tanthk = FLT_MAX ;
                }
                else if ( Ths[idx] == -90.0 ) 
                {
                    costhk = 0.0 ;
                    sinthk = -1.0 ;
                    tanthk = -FLT_MAX ;
                }
                else 
                {
                    costhk = cos( t2 );
                    sinthk = sin( t2 );
                    tanthk = sinthk/costhk;
                } 

                /* more special attention at the poles continued */
                if ( the == 90.0 ) 
                {
                    costhi = 0.0 ;
                    sinthi = 1.0 ;
                    tanthi = FLT_MAX ;
                }
                else if ( the == -90.0 ) 
                {
                    costhi = 0.0 ;
                    sinthi = -1.0 ;
                    tanthi = -FLT_MAX ;
                }
                else 
                {
                    costhi = cos( t1 );
                    sinthi = sin( t1 );
                    tanthi = sinthi/costhi;
                }

            }

            el = ec2/onemec2;
            e1 = 1. + el;
            al = tanthi/(e1*tanthk) + ec2*sqrt( (e1 + powi(tanthi,2))/
                 (e1 + powi(tanthk,2)) );
            dl = p1 - p2;
            a12top = sin( dl );
            a12bot = (al - cos( dl ))*sinthk;

            /* Rewrote these three lines with help from trig identities.  maf 990415 */
            a12 = atan2( a12top, a12bot );
            cosa12 = cos( a12 );
            sina12 = sin( a12 );

            /*cosa12 = sqrt ( a12bot*a12bot / ( a12bot*a12bot + a12top*a12top ) ) ;
              sina12 = sqrt ( a12top*a12top / ( a12bot*a12bot + a12top*a12top ) ) ; */

            e1 = el*(powi(costhk*cosa12,2) + powi(sinthk,2));
            e2 = e1*e1;
            e3 = e1*e2;
            c0 = c00 + c01*e1 + c02*e2 + c03*e3;
            c2 = c21*e1 + c22*e2 + c23*e3;
            c4 = c42*e2 + c43*e3;
            v1 = rad/sqrt( 1. - ec2*powi(sinthk,2) );
            v2 = rad/sqrt( 1. - ec2*powi(sinthi,2) );
            z1 = v1*(1. - ec2)*sinthk;
            z2 = v2*(1. - ec2)*sinthi;
            x2 = v2*costhi*cos( dl );
            y2 = v2*costhi*sin( dl );
            e1p1 = e1 + 1.;
            sqrte1p1 = sqrt( e1p1 );
            u1bot = sqrte1p1*cosa12;
            u1 = atan2( tanthk, u1bot );
            u2top = v1*sinthk + e1p1*(z2 - z1);
            u2bot = sqrte1p1*(x2*cosa12 - y2*sinthk*sina12);
            u2 = atan2( u2top, u2bot );
            b0 = v1*sqrt( 1. + el*powi(costhk*cosa12,2) )/e1p1;
            du = u2 - u1;
            pdist = b0*(c2*(sin( 2.*u2 ) - sin( 2.*u1 )) + c4*(sin( 4.*
            u2 ) - sin( 4.*u1 )));
            Dist[idx] = fabs( b0*c0*du + pdist );
            if( lxdeg && (fabs( Dist[idx] - degtokm*Xdeg[idx] )) > 100. )
            {
                *nerr = 904;
                printf("Error %d\n", *nerr);
                /*                  
                    setmsg( "ERROR", *nerr );
                    apimsg( idx );
                */
            }
        } /* end if ( ldist ) */
    } /* end for */

//L_8888:
    return;

        /*=====================================================================
         * MODIFICATION HISTORY:
         *    830603:  Fixed bug with negative station latiudes.
         *    810000:  Original version.
         *===================================================================== */

} /* end of function */


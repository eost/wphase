#!/usr/bin/env python
# *-* coding: iso-8859-1 *-*

import scipy as sp
import pylab as pyl

def main(MT):
    N = 200
    x = sp.arange(-N,N+1,1)
    y = sp.arange(-N,N+1,1)
    fig = sp.zeros((2*N+1,2*N+1),dtype='i')
    for i in x:     # North
        for j in y: # East
            rad = sp.sqrt(sp.cast['g'](i*i + j*j))
            if rad < N:
                Ain = sp.cast['g'](2.) * sp.arcsin(rad/sp.cast['g'](N)/sp.sqrt(2))
                Azi = atan2(i,j)
                x   = sin(Ain) * sin(Azi)
                y   = sin(Ain) * cos(Azi)
                z   = -cos(Ain)
                amp = MT[0]*x*x + MT[1]*y*y + MT[2]*z*z
                amp += 2* (MT[3]*x*y + MT[4]*x*z + MT[5]*y*z)
                if amp > 0:
                    fig[N+j,N+i] = 1
    Y,X = pyl.meshgrid(y,x)
    pyl.pcolor(Y,X,fig)
    pyl.show
if __name__ == "__main__":
    value = sp.cast['g'](10**28);    # dyne
    Mv    = sp.zeros(9, dtype = 'g') # actual source model
    Mv[3] = value
    Mv[6] = x0 * sp.cast['g'](10**5) 
    Mv[7] = y0 * sp.cast['g'](10**5)
    Mv[8] = z0 * sp.cast['g'](10**5)    
    main(Mv)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

############################################################################
#
#                  W phase source inversion package                 
#                               -------------
#
#        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#                      
# (c) California Institute of Technology and Universite de Strasbourg / CNRS 
#                                  April 2013
#
#    Neither the name of the California Institute of Technology (Caltech) 
#    nor the names of its contributors may be used to endorse or promote 
#    products derived from this software without specific prior written 
#    permission
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################

# Append WPHASE BIN to pythonpath
import os,sys
WPHOME = os.path.expandvars('$WPHASE_HOME')
WPBIN = os.path.join(WPHOME,'bin')
sys.path.insert(0,'./')
sys.path.insert(1,WPBIN)

# Avoid writing pyc files
sys.dont_write_bytecode = True

# W PHASE TRACES
from Arguments import *

# Customizing matplotlib
import matplotlib as mpl
mpl.use('PDF')
mpl.rcParams.update(TRACES_PLOTPARAMS)

# Import external modules
import re
import getopt as go
import shutil as sh
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs, feature
from subprocess import call
from copy import deepcopy
from traceback import format_exception

# Import internal modules
import sacpy
import utils
from EQ import EarthQuake

# Internal functions
def showmap(ax,evla,evlo,stla,stlo,coords,flagreg=False):    
    '''
    Show network map
    
    Parameters
    ----------
    ax: matplotlib.axes.Axes
        Axes of the seismogram
    evla: float or ndarray
        Event latitude
    evlo: float or ndarray
        Event longitude
    stla: float
        Station latitudes
    stlo: float
        Station longitudes
    coords: ndarray
        Array of latitudes,longitudes for the entire network
    flagreg: bool
        if False, use a global orthographic map
        if True, use a regional Lambert azimuthal equal-area map
    '''
    
    lonlat = crs.PlateCarree()    
    if flagreg: proj = crs.AzimuthalEquidistant(    central_longitude=evlo, central_latitude=evla)
    else:       proj = crs.cartopy.crs.Orthographic(central_longitude=evlo, central_latitude=evla)

    pos    = ax.get_position().get_points()
    W      = pos[1][0]-pos[0][0] ; H  = pos[1][1]-pos[0][1];
    ax2    = plt.axes([pos[1][0]-W*0.38,pos[0][1]+H*0.01,H*1.08,H],projection=proj)

    if flagreg: ax2.set_extent([evlo-DLON,evlo+DLON,evla-DLAT,evla+DLAT])
    else:       ax2.set_global()
    
    ax2.coastlines(resolution='110m')
    ax2.add_feature(feature.OCEAN, zorder=0, color=(.95,.95,1.))
    ax2.add_feature(feature.LAND,  zorder=0, color=(.8,.8,.8), edgecolor='black')
    ax2.gridlines(color=(.9,.9,.9), zorder=1)
    ax2.scatter(coords[:,1],coords[:,0], color=(1.,.5,0.), marker='o',     s=8,  zorder=2, transform=lonlat)
    ax2.scatter([stlo],[stla],           color=(1.,.0,0.), marker='o',     s=10, zorder=3, transform=lonlat)
    ax2.scatter([evlo],[evla],           color=(0.,.0,1.), marker=(5,1,0), s=20, zorder=4, transform=lonlat)

    # All done
    return


def showPolarmap(ax,az,dist,coords):
    '''
    Show polar map
    
    Parameters
    ----------
    ax: matplotlib.axes.Axes
        Axes of the seismogram
    az: float
        Station azimuth
    dist: float
        Station distance
    coords: ndarray
        Array of latitudes,longitudes for the entire network
    '''

    distlabel  = 6371.0*np.arange(30.0,120.0,30.0)*np.pi/180.0
    pos  = ax.get_position().get_points()
    W  = pos[1][0]-pos[0][0] ; H  = pos[1][1]-pos[0][1] ;        
    ax2 = plt.axes([pos[1][0]-W*0.3,pos[0][1]+H*0.38,H*0.73,H*0.73],polar=True) 
    ax2.plot(coords[:,2],coords[:,3],'yv',ms=4)
    ax2.plot(az,dist,'rv',ms=6)
    ax2.scatter([0.0],[0.0],c='r',marker=(5,1,0),s=70,zorder=1002)
    ax2.set_thetagrids([])
    ax2.set_rgrids(distlabel,labels='')
    # All done
    return;


def disphelp(cmd,solfile,syndir):
    print('Displays W phase traces\n')
    print('usage: %s [option] (for help see %s -h)'%(cmd,cmd))
    print('\nAll parameters are optional:')
    print('   -i, --icmtf          (w)cmt file name (default: %s)'%(solfile))
    print('   -d, --osyndir        output synthetic directory (default: %s)'%(syndir))
    print('   -r, --regional       plot traces for a regional network')
    print('\n   -h, --help         display this help and exit')
    print('\nReport bugs to: <zacharie.duputel@unistra.fr>')
    # All done
    return;


class InvalidOption(Exception):
    """
    Raised if invalid option
    """
    pass


def main(argv):
    # Input parameters (from Arguments.py)
    imaster = IMASTER
    length  = LENGTH_GLOBAL
    syndir  = 'SYNTH_traces'
    o_wpinversion = O_WPINVERSION
    nc = NC
    nl = NL
    solfile = None
    flagreg = False

    # Parse options
    try:
        opts, args = go.gnu_getopt(argv[1:],'i:d:rh',["icmtf=","osydir=","regional","help"])
    except go.GetoptError as err:
        sys.stderr.write('usage: %s [option] (for help see %s -h)\n'%(sys.argv[0],sys.argv[0]))            
        raise

    for o, a in opts:
        if o == '-h' or o == '--help':
            disphelp(sys.argv[0],solfile,syndir)
            sys.exit(0)
        if o == '-r' or o == '--regional':
            length  = LENGTH_REGIONAL
            flagreg = True
        if o == '-i' or o=='--icmtf':
            solfile = a
            if not os.path.exists(solfile):                
                raise IOError('No wcmtfile named %s'%(solfile))
        if o == '-d' or o == '--osyndir':
            syndir = a
    if not solfile:
        for f in ['xy_WCMTSOLUTION','ts_WCMTSOLUTION','WCMTSOLUTION']:
            if os.path.exists(f):
                solfile = f
                break
        if not solfile:
            raise IOError('No wcmtfile available, can be specified with --icmtf')
 
    eq   = EarthQuake()
    eq.rcmtfile(solfile)
    cmtla,cmtlo = eq.lat, eq.lon   

    # Title
    conf  = utils.parseConfig(imaster)
    title = '_'.join(conf['EVNAME'].split())
    title += ',  filter = (%s, %s, %s, %s)'%(conf['filt_cf1'],conf['filt_cf2'],conf['filt_order'],conf['filt_pass']) 

    # Cleanup run dir
    if os.path.exists(syndir) and syndir != '.' and syndir != './':
        utils.rm(syndir)
    if syndir != '.' and syndir != './':
        os.mkdir(syndir)
    if not os.path.exists(LOGDIR):
        os.mkdir(LOGDIR)
    for l in os.listdir('.'):
        if l[:4]=='page' and l[-4:]=='.pdf':
            utils.rm(l)
            
    # Compute synthetics
    cmd    = SYNTHS+' '+imaster+' '+solfile+' '+o_wpinversion+' '+syndir
    print(cmd)
    #status = call(cmd, shell=True, stdin=sys.stdin, stdout=sys.stdout);
    status = os.system(SYNTHS+' '+imaster+' '+solfile+' '+o_wpinversion+' '+syndir+' > '+os.path.join(LOGDIR,'_tmp_synths'))
    if status:        
        print('Error while running '+SYNTHS)
        sys.exit(1)
        
    # Create Sac Objects
    sacdata = sacpy.Sac()
    sacsynt = sacpy.Sac()
    coords = []
    L = open(o_wpinversion).readlines()
    for l in L:
        sacf = l.strip().split()[0]
        sacdata.read(sacf,datflag=0)
        coords.append([sacdata.stla,sacdata.stlo,sacdata.az,sacdata.dist])
    coords = np.array(coords)
    
    # Main loop
    print('Input (W)CMTSOLUTION file is: %s'%(solfile))
    print('Output synthetic directory is: %s'%(syndir))
    perpage = nl*nc
    ntot   = len(L)
    npages = np.ceil(float(ntot)/float(perpage))
    nchan = 1
    count = 1
    pages = 1
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.06,top=0.87,left=0.06,right=0.95,wspace=0.25,hspace=0.35)
    print('All pages will be saved in %s'%(OPDFFILE))
    pp = mpl.backends.backend_pdf.PdfPages(OPDFFILE)
    for l in L:
        # Parse line
        items = l.strip().split()
        fic1 = items[0]
        sacdata.read(fic1)
        chan = sacdata.kcmpnm[0:3]
        loc  = sacdata.khole
        fic2 = syndir+'/%s.%s.%s.%s.complete_synth.bp.sac'\
               %(sacdata.kstnm,sacdata.knetwk,chan,loc)
        sacsynt.read(fic2)        
        # pages
        if count > perpage:
            plt.suptitle(title+ ',   p %d/%d'%(pages,npages), fontsize=16, y=0.95)
            print('page %d/%d'%(pages,npages))
            #fig.set_rasterized(True)
            pp.savefig(orientation='landscape')
            plt.close()
            pages += 1
            count = 1
            fig = plt.figure()
            fig.subplots_adjust(bottom=0.06,top=0.87,left=0.06,right=0.95,wspace=0.25,hspace=0.35)
        # Time - W phase window
        t1 = np.arange(sacdata.npts,dtype='double')*sacdata.delta + sacdata.b - sacdata.o
        t2 = np.arange(sacsynt.npts,dtype='double')*sacsynt.delta + sacsynt.b - sacsynt.o        
        wnb = float(items[5])
        wne = float(items[6])
        wtb = sacdata.b - sacdata.o + wnb * sacdata.delta
        wte = sacdata.b - sacdata.o + wne * sacdata.delta        
        # Plot trace
        ax = plt.subplot(nl,nc,count)
        plt.plot(t1,sacdata.depvar*1000.,'k')
        plt.plot(t2,sacsynt.depvar*1000.,'r-')
        plt.plot([wtb,wte],[0.,0.],'ro')
        # Axes limits
        B=wtb-150.0
        if B<0:
            B = 0.0
        plt.xlim([B,B+length*sacsynt.delta])        
        if YLIM_AUTO:
            a    = np.absolute(sacsynt.depvar[:length]).max()*1000.
            ymin = -1.1*a
            ymax =  1.1*a
            ylims = [ymin,ymax]
        else:
            ylims = YLIMFIXED
        plt.ylim(ylims)        
        # Annotations
        if sacdata.kcmpnm[2] == 'Z':
            label = u'%s %s %s %s (\u03C6,\u0394) = %6.1f\u00B0, %5.1f\u00B0'
            label = label%(sacdata.knetwk,sacdata.kstnm, sacdata.kcmpnm, sacdata.khole,
            sacdata.az, sacdata.gcarc)
        else:
            label = u'%s %s %s %s (\u03C6,\u0394,\u03B1) = %6.1f\u00B0, %5.1f\u00B0, %6.1f\u00B0'
            label = label%(sacdata.knetwk,sacdata.kstnm, sacdata.kcmpnm, sacdata.khole,
            sacdata.az, sacdata.gcarc, sacdata.cmpaz)    
        plt.title(label,fontsize=10.0,va='center',ha='center')
        if not (count-1)%nc:
            plt.ylabel('mm',fontsize=10)
        if (count-1)/nc == nl-1 or nchan+nc > ntot:
            plt.xlabel('time, sec',fontsize=10) 
        plt.grid()
        try:
            showmap(ax,cmtla,cmtlo,sacdata.stla,sacdata.stlo,coords,flagreg)
            pass
        except:
            showPolarmap(ax,sacdata.az,sacdata.dist,coords)
            print('Cannot plot Cartopy map')
            print(format_exception(*sys.exc_info()))
        count += 1
        nchan += 1
    print('page %d/%d'%(pages,npages))
    #fig.set_rasterized(True)
    plt.suptitle(title + ',    p %d/%d'%(pages,npages), fontsize=16, y=0.95)
    pp.savefig(orientation='landscape')
    plt.close()
    pp.close()


if __name__=='__main__':
    main(sys.argv)

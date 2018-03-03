#!/usr/bin/env python
# *-* coding: iso-8859-1 *-*

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

# GRID SEARCH PLOTS

# Customizing matplotlib
import matplotlib as mpl
mpl.use('AGG')

# External modules
import sys
import getopt as go
import numpy as np
import matplotlib.pyplot as plt

# Avoid writing pyc files
sys.dont_write_bytecode = True

def wraplons(lons):
    for i in range(len(lons)):
        if lons[i]<0.:
            lons[i]+=360.
    # All done
    return;

def interp (j, n, a, b):
    if n == 1:
        v = a
    else:
        v = ((n-1-j)*a + j*b)/(n-1)
    # All done
    return v

def findcdep(depth,depths):
    N = len(depths)
    for i,dep in zip(range(N),depths):
        if int(depth*100) == int(dep*100):
            return i
    # All done
    return -1

def rxyzgfile(ifile):
    fid= open(ifile,'r')
    L=fid.readlines()
    fid.close()
    latopt, lonopt, depopt, rmsopt = map(float,L[0].strip('\n').split())    
    latpde, lonpde, deppde, rmspde = map(float,L[1].strip('\n').split())
    lons = [] ; lats = [] ; deps = [] ; rms  = [] ; 
    rmsdepths = [] ; depths = []    
    for l in L[2:]:
        items = l.strip().split()
        dep = float(items[6])
        err = float(items[7])
        i   = findcdep(dep,depths)
        if i>=0:
            if err < rmsdepths[i]:
                rmsdepths[i] = err
        else:
            depths.append(dep)
            rmsdepths.append(err)
        lats.append(float(items[4]))
        lons.append(float(items[5]))
        deps.append(dep)
        rms.append(err)
    lats=np.array(lats); lons=np.array(lons); deps=np.array(deps)
    rms=np.array(rms)    
    depths = np.array(depths); rmsdepths = np.array(rmsdepths)    
    idepths   = np.argsort(depths)
    depths    = depths[idepths]
    rmsdepths = rmsdepths[idepths]
    # All done
    return [latopt,lonopt,depopt,rmsopt,latpde,lonpde,deppde,rmspde,lats,lons,deps,rms,depths,rmsdepths]

def r_xy_gfile(ifile):
    fid= open(ifile,'r')
    L=fid.readlines()
    fid.close()
    latopt, lonopt, depopt, rmsopt = map(float,L[0].strip('\n').split())    
    latpde, lonpde, deppde, rmspde = map(float,L[1].strip('\n').split())
    lat = []; lon = []; rms = []
    for l in L[2:]:
            tmp = l.strip('\n').split()
            lat.append(float(tmp[4]))
            lon.append(float(tmp[5]))
            rms.append(float(tmp[7]))
    lat = np.array(lat)
    lon = np.array(lon)
    rms = np.array(rms)            
    # All done
    return [latopt,lonopt,depopt,rmsopt,latpde,lonpde,rmspde,lat,lon,rms]

def plotxyz(ifilexyz='grid_search_xyz_out',ifilexy='grid_search_xy_out',flag=False,mksmin=1.,mksmax=300.):
    if not flag:
        try:            
            from mpl_toolkits.mplot3d import Axes3D
        except:
            flag = True
    latopt,lonopt,depopt,rmsopt,latpde,lonpde,deppde,rmspde,lats,lons,deps,rmss,depths,rmsdep = rxyzgfile(ifilexyz) 
    fig = plt.figure(figsize=(7.6875, 6.125))
    if flag:
        ofile='grid_search_z.pdf'
        minrms = rmsopt
        maxrms = rmss.max()
        nrmsdep = rmsdep/minrms*100.0
        nrmsopt = 100.0
        rmsmax  = np.ceil((max(nrmsdep)+1.0)/10.)*10.
        plt.plot(depths,nrmsdep,'bo-',ms=4)
        plt.plot([deppde],[rmspde/minrms*100.],'k+',ms=14,mew=2.5,alpha=0.7)
        plt.plot([depopt],[nrmsopt],'rv',ms=14,alpha=0.7)
        plt.xlabel('Centroid depth, km')
        plt.ylabel('Normalized RMS misfit')
        plt.ylim([99.,rmsmax])
        plt.grid()
    else:
        ofile='grid_search_xyz.pdf'
        latopt,lonopt,depopt,rmsopt,latpde,lonpde,rmspde,lat,lon,rms = r_xy_gfile(ifilexy)
        minrms  = rmsopt
        maxrms  = rmss.max()
        maxrms2 = rms.max()
        if maxrms2>maxrms:
            maxrms=maxrms2    
        nrmss   = rmss/rmsopt*100.
        nrms    = rms/rmsopt*100.        
        minrms  = 100.
        maxrms  = maxrms/rmsopt * 100.        
        ax2 = plt.axes([0.825,0.2,0.1,0.6])
        cm = plt.get_cmap('jet')
        print(minrms,maxrms)
        for i in range(8):
            Bpos   = interp(i,8,0.,1.)
            Brms=interp(i,8,minrms,maxrms)
            mksize = ((Brms-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
            plt.scatter(0.5,Bpos,c=cm((Brms-minrms)/(maxrms-minrms)),marker='o',s=mksize)
            plt.text(0.9,Bpos-0.02,'%7.1f'%Brms)
        ax2.set_axis_off()
        mksize = ((np.array(nrmss)-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
        plt.ylim(-0.1,1.1)
        plt.xlim(0,1.6)
        plt.title('Normalized RMS')
        ax1 = plt.axes([0.05,0.1,0.7,0.8], projection='3d',azim=-53.,elev=13.)

        cstla = 6371.*np.pi/180.
        cstlo = cstla*np.cos(latpde*np.pi/180.)
        for la,lo,z,err,siz in zip(lats,lons,deps,nrmss,mksize):
            x = (lo - lonpde)*cstlo
            y = (la - latpde)*cstla
            col = cm((err-minrms)/(maxrms-minrms))            
            ax1.scatter(x,y,z,c=col,marker='o',s=siz)
        mksize = ((np.array(nrms)-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
        for la,lo,err,siz in zip(lat,lon,nrms,mksize):
             x = (lo - lonpde)*cstlo
             y = (la - latpde)*cstla
             col = cm((err-minrms)/(maxrms-minrms))        
             ax1.scatter(x,y,depopt,c=col,marker='o',s=siz)
        x = (lonopt-lonpde)*cstlo
        y = (latopt-latpde)*cstla
        ax1.scatter(0,0,deppde,c='k',marker='+',s=140)
        ax1.scatter(x,y,depopt,c='r',marker='v',s=140)
        ax1.set_xlabel('East, km')
        ax1.set_ylabel('North, km')
        ax1.set_zlabel('Depth, km')
    plt.savefig(ofile)
    # All done
    return;

def prepColorbar(minz,maxz,Lref):
    if -minz>=maxz:
        Laxo = Lref
        Laxc = -Lref/minz * maxz
        tico = np.linspace(minz,0.,4)
        ticc = np.arange(0.,maxz,tico[1]-tico[0])
        if len(ticc)<=1:
            if Laxc >= 0.04:
                ticc = np.array([maxz])
            else:
                Laxc = 0.
        else:
            ticc = list(ticc[1:])
            ticc.reverse()
            ticc = np.array(ticc)        
    else:
        Laxc = Lref
        Laxo = -Lref/maxz * minz
        ticc = np.linspace(0.,maxz,4)
        tico = -np.arange(0.,-minz,ticc[1]-ticc[0])
        if len(tico)<=1:
            if Laxo >= 0.04:
                tico = np.array([minz])
            else:
                Laxo = 0.
        else:
            tico = list(tico[1:])
            tico.reverse()
            tico = np.array(tico)
    # All done
    return Laxo,Laxc,tico,ticc

def concatCmap(cmaps,offs,cuts,prop): 
    nps = []
    idx = []
    Nc  = len(cmaps)    
    cmprop = 0.
    for k in range(Nc):
        nps.append(len(cmaps[k]._segmentdata['blue']) - offs[k] - cuts[k])
        tmp = np.linspace(cmprop,cmprop+prop[k],nps[-1])
        idx.extend(list(tmp))
        cmprop += prop[k]        
    cdict = {} # RGB dictionary
    for key in ['red','green','blue']:
        cdict[key] = []
    anp   = 0
    for k in range(Nc): 
        for key in ['red','green','blue']:
            for i in range(nps[k]):
                cur = cmaps[k]._segmentdata[key][i+offs[k]]
                cdict[key].append((float(int(idx[i+anp]*1000))/1000.,cur[1],cur[2]))
        anp += nps[k]
    colmap = matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)
    # All done
    return colmap

def plot_etopo(file,m,ax):
    from copy import deepcopy
    try:
        from mpl_toolkits.basemap import NetCDFFile
    except:
        from netCDF4 import Dataset as NetCDFFile
    latll = m.llcrnrlat
    latur = m.urcrnrlat
    lonll = m.llcrnrlon
    lonur = m.urcrnrlon    
    # Read NetCDF ETOPO file
    try:
        data=NetCDFFile(file)
    except:
        print('plot_etopo: Incorrect ETOPO NetCDF file')
        raise 'WARNING: error encountered while plotting ETOPO'
    lats = data.variables['lat'][:]
    lons = deepcopy(data.variables['lon'][:])
    wraplons(lons)
    ila = [] ; ilo = []
    for i in range(len(lats)):
        if lats[i] >= latll and lats[i] <=latur:
            ila.append(i)
    for i in range(len(lons)):
        if lons[i] >= lonll and lons[i] <=lonur:
            ilo.append(i)
    ila = np.array(ila) ; ilo = np.array(ilo)
    inds = np.argsort(lons[ilo]) ; ilo = ilo[inds] ;
    la  = lats[ila] ; lo  = lons[ilo] ;
    z = data.variables['z'][ila[0]:ila[-1]+1,ilo]
    lon, lat = np.meshgrid(lo, la)
    x, y = m(lon, lat)
    # Colormaps
    Lref = 0.35 ; # Colorbar half length
    minz = z.min() ; maxz = z.max() ;
    Laxo,Laxc,ticko,tickc= prepColorbar(minz,maxz,Lref)
    H = float(Laxo+Laxc)/2.0
    oceancmap = plt.cm.Blues_r # Ocean depth colormap
    if Laxc > 0.:              # Elevation colormap
        cmaps  = [plt.cm.YlGn,plt.cm.BrBG]
        offs   = [0  ,  0]
        cuts   = [0  ,  5]
        prop   = [0.5,0.5]
        elevcmap = concatCmap(cmaps,offs,cuts,prop)
    else:
        elevcmap = plt.cm.YlGn    
    # Ocean contour
    if minz < 0.:
        plt.axes(ax)
        zo = np.ma.masked_where(z >= 0.,z)
        co = m.contourf(x,y,zo,30,cmap=oceancmap)
        if Laxo>0.:   # Ocean depth colorbar
            caxo = plt.axes([0.45-H,0.04,Laxo,0.02])        
            plt.title('Ocean Depth, m',fontsize='medium')
            cbo = plt.colorbar(mappable=co,cax=caxo,ticks=ticko,format='%.0f',
                       orientation='horizontal')
    # Land contour
    if maxz >= 0.:
        plt.axes(ax)
        zc = np.ma.masked_where(z <  0.,z)
        cc = m.contourf(x,y,zc,30,cmap=elevcmap)
        if Laxc > 0.: # Elevation colorbar
            caxc = plt.axes([0.45-H+Laxo,0.04,Laxc,0.02])
            plt.title('Elevation, m',fontsize='medium')
            cbc = plt.colorbar(mappable=cc,cax=caxc,ticks=tickc,format='%.0f',
                       orientation='horizontal')
        plt.axes(ax)

def plotxy(ifile='grid_search_xy_out',ofile='grid_search_xy.pdf',basemapflag=False,mksmin=1.,
        mksmax=30.,delta=1.0,resolution = 'h'):
    # Initialize variables
    rms  = []
    lon  = []
    lat  = []
    flag = False
    # Read file
    latopt,lonopt,depopt,rmsopt,latpde,lonpde,rmspde,lat,lon,rms = r_xy_gfile(ifile)    
    # RMS Scale
    minrms = rmsopt    
    nrms = rms/minrms*100.
    maxrms = max(rms)/minrms*100.
    minrms = 100.
    plt.figure(figsize=(9.6125,  8.1))
    ax1 = plt.axes([0.04,0.13,0.8,0.85])
    ax2 = plt.axes([0.85,0.2,0.1,0.6])
    cm = plt.get_cmap('jet')
    for i in range(8):
        Bpos   = interp(i,8,0.,1.)
        Brms=interp(i,8,minrms,maxrms)
        mksize = ((Brms-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
        plt.plot([0.5],[Bpos],'ko',ms=mksize,markerfacecolor=cm((Brms-minrms)/(maxrms-minrms)))
        plt.text(0.9,Bpos-0.02,'%7.1f'%Brms)
    ax2.set_axis_off()
    mksize = ((np.array(nrms)-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin    
    plt.ylim(-0.1,1.1)
    plt.xlim(0,1.6)
    plt.title('Normalized RMS')
    # DISPLAY MAP
    plt.axes(ax1)        
    if basemapflag:
        try:
            from mpl_toolkits.basemap import Basemap
        except:
            print('WARNING: No module named basemap')
            print('   The mpl_toolkits.basemap module is necessary')
            print('   if you want to plot bathymetry and coastlines')
            basemapflag = False
    if basemapflag:
        from os.path import expandvars,exists        
        wraplons(lon)
        deltalon = delta/np.cos(np.pi*latpde/180.0)
        latll = lat.min() - delta ; latur = lat.max() + delta ;
        lonll = lon.min() - deltalon ; lonur = lon.max() + deltalon ;        
        # Basemap instance
        m = Basemap(projection='merc',llcrnrlon=lonll,llcrnrlat=latll,urcrnrlon=lonur,
                urcrnrlat=latur,resolution=resolution)
        # Bathymetry        
        ETOPO_file = expandvars('$ETOPOFILE')
        if exists(ETOPO_file):
            try:
                plot_etopo(ETOPO_file,m,ax1)
            except:
                print('WARNING: error encountered while plotting ETOPO')
                print('         Will not display topography/bathymetry')
        else:
            if ETOPO_file[0]=='$':
                print('WARNING: Undefined environment variable $ETOPOFILE')
            else:
                print('WARNING: ETOPOFILE=%s does not exists'%(ETOPO_file))
                print('         Will not display topography/bathymetry')
        # Coastlines/meridians/paralells
        plt.axes(ax1)
        m.drawcoastlines(linewidth=0.3)
        m.drawmeridians(np.arange(float(int(lonll)),lonur+delta,delta),labels=[0,0,0,1],
                dashes=[1,1],linewidth=0.5,color='k')
        m.drawparallels(np.arange(float(int(latll)),latur+delta,delta),labels=[1,0,0,0],
                dashes=[1,1],linewidth=0.5,color='k')
        # RMS misfit
        for la,lo,err,siz in zip(lat,lon,nrms,mksize):
            x,y = m(lo,la)
            col = cm((err-minrms)/(maxrms-minrms))        
            m.plot([x],[y],c=col,marker='o',ms=siz)
        l = [lonpde,lonopt]
        wraplons(l)
        xpde,ypde = m(l[0],latpde)
        xopt,yopt = m(l[1],latopt)
        m.plot([xpde],[ypde],'k+',ms=14,mew=2.5,alpha=0.7)
        m.plot([xopt],[yopt],'rv',ms=14,alpha=0.7)
    else:
        deltalon = delta/np.cos(np.pi*latpde/180.0)
        latll = lat.min() - delta ; latur = lat.max() + delta ;
        lonll = lon.min() - deltalon ; lonur = lon.max() + deltalon ;
        for la,lo,err,siz in zip(lat,lon,nrms,mksize):
            col = cm((err-minrms)/(maxrms-minrms))        
            plt.plot([lo],[la],c=col,marker='o',ms=siz)
        plt.plot([lonpde],[latpde],'k+',ms=14,mew=2.5,alpha=0.7)
        plt.plot([lonopt],[latopt],'rv',ms=14,alpha=0.7)
        plt.xlim([lonll,lonur])
        plt.ylim([latll,latur])
    plt.savefig(ofile)
    # All done
    return;
    

def plotts(ifile='grid_search_ts_out',ofile='grid_search_ts.pdf'):
    fid= open(ifile,'r')
    L=fid.readlines()
    fid.close()
    tsopt,rmsopt = map(float,L[0].strip('\n').split())
    tsini,rmsini = map(float,L[1].strip('\n').split())    
    ts  = []
    rms = []
    for l in L[2:]:
        tmp = l.strip('\n').split()
        ts.append(float(tmp[2]))
        rms.append(float(tmp[7]))
    rms = np.array(rms)
    nrms = rms/rmsopt*100.
    rmsini = rmsini/rmsopt*100.
    rmsopt = 100.
    rmsmax  = np.ceil((max(nrms)+1.0)/10.)*10.
    plt.plot(ts,nrms,'bo',ms=4)
    plt.plot([tsini],[rmsini],'k+',ms=14,mew=2.5,alpha=0.7)
    plt.plot([tsopt],[rmsopt],'rv',ms=14,alpha=0.7)
    plt.grid('on')
    plt.ylim([90.,rmsmax])
    plt.xlabel('Centroid time shift, ts (sec.)')
    plt.ylabel('Normalized RMS')
    plt.savefig(ofile)
    # All done
    return;

def usage(cmd):
    print('usage: %s [option] (for help see %s -h)'%(cmd,cmd))
    # All done
    return;


def disphelp(cmd):
    print('Displays grid search results\n')
    usage(cmd)
    print('\nAll parameters are optional:')
    print('   -t, --onlyts         centroid time-shift grid search (ts) only')
    print('   -p, --onlyxy         centroid position grid search (xy) only')
    print('   -b, --basemap        display coastlines and bathymetry')
    print('   --its \'file\'       set input ASCII file for ts (grid_search_ts_out)')
    print('   --ixy \'file\'       set input ASCII file for xy ((grid_search_xy_out))')
    print('   --ots \'file\'       set output png file for ts (grid_search_ts.pdf)')
    print('   --oxy \'file\'       set output png file for xy ((grid_search_xy.pdf))')
    print('\n   -h, --help           display this help and exit')
    print('\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>')
    # All done
    return;
        

def main(argv):
    try:
        opts, args = go.gnu_getopt(argv[1:],'tphzb',["onlyts","onlyxy","basemap","its=","ixy=","ots=","oxy=","help"])
    except go.GetoptError as err:
        usage(sys.argv[0])
        raise
    flagts  = True
    flagxy  = True
    flagxyz = False
    basemap = False
    ts_ifile='grid_search_ts_out'
    ts_ofile='grid_search_ts.pdf'
    xy_ifile='grid_search_xy_out'
    xy_ofile='grid_search_xy.pdf'
    for o, a in opts:
        if o == '-h' or o == '--help':
            disphelp(sys.argv[0])
            sys.exit(0)
        if o == '-t' or o == '--onlyts':
            if not flagts:
                usage(sys.argv[0])
                raise go.GetoptError('Options -t and -p cannot be used simultaneously')
            flagxy = False
            flagts = True
        if o == '-p' or o == '--onlyxy':
            if not flagxy:
                usage(sys.argv[0])
                raise go.GetoptError('Options -t and -p cannot be used simultaneously')                
            flagts = False
            flagxy = True
        if o == '--its':
            ts_ifile = a
        if o == '--ixy':
            xy_ifile = a
        if o == '--ots':
            ts_ofile = a
        if o == '--oxy':
            xy_ofile = a
        if o == '-b' or o=='--basemap':
            basemap = 1
        if o == '-z':
            flagxyz = True

    if flagts:
        plotts(ts_ifile,ts_ofile)
    if flagxy:
        plotxy(xy_ifile,xy_ofile,basemapflag=basemap)
    if flagxyz:
        #plotxyz(flag=False)
        plotxyz(flag=True)


if __name__=='__main__':
    main(sys.argv)

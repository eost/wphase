#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

######################################
# GRID SEARCH PLOTS
###
# Z.Duputel, L.Rivera and H.Kanamori
# 2009/09/08 -- initial version

import sys
import getopt as go
import pylab as pyl



def interp (j, n, a, b):
	if n == 1:
		v = a
	else:
		v = ((n-1-j)*a + j*b)/(n-1)
	return v

def plot_xy(ifile='grid_search_xy_out',ofile='grid_search_xy.png',mksmin=1.,mksmax=30.):
	# Initialize variables
	rms  = []
	lon  = []
	lat  = []
	flag = 0
	# Read file
	fid= open(ifile,'r')
	L=fid.readlines()
	fid.close()
	
	# get lat,lon,rms
	latopt, lonopt, rmsopt = map(float,L[0].strip('\n').split())	
	latpde, lonpde, rmspde = map(float,L[1].strip('\n').split())
	for l in L[2:]:
            tmp = l.strip('\n').split()
            lat.append(float(tmp[4]))
            lon.append(float(tmp[5]))
            rms.append(float(tmp[7]))

	# RMS Scale 
        minrms = rmsopt
	maxrms = max(rms)        
	pyl.figure(figsize=(7.6875, 6.125))
	ax1 = pyl.axes([0.1,0.1,0.7,0.8])
	ax2 = pyl.axes([0.85,0.2,0.1,0.6])
	for i in xrange(8):
		Bpos   = interp(i,8,0.,1.)
		Brms=interp(i,8,minrms,maxrms)
		mksize = ((Brms-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
		pyl.plot([0.5],[Bpos],'ko',ms=mksize)
		pyl.text(0.9,Bpos-0.02,'%7.4f'%Brms)
	ax2.set_axis_off()
	mksize = ((pyl.array(rms)-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
	pyl.ylim(-0.1,1.1)
	pyl.xlim(0,1.6)
	pyl.title('RMS (mm)')
	pyl.axes(ax1)
	for i in xrange(len(L)-2):
            pyl.plot([lon[i]],[lat[i]],'ko',ms=mksize[i])
	pyl.plot([lonpde],[latpde],'rv',ms=6)
       	pyl.plot([lonopt],[latopt],'bv',ms=6)
	pyl.axis('equal')
	pyl.savefig(ofile)
	pyl.close()
	

def plot_ts(ifile='grid_search_ts_out',ofile='grid_search_ts.png'):
	fid= open(ifile,'r')
	L=fid.readlines()
	fid.close()
	tsopt,rmsopt = map(float,L[0].strip('\n').split())
	tsini,rmsini = map(float,L[1].strip('\n').split())	
	ts  = []
	rms = []
	for l in L[2:]:
		tmp = l.strip('\n').split()
		ts.append(float(tmp[1]))
		rms.append(float(tmp[5]))
	pyl.plot(ts,rms,'bo',ms=4)
	pyl.plot([tsini],[rmsini],'rv')
	pyl.plot([tsopt],[rmsopt],'bv')
	pyl.grid('on')
	pyl.xlabel('Centroid time shift, ts (sec.)')
	pyl.ylabel('RMS (mm)')
	pyl.savefig(ofile)
	pyl.close()

def usage():
	print 'usage: make_grids [-f] [-t] [-p] [-i] ... [--help]'


def disphelp():
	print 'Display grid search results\n'
	usage()
	print '\nAll parameters are optional:'
	print '   -t, --onlyts         centroid time-shift grid search (ts) only'
	print '   -p, --onlyxy         centroid position grid search (xy) only'
	print '   --its \'file\'       set input ASCII file for ts (grid_search_ts_out)'
	print '   --ixy \'file\'       set input ASCII file for xy ((grid_search_xy_out))'
        print '   --ots \'file\'       set output png file for ts (grid_search_ts.png)'
	print '   --oxy \'file\'       set output png file for xy ((grid_search_xy.png))'
	print '\n   -h, --help           display this help and exit'
	print '\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>'
        
##### MAIN #####	
if __name__ == "__main__":
    try:
        opts, args = go.gnu_getopt(sys.argv[1:],'tph',["onlyts","onlyxy","its=","ixy=","ots=","oxy=","help"])
    except go.GetoptError, err:
        print '*** ERROR ***'
        print str(err)
        usage()
        sys.exit(1)

    flagts = 1
    flagxy = 1
    ts_ifile='grid_search_ts_out'
    ts_ofile='grid_search_ts.png'
    xy_ifile='grid_search_xy_out'
    xy_ofile='grid_search_xy.png'

    for o, a in opts:
        if o == '-h' or o == '--help':
            disphelp()
            sys.exit(0)
        if o == '-t' or o == '--onlyts':
            if flagts == 0:
                print '** ERROR (options -t and -p cannot be used simultaneously) **'
                usage()
                sys.exit(1)
            flagxy = 0
            flagts = 1
        if o == '-p' or o == '--onlyxy':
            if flagxy == 0:
                print '** ERROR (options -t and -p cannot be used simultaneously) **'
                usage()
                sys.exit(1)
            flagts = 0
            flagxy = 1
        if o == '--its':
            ts_ifile = a
        if o == '--ixy':
            xy_ifile = a
        if o == '--ots':
            ts_ofile = a
        if o == '--oxy':
            xy_ofile = a

    if flagts:
        plot_ts(ts_ifile,ts_ofile)
    if flagxy:
        plot_xy(xy_ifile,xy_ofile)

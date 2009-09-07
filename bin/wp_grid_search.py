#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

######################################
# GRID SEARCH FOR WPHASE INVERSION
###
# Z.Duputel, L.Rivera and H.Kanamori
# 2009/07/15 -- initial version for grid search
# 2009/07/19 -- optimization for time shift : B-tree sampling method
# 2009/07/26 -- optimization for centroïd position : Oct-tree sampling method

import os,re,shutil,sys,time,calendar,getopt
from EQ import *

import pylab as pyl

WPHOME = os.path.expandvars('$WPHASE_HOME')
if WPHOME[-1] != '/':
	WPHOME += '/'

BIN = WPHOME+'bin/'

REPREPARE_TS = BIN+'reprepare_wp_ts.csh'
WPINV_TS     = BIN+'wpinversion -imas ts_i_master -ofil ts_o_wpinversion -ocmtf ts_WCMTSOLUTION -ps ts_p_wpinversion '+\
	       '-wpbm ts_wpinv.pgm -log LOG/_ts_wpinversion.log -osyndir ts_SYNTH -pdata ts_predicted_data -nt '

CALCSYN      = BIN+'recalc_fast_synths.csh'
RECALCSYN_XY = BIN+'recalc_fast_synths.csh'
REPREPARE_XY = BIN+'reprepare_wp_xy.csh'
WPINV_XY     = BIN+'wpinversion -imas xy_i_master -ofil xy_o_wpinversion -ocmtf xy_WCMTSOLUTION -ps xy_p_wpinversion '+\
	       '-wpbm xy_wpinv.pgm -log LOG/_xy_wpinversion.log -osyndir xy_SYNTH -pdata xy_predicted_data -nt'


def plot_xy(file,mksmin=1.,mksmax=30.,minrms=-99.,maxrms=-99.):
	# Initialize variables
	rms  = []
	lon  = []
	lat  = []
	flag = 0
	# Read file
	fid= open(file,'r')
	L=fid.readlines()
	fid.close()
	# get lat,lon,rms
	tmp  = L[0].strip('\n').split()
	lat0 = float(tmp[0])
	lon0 = float(tmp[1])
	if (minrms < 0.):
		minrms = float(tmp[2])
		flag   = 1
	for l in L[1:]:
		tmp = l.strip('\n').split()
		lat.append(float(tmp[4]))
		lon.append(float(tmp[5]))
		rms.append(float(tmp[7]))
	if (maxrms < 0):
		maxrms = max(rms)
		flag   = 1
	mksize = ((pyl.array(rms)-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
	for i in xrange(len(L)-1):
		pyl.plot([lon[i]],[lat[i]],'ko',ms=mksize[i])
	if flag:
		return [minrms,maxrms]

def plot_ts(file,label):
	fid= open(file,'r')
	L=fid.readlines()
	fid.close()
	ts  = []
	rms = []
	for l in L[1:]:
		tmp = l.strip('\n').split()
		ts.append(float(tmp[1]))
		rms.append(float(tmp[5]))
	pyl.plot(ts,rms,label)
	pyl.xlabel('Centroid time shift, ts (sec.)')
	pyl.ylabel('RMS (mm)')

def nprocessors():
	try:
		s   = open('/proc/cpuinfo', 'r').read()
		npr = s.replace(' ', '').replace('\t', '').count('processor:')
		if npr < 1:
			return 1
		else:
			return npr
	except:
		return 1
	
def interp (j, n, a, b):
	if n == 1:
		v = a
	else:
		v = ((n-1-j)*a + j*b)/(n-1)
	return v

def grep(chaine, file):
	out = [];
	rms = re.compile(chaine)
	ps  = open(file, 'r')
	for line in ps:
		if rms.match(line):
			out.append(line)
	ps.close()
	return(out)

def grep2(list, file):
	out   = [];
	ps    = open(file, 'r')
	lines = ps.readlines()
	ps.close()
	for chaine in list:
		rexp = re.compile(chaine)
		for line in lines:
			if rexp.match(line):
				out.append(line)
				break
	return(out)

def addrefsol(cmtref,cmtfile):
	cmtf = open(cmtref,'r')
	L=cmtf.readlines()
	cmtf.close()
	cmtf = open(cmtfile,'a')
	if len(L) < 13:
		print '*** ERROR (reading reference solution) ***' 
		print 'incomplete cmtfile: %s'%(cmtref)
		sys.exit(1)
	for l in L[7:]:
		cmtf.write(l)
	cmtf.close()

def add_coor(coor,lat,lon):
	for cds in coor:
		if int(lat*100) == int(cds[0]*100) and int(lon*100) == int(cds[1]*100):
			return coor
	coor.append([lat,lon])
	return coor

def grid_search_xy(datdir,dmin,dmax,cmtref,ftable,eq,ts,hd,wpwin=[15.],flagref=0,out='stdout'):
	if out == 'stdout':
		fid = sys.stdout
	else:
		fid = open(out,'w')
	fid.write('CENTROID POSITION GRID SEARCH\n')

	# Initialize variables
 	eq_gs = EarthQuake()
 	EQcopy(eq_gs,eq)
	
	Nit  = 2
	dx   = 0.8
	lat1 = eq.lat - 2.0
	lat2 = eq.lat + 2.0
	lon1 = eq.lon - 2.0
	lon2 = eq.lon + 2.0

	n    = (lat2-lat1)/dx+1
	lats = pyl.arange(lat1,lat2+dx,dx)
	lons = pyl.arange(lon1,lon2+dx,dx)
	lons = lons[:n]
	lats = lats[:n]
	coor = []
	for lat in lats:
		for lon in lons:
			coor.append([lat,lon])
	
	format = '%03d %03d %8.2f %8.2f %8.2f %8.2f %8.2f %12.7f %12.7f\n'
	
	cmttmp = cmtref+'_xy_tmp'
	eq.wimaster(datdir,dmin,dmax,ftable,cmttmp,'xy_i_master','./xy_GF/',wpwin) 	
	if os.access('xy_SYNTH',os.F_OK):
		shutil.rmtree('xy_SYNTH')
	if os.access('xy_DATA',os.F_OK):
		shutil.rmtree('xy_DATA')
	os.mkdir('xy_SYNTH')
	os.mkdir('xy_DATA')

	Nopt   = [4,4,2,1]
	optrms = pyl.ones(Nopt[0],dtype='float64')*1.e10 #### a voir
	optlat = pyl.ones(Nopt[0],dtype='float64')*lat2
	optlon = pyl.ones(Nopt[0],dtype='float64')*lon1
	for it in xrange(Nit):
		if it != 0:
			dx = dx/2.
			coor = []
			for i in xrange(Nopt[it]):
				add_coor(coor,optlat[i]+dx,optlon[i]-dx)
				add_coor(coor,optlat[i]+dx,optlon[i])
				add_coor(coor,optlat[i]+dx,optlon[i]+dx)
				add_coor(coor,optlat[i]   ,optlon[i]+dx)
				add_coor(coor,optlat[i]-dx,optlon[i]+dx)
				add_coor(coor,optlat[i]-dx,optlon[i])
				add_coor(coor,optlat[i]-dx,optlon[i]-dx)
				add_coor(coor,optlat[i]   ,optlon[i]-dx)

		fid.write('Iteration %d:\n' % (it+1))
		tmp_table  = open('_tmp_xy_table', 'w')
		ncel = 0
		for cds in coor:
			eq_gs.lat = cds[0]
			eq_gs.lon = cds[1]
			eq_gs.wcmtfile(cmttmp,ts,hd)
			os.system(RECALCSYN_XY)
			os.system(REPREPARE_XY)
			os.system(WPINV_XY)
			out  = grep(r'^W_cmt_err:', 'LOG/_xy_wpinversion.log')
			rms  = float(out[0].strip('\n').split()[1])
			nrms = float(out[0].strip('\n').split()[2])
			fid.write('   cell %3d : lat=%8.3fdeg lon=%8.3fdeg, rms = %12.7f mm\n'% (ncel+1,eq_gs.lat,eq_gs.lon,rms))
			for i in xrange(Nopt[it]):
				if rms < optrms[i]:
					for j in xrange(Nopt[it]-1,i-1,-1):
						optrms[j] = optrms[j-1]
						optlat[j] = optlat[j-1]
						optlon[j] = optlon[j-1]
					optrms[i] = rms
					optlat[i] = eq_gs.lat
					optlon[i] = eq_gs.lon
					break
			ncel += 1
			tmp_table.write(format%(-99,-99,ts,hd,eq_gs.lat,eq_gs.lon,eq_gs.dep,rms,nrms))
			tmp_table.flush()
		tmp_table.close()
		fid.write('Optimum centroid location: %8.3f %8.3f;  rms = %12.7f mm\n'%(optlat[0], optlon[0], optrms[0]))
		tmp_table = open('_tmp_xy_table', 'r')
		out_table = open('grid_search_xy_out%d'%(it+1), 'w')
		out_table.write('%8.3f %8.3f %12.7f\n'%(optlat[0], optlon[0], optrms[0]))
		out_table.write(tmp_table.read())
		out_table.close()
		tmp_table.close()
		os.remove('_tmp_xy_table')

	eq_gs.lat = optlat[0]
	eq_gs.lon = optlon[0]
	eq_gs.wcmtfile(cmttmp,ts,hd)
	os.system(RECALCSYN_XY+'> /dev/null')
	os.system(REPREPARE_XY+'> /dev/null')
	if flagref:
		addrefsol(cmtpde,cmttmp)
		os.system(WPINV_XY+' -ref')
	else:
		os.system(WPINV_XY)

	######### USE MATPLOT LIB ######### 
	mksmin = 1
	mksmax = 30
	pyl.figure(figsize=(7.6875, 6.125))
	ax1 = pyl.axes([0.1,0.1,0.7,0.8])
	[minrms,maxrms]=plot_xy('grid_search_xy_out1',mksmin,mksmax)
	ax2 = pyl.axes([0.85,0.2,0.1,0.6])
	for i in xrange(8):
		Bpos   = interp(i,8,0.,1.)
		Brms=interp(i,8,minrms,maxrms)
		mksize = ((Brms-minrms)/(maxrms-minrms))*(mksmax-mksmin)+mksmin
		pyl.plot([0.5],[Bpos],'ko',ms=mksize)
		pyl.text(0.9,Bpos-0.02,'%7.4f'%Brms)
	ax2.set_axis_off()
	pyl.ylim(-0.1,1.1)
	pyl.xlim(0,1.6)
	pyl.title('RMS (mm)')
	pyl.axes(ax1)
	for j in xrange(1,Nit):
		plot_xy('grid_search_xy_out%d'%(it+1),mksmin,mksmax,minrms,maxrms)
	pyl.plot([eq.lon],[eq.lat],'rv',ms=6)
       	pyl.plot([optlon[0]],[optlat[0]],'bv',ms=6)
	pyl.axis('equal')
	pyl.savefig('grid_search_xy.png')
	pyl.close()		



def grid_search_ts(datdir,dmin,dmax,cmtref,ftable,eq,tsini,hdini,wpwin=[15.],flagref=0,out='stdout'):
	if out == 'stdout':
		fid = sys.stdout
	else:
		fid = open(out,'w')

	out     = grep(r'^W_cmt_err:', 'LOG/wpinversion.log')			
	rmsini  = float(out[0].strip('\n').split()[1])
	nrmsini = float(out[0].strip('\n').split()[2])
	format     = '%02d %8.2f %8.2f %8.2f %8.2f %12.8f %12.8f\n'

	fid.write('CENTROID TIME DELAY GRID SEARCH\n')
	Nit = 3
	Sts = [4.,4.,2.,1.]
	
	# if eq.mag <= 7.0:
	# 	ts1 = 4.
	# 	ts2 = 20.
	# elif eq.mag <8.0:
	# 	ts1 = 8.
	# 	ts2 = 48.
	# else:
	# 	ts1 = 14. 
	# 	ts2 = 56.
	
	ts1 = 1.
	ts2 = tsini*3.
	if ts2 > 60.:
		ts2 = 60.

	lat = eq.lat
	lon = eq.lon
	dep = eq.dep
	
	cmttmp = cmtref+'_ts_tmp'
	eq.wimaster(datdir,dmin,dmax,ftable,cmttmp,'ts_i_master','ts_GF',wpwin=[15.])
	
	if os.access('ts_SYNTH',os.F_OK):
		shutil.rmtree('ts_SYNTH')
	os.mkdir('ts_SYNTH')
	if os.access('ts_GF',os.F_OK):
		shutil.rmtree('ts_GF')
	shutil.copytree('GF','./ts_GF')
		
	ts_opt   = tsini
	valbest  = rmsini
	ts_opt2  = ts1
	valbest2 = 1.1e10	
	for j in xrange(Nit):
		tmp_table = open('_tmp_ts_table', 'w')
		sts = Sts[j]
		if j>0:
			if (ts_opt2 <= ts_opt):
				ts1 = ts_opt2 - sts/2
				ts2 = ts_opt  + sts/2
			elif(ts_opt2 > ts_opt):
				ts1 = ts_opt  - sts/2
				ts2 = ts_opt2 + sts/2
			if ts1 < 1.:
				ts1 += 2.
		fid.write('iteration %d (%f<=ts<=%f)\n'% (j+1,ts1,ts2))
		its = 0
		for ts in pyl.arange(ts1,ts2+sts,sts):
			flag = 0
			eq.wcmtfile(cmttmp,ts,ts)
			os.system(REPREPARE_TS+'> /dev/null')
			os.system(WPINV_TS+' > /dev/null')
			out  = grep(r'^W_cmt_err:', 'LOG/_ts_wpinversion.log')			
			rms  = float(out[0].strip('\n').split()[1])
			nrms = float(out[0].strip('\n').split()[2])
			tmp_table.write(format%(its, ts, lat, lon, dep, rms, nrms))
			tmp_table.flush()
			if rms < valbest:
				ts_opt2  = ts_opt
				valbest2 = valbest
				ts_opt   = ts
				valbest  = rms
			elif rms < valbest2:
				ts_opt2  = ts
				valbest2 = rms
			fid.write('   ts=hd = %4.1f sec, rms = %12.7f mm\n'% (ts,rms))
			its += 1
		tmp_table.close()
		tmp_table = open('_tmp_ts_table', 'r')
		out_table = open('grid_search_ts_out%d'%(j+1),   'w')
		out_table.write('%5.1f%12.7f\n'%(ts_opt, valbest))
		out_table.write(tmp_table.read())
		out_table.close()
		tmp_table.close()
		os.remove('_tmp_ts_table')
		fid.write('   after iteration %d : ts_opt=%4.1f sec rms =%12.7f mm\n'%(j+1,ts_opt, valbest))
		
	fid.write('\nFinal Optimum values: time_shift (=half_duration) =  %5.1f   rms = %12.7f mm\n'%(ts_opt, valbest))
	eq.wcmtfile(cmttmp,ts_opt,ts_opt)
	os.system(REPREPARE_TS)
	if flagref:
		addrefsol(cmtpde,cmttmp)
		os.system(WPINV_TS+' -ref')
	else:
		os.system(WPINV_TS)

	######### USE MATPLOTLIB #########
	lab=['bo','b+','r.']
	for j in xrange(0,Nit):
		file = 'grid_search_ts_out%d' % (j+1)
		plot_ts(file,lab[j])
	pyl.plot([tsini],[rmsini],'rv')
	pyl.plot([ts_opt],[valbest],'bv')
	pyl.grid('on')
	pyl.savefig('grid_search_ts.png')
	pyl.close()
	return [ts_opt,ts_opt]

def fast_grid_search_ts(datdir,dmin,dmax,cmtref,ftable,eq,tsini,hdini,wpwin=[15.],flagref=0,out='stdout'):
	if out == 'stdout':
		fid = sys.stdout
	else:
		fid = open(out,'w')

	out     = grep(r'^W_cmt_err:', 'LOG/wpinversion.log')			
	rmsini  = float(out[0].strip('\n').split()[1])
	nrmsini = float(out[0].strip('\n').split()[2])

	fid.write('FAST CENTROID TIME DELAY GRID SEARCH\n')
	Nit = 3
	Sts = [4.,4.,2.,1.]
	
 	# if eq.mag <= 7.0:
	#	ts1 =  4. - tsini
	#	ts2 = 20. - tsini
	# elif eq.mag <8.0:
	#      	ts1 =  8. - tsini
	#  	ts2 = 48. - tsini
	# else: 
	#  	ts1 = 14. - tsini
	#	ts2 = 56. - tsini
	
	ts1 = 1. - tsini
	ts2 = tsini*2.	
	if ts2 > 60. - tsini:
		ts2 = 60. - tsini

	lat = eq.lat
	lon = eq.lon
	dep = eq.dep
	
	cmttmp = cmtref+'_ts_tmp'
	eq.wimaster(datdir,dmin,dmax,ftable,cmttmp,'ts_i_master','ts_GF',wpwin)
	eq.wcmtfile(cmttmp,tsini,hdini)
	
	if os.access('ts_SYNTH',os.F_OK):
		shutil.rmtree('ts_SYNTH')
	os.mkdir('ts_SYNTH')
	if os.access('ts_GF',os.F_OK):
		shutil.rmtree('ts_GF')
	shutil.copytree('GF','./ts_GF')

	ts_opt  = 0.
	valbest = rmsini
	ts_opt2 = ts1 - tsini
	valbest2 = 1.1e10
	os.system(REPREPARE_TS+'> /dev/null')
	for j in xrange(Nit):
		tmp_table  = open('_tmp_ts_table', 'w')
		sts = Sts[j]
		if j>0:
			if (ts_opt2 <= ts_opt):
				ts1 = ts_opt2 - sts/2
				ts2 = ts_opt  + sts/2 
			elif(ts_opt2 > ts_opt):
				ts1 = ts_opt  - sts/2 
				ts2 = ts_opt2 + sts/2 
			if ts1 < (1. - tsini):
				ts1 += abs(2. - tsini)
		fid.write('iteration %d (%f<=ts<=%f)\n'% (j+1,ts1+tsini,ts2+tsini))
		format     = '%02d %8.2f %8.2f %8.2f %8.2f %12.8f %12.8f\n'
		its = 0
		for ts in pyl.arange(ts1,ts2+sts,sts):
			flag = 0
			os.system(WPINV_TS+' -dts %4.1f > /dev/null'% ts)
			out  = grep(r'^W_cmt_err:', 'LOG/_ts_wpinversion.log')			
			rms  = float(out[0].strip('\n').split()[1])
			nrms = float(out[0].strip('\n').split()[2])
			tmp_table.write(format%(its, ts+tsini, lat, lon, dep, rms, nrms))
			tmp_table.flush()
			if rms < valbest:
				ts_opt2  = ts_opt
				valbest2 = valbest
				ts_opt   = ts
				valbest  = rms
			elif rms < valbest2:
				ts_opt2  = ts
				valbest2 = rms
			fid.write('   ts = %4.1f sec, rms = %12.7f mm\n'% (ts+tsini,rms))
			its += 1
		tmp_table.close()
		tmp_table = open('_tmp_ts_table', 'r')
		out_table = open('grid_search_ts_out%d'%(j+1),   'w')
		out_table.write('%5.1f%12.7f\n'%(ts_opt+tsini, valbest))
		out_table.write(tmp_table.read())
		out_table.close()
		tmp_table.close()
		os.remove('_tmp_ts_table')
		fid.write('   after iteration %d : ts_opt=%4.1f sec rms =%12.7f mm\n'%(j+1,ts_opt+tsini, valbest))
		
	fid.write('\nFinal Optimum values: time_shift =  %5.1f   rms = %12.7f mm\n'%(ts_opt+tsini, valbest))

	if os.access('ts_SYNTH',os.F_OK):
		shutil.rmtree('ts_SYNTH')
	os.mkdir('ts_SYNTH')
	eq.wcmtfile(cmttmp,ts_opt+tsini,ts_opt+tsini)
	os.system(REPREPARE_TS)
	if flagref:
		addrefsol(cmtpde,cmttmp)
		os.system(WPINV_TS+' -ref')
	else:
		os.system(WPINV_TS)
		

	######### USE MATPLOTLIB #########
	#lab=['bo','r.','r+','k.']
	lab=['bo','bo','bo','bo']
	for j in xrange(0,Nit):
		file = 'grid_search_ts_out%d' % (j+1)
		plot_ts(file,lab[j])
	pyl.plot([tsini],[rmsini],'rv')
	pyl.plot([ts_opt+tsini],[valbest],'bv')
	pyl.grid('on')
	pyl.savefig('grid_search_ts.png')
	pyl.close()
	return [ts_opt+tsini,ts_opt+tsini]


def usage():
	print 'usage: wp_grid_search [-f] [-t] [-p] [-i] ... [--help]'

def disphelp():
	print 'Centroid time-shift and centroid position grid search\n'
	usage()
	print '\nAll parameters are optional:'
	print '   -f, --fast           use a fast time-shift search'
	print '   -t, --onlyts         centroid time-shift grid search only'
	print '   -p, --onlyxy         centroid position grid search only'
	print '   -i, --imas \'file\'    set i_master file (i_master)'
	print '   -r, --ref            read the reference solution in cmtfile (no ref. sol.)'
	print '\n   -h, --help           display this help and exit'
	print '\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>'

##### MAIN #####	
if __name__ == "__main__":
	try:
		opts, args = getopt.gnu_getopt(sys.argv[1:],'ftpi:rh',["fast","onlyts","onlyxy","imas=","--ref","help"])
	except getopt.GetoptError, err:
		print '*** ERROR ***'
		print str(err)
		usage()
		sys.exit(1)
	
	fastflag = 0	
	flagts   = 1
	flagxy   = 1
	flagref  = 0
	i_master = 'i_master' 
	for o, a in opts:
		if o == '-h' or o == '--help':
			usage()
			disphelp()
			sys.exit(0)
		if o == '-f' or o == '--fast':
			fastflag = 1
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
		if o == '-i' or o == '--imas':
			i_master = a
		if o == '-r' or o == '--ref':
			flagref = 1
			
	out    = grep2([r'^SEED',r'^CMTFILE',r'^filt_cf1',r'^filt_cf2',\
				 r'^WP_WIN'], i_master)
 	dat    = out[0].replace(':','').strip('\n').split()[1]
 	cmtpde = out[1].replace(':','').strip('\n').split()[1]
	ftable = []
 	ftable.append(float(out[2].replace(':','').strip('\n').split()[1]))
 	ftable.append(float(out[3].replace(':','').strip('\n').split()[1]))
	wpwin  = map(float,out[4].replace(':','').strip('\n').split()[1:])
	
	try:
		out    = grep(r'^DMIN', i_master)
		dmin   = float(out[0].replace(':','').strip('\n').split()[1])
	except:
		dmin = 5.
   
	try:
		out    = grep(r'^DMAX', i_master)
		dmax   = float(out[0].replace(':','').strip('\n').split()[1])
	except:
		dmax = 88.

 	eq   = EarthQuake()
 	eq.rcmtfile(cmtpde)

 	if flagts == 1:
 		if fastflag == 1:
 			[eq.ts,eq.hd]=fast_grid_search_ts(dat,dmin,dmax,cmtpde,ftable,eq,eq.ts,eq.hd,wpwin,flagref)
 		else:
 			[eq.ts,eq.hd]=grid_search_ts(dat,dmin,dmax,cmtpde,ftable,eq,eq.ts,eq.hd,wpwin,flagref)
 	if flagxy == 1:
 		grid_search_xy(dat,dmin,dmax,cmtpde,ftable,eq,eq.ts,eq.hd,wpwin,flagref)

#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

######################################
# GRID SEARCH FOR WPHASE INVERSION
###
# Z.Duputel, L.Rivera and H.Kanamori
# 2009/07/15 -- initial version for grid search
# 2009/07/19 -- optimization for time shift : B-tree sampling method
# 2009/07/26 -- optimization for centroïd position : Oct-tree sampling method
# 2009/09/09 -- plot routines are now in a separate script

import os,re,shutil,sys,time,getopt
from EQ import *
from make_grids import *

WPHOME = os.path.expandvars('$WPHASE_HOME')
if WPHOME[-1] != '/':
	WPHOME += '/'

BIN = WPHOME+'bin/'

REPREPARE_TS = BIN+'reprepare_wp_ts.csh'
WPINV_TS     = BIN+'wpinversion -imas ts_i_master -ifil o_wpinversion -ofil ts_o_wpinversion -ocmtf ts_WCMTSOLUTION -ps ts_p_wpinversion '+\
	       '-wpbm ts_wpinv.pgm -log LOG/_ts_wpinversion.log -osyndir ts_SYNTH -pdata ts_fort.15 -nt '

RECALCSYN_XY = BIN+'recalc_fast_synths.csh'
REPREPARE_XY = BIN+'reprepare_wp_xy.csh'
WPINV_XY     = BIN+'wpinversion -imas xy_i_master -ifil o_wpinversion -ofil xy_o_wpinversion -ocmtf xy_WCMTSOLUTION -ps xy_p_wpinversion '+\
	       '-wpbm xy_wpinv.pgm -log LOG/_xy_wpinversion.log -osyndir xy_SYNTH -pdata xy_fort.15 -nt'


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

def grid_search_xy(datdir,cmtref,ftable,eq,ts,hd,wpwin=[15.],flagref=0,dmin=0.,dmax=90.,fileout='stdout'):
	if fileout == 'stdout':
		fid = sys.stdout
		flag = 0
	else:
		fid = open(fileout,'w')
		flag = 1

	o_file = 'grid_search_xy_out'
	if os.access(o_file,os.F_OK):
		os.remove(o_file)
# 	out     = grep(r'^W_cmt_err:', 'LOG/wpinversion.log')			
# 	rmsini  = float(out[0].strip('\n').split()[1])
# 	nrmsini = float(out[0].strip('\n').split()[2])
	format = '%03d %03d %8.2f %8.2f %8.2f %8.2f %8.2f %12.7f %12.7f\n'
	
	fid.write('CENTROID POSITION GRID SEARCH\n')

	# Initialize variables
 	eq_gs = EarthQuake()
 	EQcopy(eq_gs,eq)
	
	Nit  = 2
	dx   = 0.8
	lat1 = eq.lat - 2.4
	lat2 = eq.lat + 2.4
	lon1 = eq.lon - 2.4
	lon2 = eq.lon + 2.4

	#n    = (lat2-lat1)/dx+1#######################################

	lat = lat1
	coor = []
	while lat <= lat2:
		lon = lon1
		while lon <= lon2:
			coor.append([lat,lon])
			lon += dx
		lat += dx
	
	cmttmp = cmtref+'_xy_tmp'
	eq.wimaster(datdir,ftable,cmttmp,'xy_i_master',dmin,dmax,'./xy_GF/',wpwin) 	
	if os.access('xy_SYNTH',os.F_OK):
		shutil.rmtree('xy_SYNTH')
	if os.access('xy_DATA',os.F_OK):
		shutil.rmtree('xy_DATA')
	os.mkdir('xy_SYNTH')
	os.mkdir('xy_DATA')

	Nopt   = [4,4,2,1]
	rmsopt = []
	latopt = []
	lonopt = []
	for i in xrange(max(Nopt)):
		rmsopt.append(1.e10)
		latopt.append(lat2)
		lonopt.append(lon1)
# 	rmsopt[0] = rmsini ;
# 	rmsopt[0] = rmsini ;
# 	latopt[0] = eq.lat ;
# 	latopt[0] = eq.lon ;
	ncel = 0
	tmp_table  = open('_tmp_xy_table', 'w')	
 	for it in xrange(Nit):
 		if it != 0:
 			dx = dx/2.
 			coor = []
 			for i in xrange(Nopt[it]):
 				add_coor(coor,latopt[i]+dx,lonopt[i]-dx)
 				add_coor(coor,latopt[i]+dx,lonopt[i])
 				add_coor(coor,latopt[i]+dx,lonopt[i]+dx)
 				add_coor(coor,latopt[i]   ,lonopt[i]+dx)
 				add_coor(coor,latopt[i]-dx,lonopt[i]+dx)
 				add_coor(coor,latopt[i]-dx,lonopt[i])
 				add_coor(coor,latopt[i]-dx,lonopt[i]-dx)
 				add_coor(coor,latopt[i]   ,lonopt[i]-dx)

 		fid.write('Iteration %d:\n' % (it+1))
 		for cds in coor:
 			eq_gs.lat = cds[0]
 			eq_gs.lon = cds[1]
 			eq_gs.wcmtfile(cmttmp,ts,hd)
			os.system(RECALCSYN_XY+' > LOG/_log_py_recalsyn_xy')
			os.system(REPREPARE_XY+' > LOG/_log_py_reprepare_xy')
			os.system(WPINV_XY+' > LOG/_log_py_wpinv_xy')
			out  = grep(r'^W_cmt_err:', 'LOG/_xy_wpinversion.log')
			rms  = float(out[0].strip('\n').split()[1])
			nrms = float(out[0].strip('\n').split()[2])
			fid.write('   cell %3d : lat=%8.3fdeg lon=%8.3fdeg, rms = %12.7f mm\n'% (ncel+1,eq_gs.lat,eq_gs.lon,rms))
			for i in xrange(Nopt[it]):
				if rms < rmsopt[i]:
					for j in xrange(Nopt[it]-1,i-1,-1):
						rmsopt[j] = rmsopt[j-1]
						latopt[j] = latopt[j-1]
						lonopt[j] = lonopt[j-1]
					rmsopt[i] = rms
					latopt[i] = eq_gs.lat
					lonopt[i] = eq_gs.lon
					break
			ncel += 1
			tmp_table.write(format%(-99,-99,ts,hd,eq_gs.lat,eq_gs.lon,eq_gs.dep,rms,nrms))
			tmp_table.flush()
		fid.write('Optimum centroid location: %8.3f %8.3f;  rms = %12.7f mm\n'%(latopt[0], lonopt[0], rmsopt[0]))

	tmp_table.close()
	tmp_table = open('_tmp_xy_table', 'r')
	out_table = open(o_file, 'w')
	out_table.write('%8.3f %8.3f %12.7f\n'%(latopt[0], lonopt[0], rmsopt[0]))
	out_table.write('%8.3f %8.3f %12.7f\n'%(   eq.lat,    eq.lon, rmsopt[0]))	
       	out_table.write(tmp_table.read())
       	out_table.close()
       	tmp_table.close()
       	os.remove('_tmp_xy_table')

	eq_gs.lat = latopt[0]
	eq_gs.lon = lonopt[0]
	eq_gs.wcmtfile(cmttmp,ts,hd)
	os.system(RECALCSYN_XY+'> LOG/_log_py_recalsyn_xy')
	os.system(REPREPARE_XY+'> LOG/_log_py_reprepare_xy')
	if flagref:
		addrefsol(cmtpde,cmttmp)
		if flag:
			fid.close()
			os.system(WPINV_XY+' -ref >> '+fileout)
		else:
			os.system(WPINV_XY+' -ref')

	else:
		if flag:
			fid.close()
			os.system(WPINV_XY+' >> '+fileout)
		else:
			os.system(WPINV_XY)	

	# Set Mww
	out  = grep(r'^Wmag:', 'LOG/_xy_wpinversion.log')
	eq.mag = float(out[0].split()[1]) ;


def grid_search_ts(datdir,cmtref,ftable,eq,tsini,hdini,wpwin=[15.],flagref=0,dmin=0.,dmax=90.,fileout='stdout'):
	if fileout == 'stdout':
		fid = sys.stdout
		flag = 0
	else:
		fid = open(fileout,'w')
		flag = 1

	o_file = 'grid_search_ts_out'
	if os.access(o_file,os.F_OK):
		os.remove(o_file)
	out     = grep(r'^W_cmt_err:', 'LOG/wpinversion.log')			
	rmsini  = float(out[0].strip('\n').split()[1])
	nrmsini = float(out[0].strip('\n').split()[2])
	format     = '%02d %8.2f %8.2f %8.2f %8.2f %12.8f %12.8f\n'

	fid.write('CENTROID TIME DELAY GRID SEARCH\n')
	Nit = 3
	Sts = [4.,4.,2.,1.]

	if eq.mag < 5.5:
		ts1 = 1.
		ts2 = tsini*3.
		if ts2 > 60.:
			ts2 = 60.	
	else:
		
		if eq.mag <= 7.0:
			ts1 =  1. 
			ts2 = 30. 
		elif eq.mag <8.0:
			ts1 =  8. 
			ts2 = 48. 
		else: 
			ts1 = 14. 
			ts2 = 56. 	
				
	lat = eq.lat
	lon = eq.lon
	dep = eq.dep
	
	cmttmp = cmtref+'_ts_tmp'
	eq.wimaster(datdir,ftable,cmttmp,'ts_i_master',dmin,dmax,'ts_GF',wpwin)
	
	if os.access('ts_SYNTH',os.F_OK):
		shutil.rmtree('ts_SYNTH')
	os.mkdir('ts_SYNTH')
	if os.access('ts_GF',os.F_OK):
		shutil.rmtree('ts_GF')
	shutil.copytree('GF','./ts_GF')
		
	tsopt   = 0.
	rmsopt  = rmsini
	tsopt2  = ts1
	rmsopt2 = 1.1e10
	its = 0
	tmp_table = open('_tmp_ts_table', 'w')
	for j in xrange(Nit):
		sts = Sts[j]
		if j>0:
			if (tsopt2 <= tsopt):
				ts1 = tsopt2 - sts/2
				ts2 = tsopt  + sts/2
			elif(tsopt2 > tsopt):
				ts1 = tsopt  - sts/2
				ts2 = tsopt2 + sts/2
			if ts1 < 1.:
				ts1 += 2.
		fid.write('iteration %d (%f<=ts<=%f)\n'% (j+1,ts1,ts2))
		ts = ts1
		while ts < ts2+sts:
			eq.wcmtfile(cmttmp,ts,ts)
			os.system(REPREPARE_TS+' > LOG/_log_py_reprepare_ts')
			os.system(WPINV_TS+' > LOG/_log_py_wpinv_ts')
			out  = grep(r'^W_cmt_err:', 'LOG/_ts_wpinversion.log')			
			rms  = float(out[0].strip('\n').split()[1])
			nrms = float(out[0].strip('\n').split()[2])
			tmp_table.write(format%(its, ts, lat, lon, dep, rms, nrms))
			tmp_table.flush()
			if rms < rmsopt:
				tsopt2  = tsopt
				rmsopt2 = rmsopt
				tsopt   = ts
				rmsopt  = rms
			elif rms < rmsopt2:
				tsopt2  = ts
				rmsopt2 = rms
			fid.write('   ts=hd = %4.1f sec, rms = %12.7f mm\n'% (ts,rms))
			its += 1
			ts += sts
		fid.write('   after iteration %d : tsopt=%4.1f sec rms =%12.7f mm\n'%(j+1,tsopt, rmsopt))
	fid.write('\nFinal Optimum values: time_shift (=half_duration) =  %5.1f   rms = %12.7f mm\n'%(tsopt, rmsopt))		
	tmp_table.close()
	tmp_table = open('_tmp_ts_table','r')
	out_table = open(o_file,'w')
	out_table.write('%5.1f%12.7f\n'%(tsopt, rmsopt))	
	out_table.write('%5.1f%12.7f\n'%( tsini, rmsini))
	out_table.write(tmp_table.read())
	out_table.close()
	tmp_table.close()
	os.remove('_tmp_ts_table')
	
	eq.wcmtfile(cmttmp,tsopt,tsopt)
	os.system(REPREPARE_TS)
	if flagref:
		addrefsol(cmtpde,cmttmp)
		if flag:
			fid.close()	
			os.system(WPINV_TS+' -ref >> '+fileout)
		else:
			os.system(WPINV_TS+' -ref')
	else:
		if flag:
			fid.close()	
			os.system(WPINV_TS+' >> '+fileout)
		else:
			os.system(WPINV_TS)

	# Set Mww
	out  = grep(r'^Wmag:', 'LOG/_ts_wpinversion.log')
	eq.mag = float(out[0].split()[1]) ;

	return [tsopt,tsopt]

def fast_grid_search_ts_old(datdir,cmtref,ftable,eq,tsini,hdini,wpwin=[15.],flagref=0,dmin=0.,dmax=90.,out='stdout'):
	if out == 'stdout':
		fid = sys.stdout
		flag = 0
	else:
		fid = open(out,'w')
		flag = 1

	o_file = 'grid_search_ts_out'
	if os.access(o_file,os.F_OK):
		os.remove(o_file)
		
	out     = grep(r'^W_cmt_err:', 'LOG/wpinversion.log') 
	rmsini  = float(out[0].strip('\n').split()[1])
	nrmsini = float(out[0].strip('\n').split()[2])
	format     = '%02d %8.2f %8.2f %8.2f %8.2f %12.8f %12.8f\n'
	
	fid.write('FAST CENTROID TIME DELAY GRID SEARCH\n')
	Nit = 3
	Sts = [4.,4.,2.,1.]
	
	if eq.mag < 5.5:
		ts1 = 1. - tsini
		ts2 = tsini*3.	
		if ts2 > 60. - tsini:
			ts2 = 60. - tsini
	else:
		if eq.mag <= 7.0:
			ts1 =  1. 
			ts2 = 30. 
		elif eq.mag <8.0:
			ts1 =  8. 
			ts2 = 48. 
		else: 
			ts1 = 14. 
			ts2 = 56. 


	lat = eq.lat
	lon = eq.lon
	dep = eq.dep
	
	cmttmp = cmtref+'_ts_tmp'
	eq.wimaster(datdir,ftable,cmttmp,'ts_i_master',dmin ,dmax,'ts_GF',wpwin)
	eq.wcmtfile(cmttmp,tsini,hdini)
	
	if os.access('ts_SYNTH',os.F_OK):
		shutil.rmtree('ts_SYNTH')
	os.mkdir('ts_SYNTH')
	if os.access('ts_GF',os.F_OK):
		shutil.rmtree('ts_GF')
	shutil.copytree('GF','./ts_GF')

	tsopt  = 0.
	rmsopt = rmsini
	tsopt2 = ts1 - tsini
	rmsopt2 = 1.1e10
	os.system(REPREPARE_TS+'> LOG/_log_py_reprepare_ts')
	its = 0
	tmp_table  = open('_tmp_ts_table', 'w')
	
	for j in xrange(Nit):
		sts = Sts[j]
		if j>0:
			if (tsopt2 <= tsopt):
				ts1 = tsopt2 - sts/2
				ts2 = tsopt  + sts/2 
			elif(tsopt2 > tsopt):
				ts1 = tsopt  - sts/2 
				ts2 = tsopt2 + sts/2 
			if ts1 < (1. - tsini):
				ts1 += abs(2. - tsini)
		fid.write('iteration %d (%f<=ts<=%f)\n'% (j+1,ts1+tsini,ts2+tsini))
		ts = ts1
		while ts < ts2+sts:
			os.system(WPINV_TS+' -dts %4.1f > LOG/_log_py_wpinv_ts'% ts)
			out  = grep(r'^W_cmt_err:', 'LOG/_ts_wpinversion.log')			
			rms  = float(out[0].strip('\n').split()[1])
			nrms = float(out[0].strip('\n').split()[2])
			tmp_table.write(format%(its, ts+tsini, lat, lon, dep, rms, nrms))
			tmp_table.flush()
			if rms < rmsopt:
				tsopt2  = tsopt
				rmsopt2 = rmsopt
				tsopt   = ts
				rmsopt  = rms
			elif rms < rmsopt2:
				tsopt2  = ts
				rmsopt2 = rms
			fid.write('   ts = %4.1f sec, rms = %12.7f mm\n'% (ts+tsini,rms))
			its += 1
			ts += sts
		fid.write('   after iteration %d : tsopt=%4.1f sec rms =%12.7f mm\n'%(j+1,tsopt+tsini, rmsopt))
	fid.write('\nFinal Optimum values: time_shift =  %5.1f   rms = %12.7f mm\n'%(tsopt+tsini, rmsopt))
	tmp_table.close()
	tmp_table = open('_tmp_ts_table', 'r')
	out_table = open(o_file,'w')
	out_table.write('%5.1f%12.7f\n'%(tsopt+tsini, rmsopt))	
	out_table.write('%5.1f%12.7f\n'%(tsini, rmsini))
	out_table.write(tmp_table.read())
	out_table.close()
	tmp_table.close()
	os.remove('_tmp_ts_table')

	eq.wcmtfile(cmttmp,tsopt+tsini,tsopt+tsini)
	os.system(REPREPARE_TS)
	if flagref:
		addrefsol(cmtpde,cmttmp)
		os.system(WPINV_TS+' -ref')
	else:
		os.system(WPINV_TS)
		

	if flag:
		fid.close()

	# Set Mww
	out  = grep(r'^Wmag:', 'LOG/_ts_wpinversion.log')
	eq.mag = float(out[0].split()[1]) ;

	return [tsopt+tsini,tsopt+tsini]


def fast_grid_search_ts(datdir,cmtref,ftable,eq,tsini,hdini,wpwin=[15.],flagref=0,dmin=0.,dmax=90.,fileout='stdout'):
	if fileout == 'stdout':
		fid  = sys.stdout
		flag = 0
	else:
		fid = open(fileout,'w')
		flag = 1

	o_file = 'grid_search_ts_out'
	if os.access(o_file,os.F_OK):
		os.remove(o_file)

	out     = grep(r'^W_cmt_err:', 'LOG/wpinversion.log') 
	rmsini  = float(out[0].strip('\n').split()[1])
	nrmsini = float(out[0].strip('\n').split()[2])
	format  = '%02d %8.2f %8.2f %8.2f %8.2f %12.8f %12.8f\n'
	fid.write('FAST CENTROID TIME DELAY GRID SEARCH\n')		
	Nit = 3
	sts = 4.
	if eq.mag < 5.5:
		ts1 = 1.
		ts2 = tsini*3.	
		if ts2 > 60.:
			ts2 = 60. 
	else:
		if eq.mag <= 7.0:
			ts1 =  1. 
			ts2 = 30. 
		elif eq.mag <8.0:
			ts1 =  8. 
			ts2 = 48. 
		else: 
			ts1 = 14. 
			ts2 = 56. 
	fid.write('  ts1 = %5.1f sec, step = %5.1f sec, ts2 = %5.1f sec \n'%(ts1,sts,ts2))  	
	lat = eq.lat
	lon = eq.lon
	dep = eq.dep

	print lat,lon,dep

	cmttmp = cmtref+'_ts_tmp'
	eq.wimaster(datdir,ftable,cmttmp,'ts_i_master',dmin ,dmax,'ts_GF',wpwin)
	eq.wcmtfile(cmttmp,tsini,hdini)
	
	if os.access('ts_SYNTH',os.F_OK):
		shutil.rmtree('ts_SYNTH')
	os.mkdir('ts_SYNTH')
	if os.access('ts_GF',os.F_OK):
		shutil.rmtree('ts_GF')
	shutil.copytree('GF','./ts_GF')
	if flag:
		fid.close()
		os.system(WPINV_TS+' -ts %4.1f %4.1f %4.1f -Nit 3 -ogsf %s -ifil o_wpinversion >> %s'% (ts1,sts,ts2,o_file,fileout))
	else:
		os.system(WPINV_TS+' -ts %4.1f %4.1f %4.1f -Nit 3 -ogsf %s -ifil o_wpinversion'% (ts1,sts,ts2,o_file))
	tmp_table = open(o_file, 'r')
	tsopt, rmsopt = map(float,tmp_table.readline().strip('\n').split())
 	tmp_table.close()

 	eq.wcmtfile(cmttmp,tsopt,tsopt)
 	os.system(REPREPARE_TS)
 	if flagref:
 		addrefsol(cmtpde,cmttmp)
		if flag:
			os.system(WPINV_TS+' -ref >> '+fileout)
		else:
			os.system(WPINV_TS+' -ref')
 	else:
		if flag:
			os.system(WPINV_TS+' >> '+fileout)
		else:
			os.system(WPINV_TS)

	# Set Mww
	out  = grep(r'^Wmag:', 'LOG/_ts_wpinversion.log')
	eq.mag = float(out[0].split()[1]) ;
	
	return [tsopt,tsopt]



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
		opts, args = getopt.gnu_getopt(sys.argv[1:],'ftpi:rh',["fast","onlyts","onlyxy","imas=","ref","help"])
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
		dmin   = 0.
	try:
		out    = grep(r'^DMAX', i_master)
		dmax   = float(out[0].replace(':','').strip('\n').split()[1])
	except:
		dmax   = 90.
		
 	eq   = EarthQuake()
 	eq.rcmtfile(cmtpde)

 	if flagts == 1:
 		if fastflag == 1:
 			[eq.ts,eq.hd]=fast_grid_search_ts(dat,cmtpde,ftable,eq,eq.ts,eq.hd,wpwin,flagref,dmin,dmax)
 		else:
 			[eq.ts,eq.hd]=grid_search_ts(dat,cmtpde,ftable,eq,eq.ts,eq.hd,wpwin,flagref,dmin,dmax)
 	if flagxy == 1:
 		grid_search_xy(dat,cmtpde,ftable,eq,eq.ts,eq.hd,wpwin,flagref,dmin,dmax)

#!/usr/bin/env python
# *-* coding: iso-8859-1 *-*

############################################################################
#
#	              W phase source inversion package 	            
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

# GRID SEARCH FOR WPHASE INVERSION

# Time-shift grid-search parameters
TS_NIT   = 3  # Nb of iterations
TS_DT    = 4. # Initial time step
TSBOUNDS = [] # Bounds (empty=automatically determined from mb or Ms in the PDE line)
TS_OFILE = 'grid_search_ts_out'

# Centroid Lat/Lon grid-search parameters
XY_NIT   = 3   # Nb of iterations
XY_DX    = 0.4 # Intial samp. period
XY_NX    = 3   # Half_width = XY_NX*XY_DX
XY_NOPT  = 5   # Nb of optimal-points
XY_OFILE = 'grid_search_xy_out'

# Centroid Depth grid-search parameters
XYZ_NIT   = 1    # Nb of iterations
XYZ_DX    = 0.6  # Intial samp. period
XYZ_NX    = 1    # Half_width = XYZ_NX*XYZ_DX (if XYZ_NX=0: no Lat/Lon grid-seach is performed)
XYZ_NOPT  = 4    # Nb of optimal-points
DDEP      = 50.  # Delta depth ( Z_SEARCH within Z_INITIAL +/- DDEP )
MINDEP    = 11.5 
XYZ_OFILE = 'grid_search_xyz_out'

import os,re,shutil,sys,time,getopt
from EQ import *

WPHOME = os.path.expandvars('$WPHASE_HOME')
print 'WPHASE_HOME is %s'%(WPHOME)
if WPHOME[-1] != '/':
	WPHOME += '/'

GF_PATH = os.path.expandvars('$GF_PATH')
print 'GF_PATH is %s'%(GF_PATH)

#VERSION = 'Version: '
#entfile = WPHOME+'.svn/entries'
#if os.path.exists(entfile):
#	VERSION += open(entfile).readlines()[3].strip()
VERSION = 'Version: r244N'

BIN = WPHOME+'bin/'

WPINV_XY     = BIN+'wpinversion_gs -imas i_master -ifil o_wpinversion'

  
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

def parse_config(cfg_file):
	config = {}
	try:
		config_lines = open(cfg_file, 'r').readlines()
		for line in config_lines:
			if line.find('#')==0:
				continue
			if line.rstrip():
				key,value = line.strip().split(':')
				config[key.strip()]=value.strip()
	except:
		sys.stderr.write('Error: format  %s\n'%cfg_file)
		sys.exit(1)
	return config

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

def rm(fd):
        if os.path.islink(fd) or os.path.isfile(fd):
                os.remove(fd)
                return 0
        elif os.path.isdir(fd):
                shutil.rmtree(fd)
                return 0
        return 1

def addslash(direc):
	if len(direc) > 0:
		if direc[-1] != '/':
			direc += '/'
	return direc

def grid_search(eq,cmtref,ts_Nit,ts_dt,tsb,xy_Nit,xy_dx,xy_Nx,xy_Nopt,fastflag,flagts,flagxy,sdrM0={},dz=0.,
		minz=3.5,ts_ofile='grid_search_ts_out',xy_ofile='grid_search_xy_out',stdoutput='stdout',
		logfile='LOG/gs_o_wpinversion.log', comments = []):
	if stdoutput == 'stdout':
		fid = sys.stdout
		flag = 0
	else:
		fid = open(stdoutput,'a+')
		flag = 1
	EXE = WPINV_XY		
	fid.write('CENTROID GRID SEARCH\n')
	# Setting parameters ########
	cmttmp = cmtref
	optpar = ' -log %s -osyndir gs_SYNTH -icmtf %s '%(logfile,cmtref)
	for o,a in sdrM0.items():
		if len(a):
			optpar += ' %s %s '%(o,a)
		else:
			optpar += ' %s '%(o)
	if not os.access('gs_SYNTH',os.F_OK):
		os.mkdir('gs_SYNTH')	
	# time-shift
	if flagts:
		if len(tsb) == 2:
			ts1 = tsb[0]
			ts2 = tsb[1]
		else:		
			if eq.mag < 5.5:
				ts1 = 1.
				ts2 = eq.ts*3.	
				if ts2 > 100.:
					ts2 = 100.	
			else:
				ts1 =  1. 
				if eq.mag <= 7.0:
					ts2 = 30. 
				elif eq.mag <= 8.0:
					ts2 = 100. 
				else: 
					ts2 = 168.
		optpar += ' -ts %10.4f %10.4f %10.4f -ts_Nit %d -otsgsf %s'%(ts1,ts_dt,ts2,ts_Nit,ts_ofile)
		if not fastflag:
			optpar += ' -hdsafe '
	else:
		optpar += ' -nots '
	if flagxy:
		optpar += ' -xy_Nit %d -dx %.2f -Nx %d -Nopt %d -oxygsf %s'%(xy_Nit,xy_dx,xy_Nx,xy_Nopt,xy_ofile)
		if dz>0.:
			optpar += ' -dz %.2f -minz %.2f '%(dz,minz)
		wcmtfile = 'xy_WCMTSOLUTION'		
	else:
		optpar += ' -noxy '
		if flagts:
			wcmtfile = 'ts_WCMTSOLUTION'
		else:
			wcmtfile = 'gs_WCMTSOLUTION'
	if flag:
		fid.close()
		optpar += ' > %s '%stdoutput
	for c in comments:
		optpar += ' -comments "'+c+'"'
	print 'Command_line:'+EXE+optpar
	
	os.system(EXE+optpar)
	# Update eq
	eq.rcmtfile(wcmtfile)
	out = grep(r'^Wmag:',logfile)
	eq.mag = float(out[-1].split()[1]) ;

def usage():
	print 'usage: wp_grid_search [-s] [-t] [-p] [-i] ... [--help]'

def disphelp():
	print 'Centroid time-shift and centroid position grid search\n'
	usage()
	print '\nAll parameters are optional:'
	print '   -s, --hdsafe         Use a  time grid-search considering ts=fd'
	print '   -t, --onlyts         Centroid time-shift grid search only'
	print '   -p, --onlyxy         Centroid position grid search only'
	print '   -S, --npar           Do not use the parallelized grid-search and use '
	print '                          the sequential version instead (parallelized version is used)'
	print '   -i, --imas \'file\'    Set i_master file (i_master)'
	print '   -n, --noref          Do not use the reference solution in cmtfile (ref. sol. used)'
	print '   --nont               Full moment tensor inversion (no null trace)'	
	print '   --dc                 Double-couple inversion'
	print '   --strike \'strike\'    Double-couple inversion with fixed strike'
	print '   --dip \'dip\'          Double-couple inversion with fixed dip'
	print '   --rake \'rake\'        Double-couple inversion with fixed rake'
	print '   --mom \'mom\'          Double-couple inversion with fixed scalar moment'
	print '\n   -h, --help           Display this help and exit'
	print '\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>'

##### MAIN #####	
if __name__ == "__main__":
	try:
		opts, args = getopt.gnu_getopt(sys.argv[1:],'stpSdi:nhz',["hdsafe","onlyts","onlyxy","npar",
									  "imas=","strike=","dc","nont","dip=",
									  "rake=","mom=","noref","xyz","old",
									  "help"])
	except getopt.GetoptError, err:
		print '*** ERROR ***'
		print str(err)
		usage()
		sys.exit(1)
	
	i_master = 'i_master' 
	fastflag = 1	
	flagts   = 1
	flagxy   = 1
	flagxyz  = 0
	flagref  = 1
	sdrM0    = {}
	for o, a in opts:
		if o == '-h' or o == '--help':
			disphelp()
			sys.exit(0)
		if o == '-s' or o == '--hdsafe':
			fastflag = 0
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
			flagts   = 0
			fastflag = 0
			flagxy = 1
		if o == '--dc':
			sdrM0['-dc']=''
		if o == '--nont':
			sdrM0['-nont']=''			
		if o == '--strike':
			sdrM0['-strike']=a
		if o == '--dip':
			sdrM0['-dip']=a
		if o == '--rake':
			sdrM0['-rake']=a
		if o == '--mom':
			sdrM0['-mom']=a
		if o == '-i' or o == '--imas':
			i_master = a
		if o == '-n' or o == '--noref':
			flagref = 0
		if o == '-z' or o == '--xyz':
			flagxyz = 1
		if o == '--old':
			WPINV_XY += ' -old'

	# Read i_master
	iconfig = parse_config(i_master)
	cmtref  = iconfig['CMTFILE']
	evname  = iconfig['EVNAME'].replace(' ','_').replace(',','')
	# Set comments
	Median    = '-med '
	if iconfig.has_key('P2P_SCREENING'):
		if iconfig['P2P_SCREENING'] != 'YES':
			Median = ' '
	ths = '5.0 3.0 0.9'
	if iconfig.has_key('RMS_SCREENING'):
		ths = iconfig['RMS_SCREENING']
	comments = [VERSION,'GF_PATH: '+GF_PATH,'Screening: '+Median+ths]
	# Read CMTFILE
 	eq   = EarthQuake()
 	eq.rcmtfile(cmtref)
	eq.title = evname.strip().replace(' ','_').replace(',','')
	cmtf = open(cmtref,'r')
	L=cmtf.readlines()
	cmtf.close()
	if len(L) < 13:
		print '*** WARNING : no reference solution in %s'%(cmtref)
		flagref = 0

	i_cmtfile = cmtref
	if (flagts or flagxy) and not flagxyz: # LAT/LON Grid-search
		grid_search(eq,i_cmtfile,TS_NIT,TS_DT,TSBOUNDS,XY_NIT,XY_DX,XY_NX,XY_NOPT,fastflag,
					flagts,flagxy,sdrM0,ts_ofile=TS_OFILE,xy_ofile=XY_OFILE,comments=comments)
	if flagxyz:                              # 3D Grid-search
		grid_search(eq,i_cmtfile,TS_NIT,TS_DT,TSBOUNDS,XYZ_NIT,XYZ_DX,XYZ_NX,XYZ_NOPT,fastflag,
					flagts,flagxyz,sdrM0,dz=DDEP,minz=MINDEP,ts_ofile=TS_OFILE,xy_ofile=XYZ_OFILE,
					comments=comments)
		if flagxy:
			eq.wcmtfile('_tmp_CMTSOLUTION.xyz')
			if flagref:
				addrefsol(cmtref,'_tmp_CMTSOLUTION.xyz')
			grid_search(eq,'_tmp_CMTSOLUTION.xyz',TS_NIT,TS_DT,TSBOUNDS,XY_NIT,XY_DX,XY_NX,XY_NOPT,
				    0,0,1,sdrM0,ts_ofile=TS_OFILE,xy_ofile=XY_OFILE,comments=comments)
			rm('_tmp_CMTSOLUTION.xyz')
	if os.path.exists('_tmp_ts_table'):		
		rm('_tmp_ts_table')
	if os.path.exists('_tmp_xy_table'):		
		rm('_tmp_xy_table')

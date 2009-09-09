#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

######################################
# TRIGGERING FOR WPHASE INVERSION
###
# Z.Duputel, L.Rivera and H.Kanamori
# 2009/07/08 -- initial version
# 2009/07/26 -- Centoid grid searches

import os,shutil,sys,time,calendar
import scipy as sp

from EQ import *
from wp_grid_search import *


WPHOME = os.path.expandvars('$WPHASE_HOME')
if WPHOME[-1] != '/':
	WPHOME += '/'

RUNQRT   = WPHOME+'run_qrt/'
BIN      = WPHOME+'bin/'

DATDIR   = '/home/zac/DATA/BUD/'
DATALESS = '/home/zac/DATA/dataless/'

RUNALL = BIN+'RUNA_qrt.csh >& RUNA_log'



def wpinversion_50(ftable,eq):

	dir = '%s%s_M%d_%s_50deg/'% (RUNQRT,time.strftime('%Y_%m_%d',eq.Otime),eq.mag,eq.title[7:].replace(' ','_').replace(',','_'))
	if os.access(dir,os.F_OK):
		shutil.rmtree(dir)
	os.mkdir(dir)

	cmtpde = 'CMTSOLUTION'

	# Half duration
	Mo = 10.**(1.5*eq.mag + 16.1)
	hd = 1.2*10**(-8.)*(Mo)**(1./3.)

	# Frequency band
	for i in xrange(len(ftable[0])):
		if eq.mag >= ftable[0][i] and eq.mag < ftable[1][i]:
			break;	
	freqs = [ftable[2][i],ftable[3][i]]
	
	eq.wcmtfile(dir+cmtpde,hd,hd)
	eq.wimaster(DATDIR,freqs,cmtpde,dir+'i_master',10.,50.,DATALESS=DATALESS)
	os.chdir(dir)
	eq.affiche('PDErsstrig50')
	os.system(RUNALL)
	[ts_opt,hd_opt]=fast_grid_search_ts(DATDIR,10.,50.,cmtpde,freqs,eq,hd,hd,out='ts_gs_log')
	grid_search_xy(dir,cmtpde,freqs,eq,ts_opt,hd_opt,10.,50.,out='xy_gs_log')
	os.chdir('../')


def wpinversion_90(ftable,eq):
	
	dir = '%s%s_M%d_%s_90deg/'% (RUNQRT,time.strftime('%Y_%m_%d',eq.Otime),eq.mag,eq.title[7:].replace(' ','_').replace(',','_'))
	if os.access(dir,os.F_OK):
		shutil.rmtree(dir)
	os.mkdir(dir)

	cmtpde = 'CMTSOLUTION'

	# Half duration
	Mo = 10.**(1.5*eq.mag + 16.1)
	hd = 1.2*10**(-8.)*(Mo)**(1./3.)

	# Frequency band
	for i in xrange(len(ftable[0])):
		if eq.mag >= ftable[0][i] and eq.mag < ftable[1][i]:
			break;	
	freqs = [ftable[2][i],ftable[3][i]]
	
	eq.wcmtfile(dir+cmtpde,hd,hd)
	eq.wimaster(DATDIR,10.,88.,freqs,cmtpde,dir+'i_master',10.,50.,DATALESS=DATALESS,)
	os.chdir(dir)
	eq.affiche('PDErsstrig90')
	os.system(RUNALL)
	[ts_opt,hd_opt]=fast_grid_search_ts(dir,10.,90.,cmtpde,freqs,eq,hd,hd,out='ts_gs_log')
	grid_search_xy(DATDIR,cmtpde,freqs,eq,ts_opt,hd_opt,10.,90.,out='xy_gs_log')
	os.chdir('../')		

def main(argv=None):

	# INPUT PARAMETERS #################################################
	mintime = time.struct_time((2009, 5, 27, 9, 40, 53, 2, 147, 0)) # minimum time
	
	db_fil_50 = RUNQRT+'EQ_DB_inv50'
	db_fil_90 = RUNQRT+'EQ_DB_inv90'
	minmag_inv   = 5.9
	
	db_fil    = RUNQRT+'EQ_DB4' # data base	
	minmag    = 5.0
	
	urlemsc   = 'http://www.emsc-csem.org/rss.php'                         # RSS feeds
	urlneic   = 'http://earthquake.usgs.gov/eqcenter/catalogs/7day-M5.xml' # RSS feeds
	
	#wallmsgfil = '/home/zac/DATA/read_rss/merge/wallmsg'
	
	O_scr  = 1 * 60 # Otime screening window (width in seconds)
	Ep_scr = 1      # Epicenter screening window (width in degrees)
	M_scr  = 1      # Magnitude screening window
	flag   = 1
	
	mmin   = [0.0   ,6.5  ,7.0    ,7.5    ,8.0]
	mmax   = [6.5   ,7.0  ,7.5    ,8.0    ,99.]
	fmin   = [0.0067,0.002,0.00167,0.00167,0.001]
	fmax   = [0.02  ,0.01 ,0.01   ,0.005  ,0.005]
	ftable = [mmin,mmax,fmin,fmax]

	######################################################################
	
	O_scr2  = O_scr**2
	Ep_scr2 = Ep_scr**2
	M_scr2  = M_scr**2
	
	eqs    = {}
	eqs_50 = {}
	eqs_90 = {}	
	
	# Reading databases
	eqs    = read_db(db_fil,mintime,minmag,O_scr,Ep_scr2,M_scr2,flag,eqs)
	eqs_50 = read_db(db_fil_50,mintime,minmag,O_scr,Ep_scr2,M_scr2,0,eqs_50)	
	eqs_90 = read_db(db_fil_90,mintime,minmag,O_scr,Ep_scr2,M_scr2,0,eqs_90)
	
	# if flag == 0 : screening by eventid (one entry per id)
	# if flag == 1 : screening by otime/Epicenter/Magnitude
	#                (one entry per earthquake).
	# if flag == 2 : as for flag==1 but allows pde updates.
	#                (database 'eqs' is useless?)
 	flag = 2 
 	while(flag):
 	    rd = sp.rand()
 	    if rd <= 0.5:
  	        eqs = r_emsc_feeds(urlemsc,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag,eqs)
 		eqs = r_neic_feeds(urlneic,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag,eqs)
 		flag = 0
 	    elif rd > 0.5:
 		eqs = r_neic_feeds(urlneic,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag,eqs)
		eqs = r_emsc_feeds(urlemsc,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag,eqs)		
 		flag = 0

	# Write databases
	fid = open(db_fil,'wt')
	fid.write(str(len(eqs))+'\n')
	eqsort  = sort_EQ(eqs,1)
	NOW     = time.gmtime()
	eqtoinv50=[]
	eqtoinv90=[]
	for eq in eqsort:
	    dbflag  = 0
	    eq.write(fid)
	    delta = (calendar.timegm(NOW) - calendar.timegm(eq.Otime))/60.
	    if eqs_90.has_key(eq.id):
	        print '%s with evid : %s already inverted' % (eq.title,eq.id)
		continue
	    elif eqs_50.has_key(eq.id):
	        print '%s with evid : %s already inverted for stations until 50 deg' % (eq.title,eq.id)
		dbflag = 1
	    if (eq.mag >= minmag_inv):
	        if (delta > 37. and dbflag != 1): #24+13min (13min=approximative acquisition delay)
		   eqs_50[eq.id] = EarthQuake()
		   EQcopy(eqs_50[eq.id],eq)
		   print 'to be inverted at 50 deg: %s (id:%s)'% (eq.title,eq.id)
		   eqtoinv50.append(eq)
		if (delta > 61.): #48+13min (13min=approximative acquisition delay)
		   eqs_90[eq.id] = EarthQuake()
		   EQcopy(eqs_90[eq.id],eq)
		   print 'to be inverted at 90 deg: %s (id:%s)'% (eq.title,eq.id)
		   eqtoinv90.append(eq)
	fid.close()

	# Invertion at 50 deg
	fid = open(db_fil_50,'wt')
	fid.write(str(len(eqs_50))+'\n')	
  	for EQ in eqs_50.values():
 		EQ.write(fid)
 	fid.close()

	# Invertion at 90 deg
	fid = open(db_fil_90,'wt')
	fid.write(str(len(eqs_90))+'\n')	
  	for EQ in eqs_90.values():
 		EQ.write(fid)
 	fid.close()	

	if len(eqtoinv50)>0 or len(eqtoinv90)>0:
		
		for eq in eqtoinv50:
			print '\nWPINVERSSION (50 deg): '
			eq.affiche()
			wpinversion_50(ftable,eq)

		for eq in eqtoinv90:
			print '\nWPINVERSSION (90 deg): '
			eq.affiche()
			wpinversion_90(ftable,eq)
		
	
if __name__ == "__main__":
	main(sys.argv)
	

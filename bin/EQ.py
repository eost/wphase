#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

import os, sys, urllib2
import re
import time, calendar

class EarthQuake:
	"Store Earthquake information"
 	def __init__(self):
		self.id    = '--'
 		self.title = '--'
		self.Otime = time.gmtime(0)
		self.date  = time.gmtime(0)
		self.mag   = 0.0		
		self.lat   = 0.0 
		self.lon   = 0.0
		self.dep   = 0.0
		self.org   = '--'
		self.ts    = 0.0
		self.hd    = 0.0		
	def affiche(self,out='stdout'):
		if out == 'stdout':
		    fid = sys.stdout
		else:
		    fid = open(out,'w')
		fid.write('%s \n  id : %s (from %s)\n' % (self.title,self.id,self.org))
		fid.write('  origin time : %s UTC\n' % time.strftime('%Y/%m/%d  --  %H:%M:%S',self.Otime))
		if calendar.timegm(self.date) != 0:
		   fid.write('  received    : %s UTC\n' % time.strftime('%Y/%m/%d  --  %H:%M:%S',self.date))
		   delta = calendar.timegm(self.date) - calendar.timegm(self.Otime)
		   m = int(delta/60.)
		   s = int((delta/60.-m)*60)
		   fid.write('  delay       : %d min %d sec\n' % (m,s))
		fid.write('  Lat : %.3f\n' % self.lat)
		fid.write('  Lon : %.3f\n' % self.lon)
		fid.write('  Dep : %.2f\n' % self.dep)
		if out != 'stdout':
		   fid.close()		
	def read(self,fid):
		self.id    = fid.readline().strip('\n')
		self.org   = fid.readline().strip('\n')
		self.mag   = float(fid.readline().strip('\n'))
		self.title = fid.readline().strip('\n')
		self.Otime = time.strptime(fid.readline().strip('\n'),'%Y/%m/%d  --  %H:%M:%S')
		self.date  = time.strptime(fid.readline().strip('\n'),'%Y/%m/%d  --  %H:%M:%S')
		self.lat   = float(fid.readline().strip('\n'))
		self.lon   = float(fid.readline().strip('\n'))
		self.dep   = float(fid.readline().strip('\n'))
		a=self.id[0]
		if a.isalpha():
		    self.org = 'neic'
	        else:
		    self.org = 'emsc'
		    
	def write(self,fid):
		fid.write('%s\n' % self.id)
		fid.write('%s\n' % self.org)
		fid.write('%f\n' % self.mag)		
		fid.write('%s\n' % self.title)
		fid.write('%s\n' % time.strftime('%Y/%m/%d -- %H:%M:%S',self.Otime))
		fid.write('%s\n' % time.strftime('%Y/%m/%d -- %H:%M:%S',self.date))
		fid.write('%.3f\n' % self.lat)
		fid.write('%.3f\n' % self.lon)
		fid.write('%.2f\n' % self.dep)
	def warnmsg(self,fil = 'wallmsg'):
		self.affiche(fil)
		os.system("wall < %s" % fil)
	def wcmtfile(self,cmtfil,ts=-99.,hd=-99.):
		if ts!=-99.:
			self.ts = ts
		if hd!=-99.:
			self.hd = hd		
		fid = open(cmtfil, 'wt')
		fid.write(' PDE %4d %2d %2d %2d %2d %2d.%02d %8.4f %9.4f %5.1f %3.1f %3.1f %s\n' % \
			  (self.Otime[0],self.Otime[1],self.Otime[2],self.Otime[3],self.Otime[4],
			   self.Otime[5],0,self.lat,self.lon,self.dep,self.mag,self.mag,self.title[7:]))
		fid.write('event name:      %s\n' % self.id)
		fid.write('time shift:%12.4f\n'   % self.ts)
		fid.write('half duration:%9.4f\n' % self.hd)
		fid.write('latitude:%14.4f\n'     % self.lat)
		fid.write('longitude:%13.4f\n'    % self.lon)
		fid.write('depth:%17.4f\n'        % self.dep)
		fid.close()
	def rcmtfile(self,cmtfil):
		fid = open(cmtfil)
		tmp = fid.readline().strip('\n').split()
		dat = '%d/%0d/%0d  --  %0d:%0d:%0d' % \
			  (int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]),int(float(tmp[6])))
		self.Otime = time.strptime(dat,'%Y/%m/%d  --  %H:%M:%S')
		self.lat   = float(tmp[7])
		self.lon   = float(tmp[8])
		self.dep   = float(tmp[9])
		self.mag   = float(tmp[10])
		self.title = '%s '*len(tmp[12:])%tuple(tmp[12:])
		self.title ='M %3.1f, %s'%(self.mag,self.title)
		self.id    = fid.readline().strip('\n').split()[2]
		self.ts    = float(fid.readline().strip('\n').split()[2])
		self.hd    = float(fid.readline().strip('\n').split()[2])
		fid.close
		
	def wimaster(self,DATADIR,DMIN,DMAX,ftable,cmtfile,i_master,gf_dir='./GF',wpwin=[15.],DATALESS=''):
		#for i in xrange(len(ftable[0])):
		#	if self.mag >= ftable[0][i] and self.mag < ftable[1][i]:
		#		break;
		fid = open(i_master, 'wt')
		fid.write('EVNAME:     %s\n'% self.title[7:].replace(' ','_').replace(',',''))
		fid.write('SEED:       %s\n' %  DATADIR)
		if len(DATALESS)>0:
			fid.write('DATALESS:   %s\n' %  DATALESS)
		fid.write('DMIN:       %5.2f\n' % DMIN)
		fid.write('DMAX:       %5.2f\n' % DMAX)
		fid.write('CMTFILE:    %s\n\n' % cmtfile)
		fid.write('filt_order: 4\n')
		fid.write('filt_cf1:   %7.5f\n' % ftable[0])
		fid.write('filt_cf2:   %7.5f\n' % ftable[1])
		fid.write('filt_pass:  1\n')
		fid.write('IDEC_2:  2  200  0.1\n')
		fid.write('IDEC_3:  0.001  0.1  100  0.03\n\n')
		fid.write('GFDIR:   %s\n' % gf_dir)
		fid.write('WP_WIN:  ') 
		for p in wpwin:
			fid.write('%6.1f' % p) 
		fid.write('\n') 
		fid.close()
		
		

def EQcopy(EQ2,EQ1):
	"Deepcopy EarthQuake objects"
	EQ2.id    = EQ1.id
	EQ2.title = EQ1.title
	EQ2.Otime = EQ1.Otime
	EQ2.date  = EQ1.date
	EQ2.lat   = EQ1.lat
	EQ2.lon   = EQ1.lon
	EQ2.dep   = EQ1.dep
	EQ2.mag   = EQ1.mag
	EQ2.org   = EQ1.org
	
def warn_msg(EQ,eq):
	"Display a redundancy warning message"
	print '#############################################'
	print 'Warning : Event redundancy with different id'
	print '  --> %s : ' % eq.id 
	eq.affiche()		
	print '  --> %s (rejected) : ' % EQ.id
	EQ.affiche()
	print '#############################################\n'

def warn_msg2(field,id,val):
	    if val == -99.:
		    addmsg = '("%s" is not a float in rss field)' % field
	    else:
		    addmsg = '(negative depth)'
	    print '#############################################'
	    print 'Warning : Invalid %s' % field
	    print '  -->ev-id : %s , dep = %f %s' % (id,val,addmsg)
	    print '#############################################'	    

def isdiff(EQ1,EQ2):
	if EQ2.id   == EQ1.id and EQ2.title == EQ1.title and EQ2.Otime == EQ1.Otime \
	   and EQ2.date == EQ1.date and EQ2.lat  == EQ1.lat and EQ2.lon  == EQ1.lon \
	   and EQ2.dep  == EQ1.dep and EQ2.mag  == EQ1.mag and EQ2.org  == EQ1.org:
		return 0
	else:
		return 1
    
def screening(EQs,EQ,mintime,minmag,O_scr2,Ep_scr2,M_scr2,count,flag):
	"Screen events"
	# Check focal depth
	if EQ.dep < 0:
	    warn_msg2('depth',EQ.id,EQ.dep)
	    return 0
	# Event id screening
	if EQs.has_key(EQ.id):
	    if flag == 2 and isdiff(EQs[EQ.id],EQ) and EQ.org == 'neic':
	       return 1
            else:
	       return 0
    
	# Date screening
	Otime = calendar.timegm(EQ.Otime)
	delta = Otime - calendar.timegm(mintime)	
	if delta < 0:          
	    return 0

	# Magnitude screening
	if EQ.mag < minmag:
	    return 0

        # Otime and Epicenter screening
	if count > 0 and flag != 0:
	    for key,eq in EQs.items():
	        Otmp = calendar.timegm(eq.Otime)
		if (Otime - Otmp)**2<=O_scr2 and (EQ.lat-eq.lat)**2<=Ep_scr2\
		       and (EQ.lon-eq.lon)**2<=Ep_scr2 and (EQ.mag-eq.mag)**2<=M_scr2:
			datedb = calendar.timegm(eq.date)
			datenw = calendar.timegm(EQ.date)
			if datedb < datenw and EQ.org == 'neic' and flag == 2:
				print 'priority : neic'
				warn_msg(eq,EQ)
				del EQs[key]
				return 1
			warn_msg(EQ,eq)
			return 0
	return 1


def sort_EQ(EQs,sig):
	"Sort by increasing(if sig=1) or decreasing (if sig=-1) Otime"
	eqs = EQs.values()
	N = len(eqs)
	for i in xrange(N-1):
	   date1 = calendar.timegm(eqs[i].Otime)
	   for j in xrange(i+1,N):
	      date2 = calendar.timegm(eqs[j].Otime)
	      if sig*date2 < sig*date1:
		 tmp    = eqs[i]
		 eqs[i] = eqs[j]
		 eqs[j] = tmp
		 date1 = date2
	return eqs
	

def read_db(db_fil,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag=0,EQs={}):
	"Read events database and return a dictionary of EarthQuake objects"
	print db_fil
	EQ  = EarthQuake()  # EarthQuake instantiation	
	count = 0
	if os.access(db_fil,os.R_OK):
 		fid = open(db_fil,'rt')
 		N   = int(fid.readline().strip('\n'))
 		for j in xrange(N):
		   EQ = EarthQuake()
		   EQ.read(fid)
		   if screening(EQs,EQ,mintime,minmag,O_scr2,Ep_scr2,M_scr2,count,flag):
		      EQs[EQ.id] = EarthQuake()
		      EQcopy(EQs[EQ.id],EQ)
		      count += 1 
 		fid.close()
	return EQs

def emsc_feeds2float(line,patern,beg=0):
	try:
	   out = float(re.search(patern,line).group()[-beg:])
	   return out
	except ValueError:
	   return -99.
		
def r_emsc_feeds(url,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag=0,EQs={}):
	"Read emsc feeds and return a dictionary of EarthQuake objects"
	begfeed = re.compile('<item>')  
	endfeed = re.compile('</item>')
	EQ  = EarthQuake()
	# read all feeds
	f     = urllib2.urlopen(url)
	out = ''
	for line in f:
		out = out+line.strip('\n')
	f.close()		
	lines = out.replace('</webMaster>','</webMaster>\n').replace('</item>','</item>\n').split('\n')
	count = len(EQs)		
	for line in lines:
		if begfeed.search(line) and endfeed.search(line):
			EQ.title   = re.search('(?<=<title>).+(?=</title>)',line).group()						
			EQ.id      = re.search('(?<=id).+(?=</link>)',line).group()\
				     .replace('=','')
			EQ.org = 'emsc'

			date = re.search('(?<=<pubDate>).+(?=</pubDate>)',line).group()
			sec  = int(float(date[-4:])+0.5)
			date = '%s -- %s%d' % (date[:10],date[22:-4],sec)

			EQ.mag = emsc_feeds2float(line,'(?<=<mag>).+(?=</mag>)',4)
			EQ.Otime   = time.strptime(date,'%Y-%m-%d -- %H:%M:%S')
			
			# Date screening
			delta = calendar.timegm(EQ.Otime) - calendar.timegm(mintime)	
			if delta < 0:          
				break			
			# Magnitude screening
			if EQ.mag < minmag:
				if EQ.mag == -99.:
					warn_msg2('magnitude',EQ.id,EQ.mag)
				continue

			EQ.lat = emsc_feeds2float(line,'(?<=<geo\:lat>).+(?=</geo\:lat>)')
			EQ.lon = emsc_feeds2float(line,'(?<=<geo\:long>).+(?=</geo\:long>)')
			EQ.dep = emsc_feeds2float(line,'(?<=<depth>).+(?=</depth>)')			
			
			EQ.date    = time.gmtime()
			
			
			if screening(EQs,EQ,mintime,minmag,O_scr2,Ep_scr2,M_scr2,count,flag):
			   EQs[EQ.id] = EarthQuake()
			   EQcopy(EQs[EQ.id],EQ)   		  
			   count += 1
	return EQs


def r_neic_feeds(url,mintime,minmag,O_scr2,Ep_scr2,M_scr2,flag=0,EQs={}):
	"Read neic feeds and return a dictionary of EarthQuake objects"
	count = len(EQs)
	begfeed = re.compile('<entry>')  
	endfeed = re.compile('</entry>') 
	EQ  = EarthQuake()  
	f     = urllib2.urlopen(url)
	# read all feeds	
	for lines in f:
	        line = lines.strip('\n')
		if begfeed.search(line) and endfeed.search(line):
			EQ.id    = re.search('(?<=<id>).+(?=</id>)',line).group()\
				   .replace('urn:earthquake-usgs-gov:','').replace(':','')
			
			EQ.title = re.search('(?<=<title>).+(?=</title>)',line).group()
			EQ.org = 'neic'
			EQ.mag   = float(EQ.title[1:5])			
			date     = re.search('(?<=<updated>).+(?=</updated>)',line).group()
			EQ.Otime = time.strptime(date,'%Y-%m-%dT%H:%M:%SZ')

			EQ.lat,EQ.lon = map(float,re.search('(?<=<georss\:point>).+(?=</georss\:point>)',line).group().split(' '))
			EQ.dep   = float(re.search('(?<=\<georss\:elev\>).+(?=\</georss\:elev\>)',line).group())/-1000
			EQ.date  = time.gmtime()
			
			if screening(EQs,EQ,mintime,minmag,O_scr2,Ep_scr2,M_scr2,count,flag):
			   EQs[EQ.id] = EarthQuake()
			   EQcopy(EQs[EQ.id],EQ)			  
			   count += 1
	f.close()
	return EQs


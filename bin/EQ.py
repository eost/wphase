#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

import os, sys, urllib2
import re
import time, calendar

class EarthQuake:
	"Store Earthquake information"
 	def __init__(self):
		self.id      = '--'
 		self.title   = '--'
		self.pdeline = '--'
		self.org     = '--'
		self.pdebeg  = ' PDE '
		self.Otime   = time.gmtime(0)
		self.date    = time.gmtime(0)
		self.csec    = 0
		self.mag     = 0.0
		self.pdelat  = 0.0 
		self.pdelon  = 0.0
		self.pdedep  = 0.0		
		self.lat     = 0.0 
		self.lon     = 0.0
		self.dep     = 0.0
		self.ts      = 0.0
		self.hd      = 0.0   
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
	def warnmsg(self,fil = 'wallmsg'):
		self.affiche(fil)
		os.system("wall < %s" % fil)
	def pop_pdeline(self):
		self.pdeline = '%4s%4d %2d %2d %2d %2d %2d.%02d %8.4f %9.4f %5.1f %3.1f %3.1f %s' % \
		    (self.pdebeg,self.Otime[0],self.Otime[1],self.Otime[2],self.Otime[3],self.Otime[4],
		     self.Otime[5],self.csec,self.pdelat,self.pdelon,self.pdedep,self.mag,self.mag,self.title)
	def read(self,fid):
		L = fid.readline().strip().split('\t')
		self.id     = L[0].strip()
		self.org    = L[1].strip()
		self.mag    = float(L[2])
		self.pdelat = float(L[3])
		self.pdelon = float(L[4])
		self.pdedep = float(L[5])
		datstr,cstr = L[6].split('.')
		self.Otime  = time.strptime(datstr,'%Y-%m-%d %H:%M:%S')
		self.date   = time.strptime(L[7].split('.')[0],'%Y-%m-%d %H:%M:%S')
		self.csec   = int(cstr[:2])
		if len(L) > 8:
			self.title  = L[8].strip()
		self.lat    = self.pdelat
		self.lon    = self.pdelon
		self.dep    = self.pdedep
		a=self.id[0]
		if a.isalpha():
		    self.org    = 'neic'
		    self.pdebeg = ' NEIC'
	        else:
		    self.org = 'emsc'
		    self.pdebeg = ' EMSC'
		self.pop_pdeline()

	def write(self,fid):
		fid.write('%-15s\t' % self.id)
		fid.write('%-6s\t' % self.org)
		fid.write('%3.1f\t' % self.mag)		
		fid.write('%7.3f\t' % self.lat)
		fid.write('%8.3f\t' % self.lon)
		fid.write('%7.2f\t' % self.dep)
		fid.write('%s.%02d\t' % (time.strftime('%Y-%m-%d %H:%M:%S',self.Otime),self.csec))
		fid.write('%s.00\t' % time.strftime('%Y-%m-%d %H:%M:%S',self.date))
		fid.write('%s\n' % self.title.strip())
	def wcmtfile(self,cmtfil,ts=-99.,hd=-99.):
		if ts!=-99.:
			self.ts = ts
		if hd!=-99.:
			self.hd = hd		
		fid = open(cmtfil, 'wt')
		fid.write('%s\n'%self.pdeline)
		fid.write('event name:      %s\n' % self.id)
		fid.write('time shift:%12.4f\n'   % self.ts)
		fid.write('half duration:%9.4f\n' % self.hd)
		fid.write('latitude:%14.4f\n'     % self.lat)
		fid.write('longitude:%13.4f\n'    % self.lon)
		fid.write('depth:%17.4f\n'        % self.dep)
		fid.close()
	def rcmtfile(self,cmtfil):
		fid          = open(cmtfil)
		self.pdeline = fid.readline().strip('\n')
		try:
			self.mag   = float(self.pdeline[57:60])
		except:
			self.mag   = float(self.pdeline[53:56])
		if self.mag <= 2.:
			self.mag   = float(self.pdeline[53:56])
			if self.mag <= 2.:
				print '**** Warning : preliminary magnitude is very small'
		self.title = self.pdeline[57:].strip().replace(' ','_')
		self.id    = fid.readline().strip('\n').split()[2]
		self.ts    = float(fid.readline().strip('\n').split()[2])
		self.hd    = float(fid.readline().strip('\n').split()[2])
		self.lat   = float(fid.readline().strip('\n').split()[1])
		self.lon   = float(fid.readline().strip('\n').split()[1])
		self.dep   = float(fid.readline().strip('\n').split()[1])
		fid.close()
	def wimaster(self,DATADIR,ftable,cmtfile,i_master,DMIN=0.0,DMAX=90.0,gf_dir='./GF',wpwin=[15.],DATALESS=''):
		fid = open(i_master, 'wt')
		fid.write('EVNAME:     %s\n'% self.title.strip().replace(' ','_').replace(',',''))
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
		fid.write('IDEC_2:  2 1000  0.1\n')
		fid.write('IDEC_3:  0.001  0.1  100  0.01\n\n')
		fid.write('GFDIR:   %s\n' % gf_dir)
		fid.write('WP_WIN:  ') 
		for p in wpwin:
			fid.write('%6.1f' % p) 
		fid.write('\n') 
		fid.close()


def EQcopy(EQ2,EQ1):
	"Deepcopy EarthQuake objects"
	EQ2.id      = EQ1.id
	EQ2.title   = EQ1.title
	EQ2.pdeline = EQ1.pdeline
	EQ2.Otime   = EQ1.Otime
	EQ2.csec    = EQ1.csec
	EQ2.date    = EQ1.date
	EQ2.pdelat  = EQ1.pdelat
	EQ2.pdelon  = EQ1.pdelon
	EQ2.pdedep  = EQ1.pdedep
	EQ2.lat     = EQ1.lat
	EQ2.lon     = EQ1.lon
	EQ2.dep     = EQ1.dep
	EQ2.mag     = EQ1.mag
	EQ2.org     = EQ1.org

def isdiff(EQ1,EQ2):
	if EQ2.id   == EQ1.id and EQ2.title == EQ1.title and EQ2.Otime == EQ1.Otime \
	   and EQ2.date == EQ1.date and EQ2.lat  == EQ1.lat and EQ2.lon  == EQ1.lon \
	   and EQ2.dep  == EQ1.dep and EQ2.mag  == EQ1.mag and EQ2.org  == EQ1.org  \
	   and EQ2.pdelat  == EQ1.pdelat and EQ2.pdelon == EQ1.pdelon and EQ2.pdedep  == EQ1.pdedep\
	   and EQ2.csec == EQ1.csec and EQ2.pdeline == EQ1.pdeline :
		return 0
	else:
		return 1

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


def screening(EQs,EQ,mintime,minmag,O_scr2,Ep_scr2,M_scr2,count,flag):
	"Screen events"
	# Check focal depth
	if EQ.pdedep < 0:
	    warn_msg2('depth',EQ.id,EQ.pdedep)
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
		if (Otime - Otmp)**2<=O_scr2 and (EQ.pdelat-eq.pdelat)**2<=Ep_scr2\
		       and (EQ.pdelon-eq.pdelon)**2<=Ep_scr2 and (EQ.mag-eq.mag)**2<=M_scr2:
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
		L   = fid.readlines()
		N = 0
		for l in L:
			if len(l) == 0:
				break
			N += 1
		fid.seek(0,0)
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

			EQ.pdelat = emsc_feeds2float(line,'(?<=<geo\:lat>).+(?=</geo\:lat>)')
			EQ.pdelon = emsc_feeds2float(line,'(?<=<geo\:long>).+(?=</geo\:long>)')
			EQ.pdedep = emsc_feeds2float(line,'(?<=<depth>).+(?=</depth>)')			
			EQ.lat    = EQ.pdelat
			EQ.lon    = EQ.pdelon
			EQ.dep    = EQ.pdedep
			EQ.date   = time.gmtime()
			EQ.pop_pdeline()
			
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
				   .replace('urn:earthquake-usgs-gov:','').replace(':','').strip()
			
			EQ.title = re.search('(?<=<title>).+(?=</title>)',line).group().strip()
			EQ.org = 'neic'
			EQ.mag   = float(EQ.title[1:5])			
			date     = re.search('(?<=<updated>).+(?=</updated>)',line).group()
			EQ.Otime = time.strptime(date,'%Y-%m-%dT%H:%M:%SZ')
			EQ.pdelat,EQ.pdelon = map(float,re.search('(?<=<georss\:point>).+(?=</georss\:point>)',line).group().split(' '))
			EQ.pdedep   = float(re.search('(?<=\<georss\:elev\>).+(?=\</georss\:elev\>)',line).group())/-1000
			EQ.lat      = EQ.pdelat
			EQ.lon      = EQ.pdelon
			EQ.dep      = EQ.pdedep
			EQ.date     = time.gmtime()
			EQ.pop_pdeline()
			if screening(EQs,EQ,mintime,minmag,O_scr2,Ep_scr2,M_scr2,count,flag):
			   if EQ.title[0] == 'M' or EQ.title[0] == 'm':
				   EQ.title = EQ.title[7:]
			   EQs[EQ.id] = EarthQuake()
			   EQcopy(EQs[EQ.id],EQ)			  
			   count += 1
	f.close()
	return EQs

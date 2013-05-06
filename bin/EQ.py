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

class EarthQuake:
	"Store Earthquake information"
 	def __init__(self):
		self.id      = '--'
 		self.title   = '--'
		self.pdeline = '--'
		self.mag     = 0.0		
		self.lat     = 0.0 
		self.lon     = 0.0
		self.dep     = 0.0
		self.ts      = 0.0
		self.hd      = 0.0   
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
		self.title = self.pdeline[61:].strip().replace(' ','_')
		self.id    = fid.readline().strip('\n').split()[2]
		self.ts    = float(fid.readline().strip('\n').split()[2])
		self.hd    = float(fid.readline().strip('\n').split()[2])
		self.lat   = float(fid.readline().strip('\n').split()[1])
		self.lon   = float(fid.readline().strip('\n').split()[1])
		self.dep   = float(fid.readline().strip('\n').split()[1])
		fid.close()
	def wimaster(self,DATADIR,ftable,cmtfile,i_master,DMIN=0.0,DMAX=90.0,gf_dir='./GF',IDEC2=[2,290,0.1],wpwin=[15.],DATALESS=''):
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
		fid.write('IDEC_2:  %d %3d  %.1f\n'%tuple(IDEC2))
		fid.write('IDEC_3:  0.001  0.1  100  0.03\n\n')
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
	EQ2.mag     = EQ1.mag
	EQ2.lat     = EQ1.lat
	EQ2.lon     = EQ1.lon
	EQ2.dep     = EQ1.dep
	EQ2.ts      = EQ1.ts
	EQ2.hd      = EQ1.hd

def isdiff(EQ1,EQ2):
	if EQ2.id   == EQ1.id and EQ2.title == EQ1.title and EQ2.pdeline == EQ1.pdeline \
	   and EQ2.lat  == EQ1.lat and EQ2.lon  == EQ1.lon and EQ2.dep  == EQ1.dep  \
	   and EQ2.mag  == EQ1.mag and EQ2.ts == EQ1.ts and EQ2.hd == EQ1.hd:
		return 0
	else:
		return 1
